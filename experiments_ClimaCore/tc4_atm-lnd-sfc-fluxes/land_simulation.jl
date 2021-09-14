# add https://github.com/CliMA/ClimaCore.jl
# add https://github.com/CliMA/ClimaAtmos.jl
# add https://github.com/CliMA/LandHydrologyjl

import ClimaCore.Geometry, LinearAlgebra, UnPack
import ClimaCore:
    Fields,
    Domains,
    Topologies,
    Meshes,
    DataLayouts,
    Operators,
    Geometry,
    Spaces

using CLIMAParameters
using CLIMAParameters.Planet: ρ_cloud_liq, ρ_cloud_ice, cp_l, cp_i, T_0, LH_f0
struct EarthParameterSet <: AbstractEarthParameterSet end
const param_set = EarthParameterSet()

using RecursiveArrayTools
using OrdinaryDiffEq: ODEProblem, solve, SSPRK33,Rosenbrock23, Tsit5,SSPRK432, Feagin14, TsitPap8,CarpenterKennedy2N54
using DifferentialEquations
using UnPack
using LandHydrology
using LandHydrology.Domains
using LandHydrology.SoilInterface
using LandHydrology.SoilInterface.SoilWaterParameterizations
using LandHydrology.SoilInterface.SoilHeatParameterizations



abstract type BC{FT <: AbstractFloat} end

mutable struct FluxBC{FT} <: BC{FT}
    top_heat_flux::FT
    top_water_flux::FT
    btm_heat_flux::FT
    btm_water_flux::FT
end

function compute_soil_rhs!(dY, Y, t, p)
    
    sp = p[1]
    param_set = p[2]
    zc = p[3]
    @unpack top_heat_flux, btm_heat_flux, top_water_flux, btm_water_flux = p[4]
    @unpack ν,vgn,vgα,vgm,ksat,θr,ρc_ds, κ_sat_unfrozen, κ_sat_frozen = sp
    # ϑ_l = Y.x[1]
    # θ_i = Y.x[2]
    # ρe_int = Y.x[3]

    # dϑ_l = dY.x[1]
    # dθ_i = dY.x[2]
    # dρe_int = dY.x[3]
    dθ_l = dY.θ_l
    dθ_i = dY.θ_i
    dρe_int = dY.ρe_int
    θ_l = Y.θ_l
    θ_i = Y.θ_i
    ρe_int = Y.ρe_int

    # update water content based on prescribed profiles, set RHS to zero.
    # ϑ_l = hydrology.ϑ_l_profile.(zc, t)
    # θ_i = hydrology.θ_i_profile.(zc, t)
    # dθ_i = Fields.zeros(FT, cspace)
    # dρe_int = Fields.zeros(FT, space)

    # Compute center values of everything
    ρc_s = volumetric_heat_capacity.(θ_l, θ_i, ρc_ds, Ref(param_set))
    T = temperature_from_ρe_int.(ρe_int, θ_i, ρc_s, Ref(param_set))
    κ_dry = k_dry(param_set, sp)
    S_r = relative_saturation.(θ_l, θ_i, ν)
    kersten = kersten_number.(θ_i, S_r, Ref(sp))
    κ_sat =
        saturated_thermal_conductivity.(
            θ_l,
            θ_i,
            κ_sat_unfrozen,
            κ_sat_frozen,
        )
    κ = thermal_conductivity.(κ_dry, kersten, κ_sat)

    # rhs operators
    interpc2f = Operators.InterpolateC2F()
    gradc2f_heat = Operators.GradientC2F()
    divf2c_heat = Operators.DivergenceF2C(
        top = Operators.SetValue(Geometry.Cartesian3Vector(zero(FT))),
        bottom = Operators.SetValue(Geometry.Cartesian3Vector(zero(FT))),
    )
    @. dρe_int = -divf2c_heat(-interpc2f(κ) * gradc2f_heat(T))
    return dY
end

# General composition
ν = FT(0.395);
ν_ss_quartz = FT(0.92)
ν_ss_minerals = FT(0.08)
ν_ss_om = FT(0.0)
ν_ss_gravel = FT(0.0);

#Water specific
Ksat = FT(4.42 / 3600 / 100) # m/s
S_s = FT(1e-3) #inverse meters
vg_n = FT(1.89)
vg_α = FT(7.5); # inverse meters
vg_m = FT(1) -FT(1)/vg_n
θ_r = FT(0)

# Heat specific
κ_quartz = FT(7.7) # W/m/K
κ_minerals = FT(2.5) # W/m/K
κ_om = FT(0.25) # W/m/K
κ_liq = FT(0.57) # W/m/K
κ_ice = FT(2.29); # W/m/K
ρp = FT(2700); # kg/m^3
κ_solid = k_solid(ν_ss_om, ν_ss_quartz, κ_quartz, κ_minerals, κ_om)
κ_sat_frozen = ksat_frozen(κ_solid, ν, κ_ice)
κ_sat_unfrozen = ksat_unfrozen(κ_solid, ν, κ_liq);
ρc_ds = FT((1 - ν) * 1.926e06); # J/m^3/K
a = FT(0.24)
b = FT(18.1)
κ_dry_parameter = FT(0.053)

#collect all params
msp = SoilParams{FT}(ν,vg_n,vg_α,vg_m, Ksat, θ_r, S_s,
                     ν_ss_gravel,
                     ν_ss_om,
                     ν_ss_quartz,
                     ρc_ds,
                     κ_solid,
                     ρp,
                     κ_sat_unfrozen,
                     κ_sat_frozen,
                     a,
                     b,
                     κ_dry_parameter)


#Simulation and domain info
t0 = FT(0)
tf = FT(60 * 60 * 72)
dt = FT(0.02)

n = 50

zmax = FT(0)
zmin = FT(-1)
domain = Domains.IntervalDomain(zmin, zmax, x3boundary = (:bottom, :top))
mesh = Meshes.IntervalMesh(domain, nelems = n)

cs = Spaces.CenterFiniteDifferenceSpace(mesh)
fs = Spaces.FaceFiniteDifferenceSpace(cs)
zc = Fields.coordinate_field(cs)

# Boundary conditions
top_water_flux = FT(0)
top_heat_flux = FT(0)
bottom_water_flux = FT(0)
bottom_heat_flux = FT(0)
bc = FluxBC(top_heat_flux,
            top_water_flux,
            bottom_heat_flux,
            bottom_water_flux)

# Parameter structure

# initial conditions
T_max = FT(289.0)
T_min = FT(288.0)
c = FT(20.0)
T = @.  T_min + (T_max - T_min) * exp(-(zc - zmax) / (zmin - zmax) * c)
Tsfc = parent(T)[end]
p = [msp, param_set, zc,bc,Tsfc]
θ_i = Fields.zeros(FT,cs)

theta_max = FT(ν * 0.5)
theta_min = FT(ν * 0.4)
θ_l = @. theta_min + (theta_max - theta_min) * exp(-(zc - zmax) / (zmin - zmax) * c)

ρc_s = volumetric_heat_capacity.(θ_l, θ_i, ρc_ds, Ref(param_set))
ρe_int = volumetric_internal_energy.(θ_i, ρc_s, T, Ref(param_set))

#Y = ArrayPartition(θ_l, θ_i, ρe_int)

Y = Fields.FieldVector(θ_l = θ_l, θ_i = θ_i, ρe_int = ρe_int)

function ∑land_tendencies!(dY, Y, p, t)
    # Intermediate step to be added if needed
    compute_soil_rhs!(dY, Y, t, p)
end

land_prob = ODEProblem(∑land_tendencies!, Y, (t0, tf), p)
algorithm = CarpenterKennedy2N54()

land_simulation() = init(land_prob, algorithm, dt = dt, saveat = 1 * dt) # dt is the land model step

# function  surface_temperature_from_soil(land_sim)
#     θ_l = parent(land_sim.u.x[1])
#     ρe = parent(land_sim.u.x[3])
#     #convert energy to temp
#     ρc_s = volumetric_heat_capacity.(θ_l, parent(θ_i), Ref(msp.ρc_ds), Ref(param_set))
#     tend_T = temperature_from_ρe_int.(ρe, parent(θ_i),ρc_s, Ref(param_set))
# end