using CLIMAParameters
using CLIMAParameters.Planet: ρ_cloud_liq, ρ_cloud_ice, cp_l, cp_i, T_0, LH_f0
struct EarthParameterSet <: AbstractEarthParameterSet end
const param_set = EarthParameterSet()

using RecursiveArrayTools
using OrdinaryDiffEq: ODEProblem, solve, SSPRK33,Rosenbrock23, Tsit5,SSPRK432, Feagin14, TsitPap8,CarpenterKennedy2N54
using DifferentialEquations
using UnPack
using LandHydrology
using LandHydrology.SoilInterface
using LandHydrology.Domains: Column
using LandHydrology.SoilInterface.SoilHeatParameterizations
using LandHydrology.SoilInterface.SoilWaterParameterizations

const FT = Float64

# General composition
ν = FT(0.5);
ν_ss_quartz = FT(0.92)
ν_ss_minerals = FT(0.08)
ν_ss_om = FT(0.0)
ν_ss_gravel = FT(0.0);

#Water specific
Ksat = FT(0.0443 / 3600 / 100) # m/s
S_s = FT(1e-3) #inverse meters
vg_n = FT(2.0)
vg_α = FT(2.6); # inverse meters
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
msp = SoilParams{FT}(
        ν,
        S_s,
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
        κ_dry_parameter,
    )

#Simulation and domain info
t0 = FT(0)
tf = FT(60 * 60 * 72)
dt = FT(0.02)

n = 20

zmax = FT(0)
zmin = FT(-2)
domain =  Column(FT, zlim = (zmin, zmax), nelements = n)

# Boundary conditions
top_water_flux = FT(0)
top_heat_flux = FT(0)
bottom_water_flux = FT(0)
bottom_heat_flux = FT(0)
bc = SoilDomainBC(
    domain;
    top = SoilComponentBC(
        hydrology = VerticalFlux(top_water_flux),
        energy = VerticalFlux(top_heat_flux),
    ),
    bottom = SoilComponentBC(
        hydrology = VerticalFlux(bottom_water_flux),
        energy = VerticalFlux(bottom_heat_flux),
    ),
)

# creat model
hydraulics_model = vanGenuchten{FT}(n = vg_n, α = vg_α, Ksat = Ksat, θr = θ_r)

soil_model = SoilModel(
        domain = domain,
        energy_model = SoilEnergyModel(),
        hydrology_model = SoilHydrologyModel{FT}(
            hydraulic_model = hydraulics_model,
        ),
        boundary_conditions = bc,
        soil_param_set = msp,
        earth_param_set = param_set,
    )

# initial conditions
function initial_conditions(z::FT, t0::FT, model::SoilModel)
    param_set = model.earth_param_set
    T = 289.0 + 5.0 * z
    θ_i = 0.0
    θ_l = 0.495
    ρcds = model.soil_param_set.ρc_ds
    ρc_s = volumetric_heat_capacity(θ_l, θ_i, ρcds, param_set)
    ρe_int = volumetric_internal_energy(θ_i, ρc_s, T, param_set)
    return (ϑ_l = θ_l, θ_i = θ_i, ρe_int = ρe_int)
end
Y = set_initial_state(soil_model, initial_conditions, 0.0)

soil_rhs! = make_rhs(soil_model)
land_prob = ODEProblem(soil_rhs!, Y, (t0, tf), [])
algorithm = CarpenterKennedy2N54()

land_simulation() = init(land_prob, algorithm, dt = dt, saveat = 1 * dt) # dt is the land model step

function coupler_push!(coupler::CouplerState, sim)
    # only need to do this calculation at the boundary.
    T = temperature_from_ρe_int.(Y.ρe_int, Y.θ_i, Y.ρc_s, Ref(param_set))
    T_sfc = parent(T)[end]
    coupler_put!(coupler, :T_sfc, T_sfc, nothing, DateTime(0), u"K")

end

function coupler_pull!(sim, coupler::CouplerState)
    ∫surface_heat_flux = atmos.u.x[3][3]
    land.p[4].top_heat_flux = ∫surface_heat_flux / coupling_Δt # [W/m^2] same BC across land Δt
end