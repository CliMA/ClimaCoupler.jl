using Oceananigans.TimeSteppers: time_step!, Clock, tick!

import SciMLBase: step!

using Printf

const FT = Float64
include("SurfaceFluxes.jl")
include("dummy_surface_fluxes.jl") 
include("land_simulation.jl") 
include("atmos_simulation.jl")

using CLIMAParameters
using CLIMAParameters.Planet: cp_d, cv_d, grav, T_surf_ref
using CLIMAParameters.Atmos.SubgridScale: C_smag, C_drag
struct EarthParameterSet <: AbstractEarthParameterSet end
const CLIMAparam_set = EarthParameterSet()
import CLIMAParameters


struct CoupledSimulation{A, L, C}
    atmos :: A
    land :: L
    clock :: C
end


# abstract type BC{FT <: AbstractFloat} end
# mutable struct FluxBC{FT} <: BC{FT}
#     bottom_heat_flux::FT
#     top_heat_flux::FT
# end

function step!(coupled_sim::CoupledSimulation, coupling_Δt)

    atmos_sim = coupled_sim.atmos
    land_sim = coupled_sim.land
    clock = coupled_sim.clock
    next_time = clock.time + coupling_Δt
    @info("Coupling cycle", time = clock.time)

    # Coupler _pull: Extract states and parameters at coupled boundaries for flux calculations
    land_surface_T = land_sim.p[5]

    # Coupler_push: Update atmos variables for this coupling_Δt 
    atmos_sim.u.x[3] .= [0.0, 0.0, 0.0] # reset surface flux to be accumulated during each coupling_Δt 
    atmos_sim.p[2] .= land_surface_T # update the T_sfc aux variable

    # Advance atmos
    step!(atmos_sim, next_time - atmos_sim.t, true)

    # Coupler_pull: Extract surface fluxes from atmos
    ∫surface_x_momentum_flux = atmos_sim.u.x[3][1] # kg / m s^2
    ∫surface_y_momentum_flux = atmos_sim.u.x[3][2] # kg / m s^2
    ∫surface_heat_flux = atmos_sim.u.x[3][3]       # W / m^2
        
    # Coupler_push: Update land variables for this coupling_Δt 
    land_sim.p[4].top_heat_flux = ∫surface_heat_flux / coupling_Δt # the F_accum aux variable 

    # Advance land
    step!(land_sim, next_time - land_sim.t, true)

    tick!(clock, coupling_Δt)
    
    return nothing
end

function solve!(coupled_sim::CoupledSimulation, coupling_Δt, stop_time)
    @info("Coupler:", models = fieldnames(typeof(coupled_sim))[1:end-1])
    while coupled_sim.clock.time < stop_time
        step!(coupled_sim, coupling_Δt)
    end
end

# f = 1e-4 # Coriolis parameter
# g = 9.81 # Gravitational acceleration

atmos_Nz = 30  # Number of vertical grid points
atmos_Lz = 200 # Vertical extent of domain

# land domain params (need to abstract)
# n = 50
# zmax = FT(0)
# zmin = FT(-1)
land_Nz = n
land_Lz = zmax - zmin # this needs abstraction

start_time = 0.0 
stop_time = 1#coupling_Δt*100#60*60

# Build the respective models
land_sim = land_simulation()
atmos_sim = atmos_simulation(land_sim, Nz=atmos_Nz, Lz=atmos_Lz, start_time = start_time, stop_time = stop_time,)

# Build a coupled simulation
clock = Clock(time=0.0)
coupled_sim = CoupledSimulation(atmos_sim, land_sim, clock)

# Run it!
coupling_Δt = 0.02
solve!(coupled_sim, coupling_Δt, stop_time)

# Visualize
using Plots
Plots.GRBackend()

dirname = "heat"
path = joinpath(@__DIR__, "output", dirname)
mkpath(path)

# Atmos plots
sol_atm = coupled_sim.atmos.sol
t0_ρθ = parent(sol_atm.u[1].x[1])[:,4]
tend_ρθ = parent(sol_atm.u[end].x[1])[:,4]
t0_u = parent(sol_atm.u[1].x[1])[:,2]
tend_u = parent(sol_atm.u[end].x[1])[:,2]
t0_v = parent(sol_atm.u[1].x[1])[:,3]
tend_v = parent(sol_atm.u[end].x[1])[:,3]
z_centers =  collect(1:1:length(tend_u))#parent(Fields.coordinate_field(center_space_atm))[:,1]
Plots.png(Plots.plot([t0_ρθ tend_ρθ],z_centers, labels = ["t=0" "t=end"]), joinpath(path, "T_atm_height.png"))
Plots.png(Plots.plot([t0_u tend_u],z_centers, labels = ["t=0" "t=end"]), joinpath(path, "u_atm_height.png"))
Plots.png(Plots.plot([t0_v tend_v],z_centers, labels = ["t=0" "t=end"]), joinpath(path, "v_atm_height.png"))

# Land plots
sol_lnd = coupled_sim.land.sol
t0_θ_l = parent(sol_lnd.u[1].x[1])
tend_θ_l = parent(sol_lnd.u[end].x[1])
t0_ρe = parent(sol_lnd.u[1].x[3])
tend_ρe = parent(sol_lnd.u[end].x[3])
t0_ρc_s = volumetric_heat_capacity.(t0_θ_l, parent(θ_i), Ref(msp.ρc_ds), Ref(param_set)) #convert energy to temp
t0_T = temperature_from_ρe_int.(t0_ρe, parent(θ_i),t0_ρc_s, Ref(param_set)) #convert energy to temp
tend_ρc_s = volumetric_heat_capacity.(tend_θ_l, parent(θ_i), Ref(msp.ρc_ds), Ref(param_set)) #convert energy to temp
tend_T = temperature_from_ρe_int.(tend_ρe, parent(θ_i),tend_ρc_s, Ref(param_set)) #convert energy to temp
z_centers =  collect(1:1:length(tend_ρe))#parent(Fields.coordinate_field(center_space_atm))[:,1]
Plots.png(Plots.plot([t0_θ_l tend_θ_l],parent(zc), labels = ["t=0" "t=end"]), joinpath(path, "Th_l_lnd_height.png"))
Plots.png(Plots.plot([t0_T tend_T],parent(zc), labels = ["t=0" "t=end"]), joinpath(path, "T(K)_lnd_height.png"))

# Time evolution of all enthalpies (do for total enery = p / (γ - 1) + dot(ρu, ρu) / 2ρ + ρ * Φ, which is actually conserved)
atm_sum_u_t = [sum(parent(u.x[1])[:,4]) for u in sol_atm.u] ./ atmos_Nz .* atmos_Lz * parameters.C_p # J / m2
lnd_sfc_u_t = [sum(parent(u.x[3])[:]) for u in sol_lnd.u] ./ land_Nz .* land_Lz # J / m2

v1 = lnd_sfc_u_t .- lnd_sfc_u_t[1]
v2 = atm_sum_u_t .- atm_sum_u_t[1]
Plots.png(Plots.plot(sol_lnd.t, [v1 v2 v1+v2], labels = ["lnd" "atm" "tot"]), joinpath(path, "enthalpy_both_surface_time_atm_lnd.png"))

# relative error of enthapy
using Statistics
total = atm_sum_u_t + lnd_sfc_u_t
rel_error = (total .- total[1]) / mean(total)
Plots.png(Plots.plot(sol_lnd.t, rel_error, labels = ["tot"]), joinpath(path, "rel_error_surface_time.png"))

# TODO
# -