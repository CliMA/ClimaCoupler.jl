using Oceananigans.TimeSteppers: time_step!, Clock, tick!

import SciMLBase: step!

using Printf

include("dummy_surface_fluxes.jl") # placeholder for SurfaceFluxes.jl

include("land_simulation.jl") #refactoring of land interface to come
include("atmos_simulation.jl")

struct CoupledSimulation{A, L, C}
    atmos :: A
    land :: L
    clock :: C
end

mutable struct SensibleHeatFlux
    Cd
    T_land
end


#playback = les boundaries driving GCM (x,y,z,t), coupled = atmos accumulating flux, coupled to land
Y_atmos = (Y_gcm, Y_playback)

needed_value = model.prescribed_field(x,y,z,t) # prescribed clouds

needed_value = Y.cloud
Y.cloud (x,y,z,t, RHS = 0) Y.cloud from dY/dt

Y = (Y_energy, Y_water)
Y = (Y_energy, Y_water) 
dY_energy = 0
Y_energy = f(x,y,z,t)
get_temperature(model::EnergyModel (prescribed or not), Y,z, t)

Y_initial = ()
push!(Y_initial, Y_atmos)
push!(Y_inital, Y_coupled)

dY_coupled/dt = flux


∑sensible = ∑_layers a*(T_amos(layer)-T_leaf(layer))
make_rhs!(atmos_model)
    rhs_atmos! = make_rhs(intrinsic_atmos)
    rhs_coupler! = make_rhs(flux_accum_model)
    return rhs!
end

    (Y = rho, u, v, w, Y_coupler = flux_accum)


function calc_flux_in_atmos(Y, Ya, t, bc::BoundaryCondition)
    wind_speed = f(Y)
    T_atm = f(Y)
    T_sfc = bc.T_sfc

    return bc.Cd * wind_speed * (T_atm - T_sfc)
end

# .. meanwhile in the rhs we need
btm_flux = calc_flux_in_atmos(Y,D,t,model.bcs.ρθ.btm)

#= TODO
Interface
- coupler function that overwrites T_sfc in SensibleHeatFlux
- enable appending variables in Y (make_rhs!) for flux accumulation and prescribed/ dynamic flux comm; 

SurfaceFluxes.jl 
- 
- 
- 

=#

function step!(coupled_sim::CoupledSimulation, coupling_Δt)

    atmos_sim = coupled_sim.atmos
    land_sim = coupled_sim.land

    clock = coupled_sim.clock
    next_time = clock.time + coupling_Δt
    @info("Coupling cycle", time = clock.time)

    # Extract states and parameters at coupled boundaries for flux calculations
    land_surface_T = land_sim.p[5]

    # Step forward atmosphere
    atmos_sim.u.x[3] .= [0.0, 0.0, 0.0] # reset surface flux to be accumulated during each coupling_Δt 
    atmos_sim.p[2] .= land_surface_T # get land temperature and set on atmosphere (Tland is prognostic)

    # TODO: determine if this will be useful (init step not ran w/o this but same outcome)
    # u_atmos = atmos_sim.u 
    # u_atmos.x[3] .= u_atmos.x[3] .* -0.0
    # set_u!(atmos_sim, u_atmos)

    step!(atmos_sim, next_time - atmos_sim.t, true)

    # Extract surface fluxes for land boundaries
    ∫surface_x_momentum_flux = atmos_sim.u.x[3][1] # kg / m s^2
    ∫surface_y_momentum_flux = atmos_sim.u.x[3][2] # kg / m s^2
    ∫surface_heat_flux = atmos_sim.u.x[3][3]       # W / m^2
        
    # Advance land
    @show(∫surface_x_momentum_flux, ∫surface_y_momentum_flux, ∫surface_heat_flux)
    land_sim.p[4].top_heat_flux = ∫surface_heat_flux / coupling_Δt # [W/m^2] same BC across land Δt
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

f = 1e-4 # Coriolis parameter
g = 9.81 # Gravitational acceleration

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