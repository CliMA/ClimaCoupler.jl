const FT = Float64
import SciMLBase: step!
using OrdinaryDiffEq
using OrdinaryDiffEq: ODEProblem, solve, SSPRK33, savevalues!, Euler
using LinearAlgebra
using ClimaCore

# pending PR merge:
# import Pkg; Pkg.add(url="https://github.com/CliMA/ClimaCore.jl",rev="sb/online-sphere-remap", subdir = "lib/ClimaCoreTempestRemap")

# using old interface - before the new one is consolidated:
import Pkg;
Pkg.add(url = "https://github.com/CliMA/ClimaAtmos.jl", rev = "jh_ln/external_boundary_conditons");
using ClimaCoreTempestRemap

import Test: @test
import ClimaAtmos.BoundaryConditions: get_boundary_flux, AbstractBoundary
export get_boundary_flux

include("atmos.jl")
include("slab.jl")
include("coupled_bc.jl")
include("conservation.jl")

# initiate spatial and temporal info
tspan = (0.0, 1300.0) # 172800.0)
cpl_Δt = 5.0 # seconds
atmos_dt, slab_dt = (5.0, 5.0) # seconds

slab_sim = slab_init(FT, tspan, dt = atmos_dt, npolynomial = 3)
atmos_sim = atmos_init(FT, tspan, dt = slab_dt, npolynomial = 3)

# parameter exchange
slab_sim.integrator.p.cp_d .= CLIMAParameters.Planet.cp_d(atmos_sim.model.parameters)

# regridding init
atmos_energy_bc = atmos_sim.model.boundary_conditions.ρe_tot.bottom
atmos_momentum_bc = atmos_sim.model.boundary_conditions.uh.bottom

F_a_space = axes(atmos_energy_bc.flux)
F_s_space = axes(slab_sim.integrator.u.T_sfc)

# weightfile = tempname()
# R_atm2slab = ClimaCoreTempestRemap.generate_map( # bring back once CC #614 is resolved/ merged 
#     F_s_space, #target
#     F_a_space, #source
#     weightfile = weightfile,
# )
# R_slab2atm = ClimaCoreTempestRemap.generate_map(# bring back once CC #614 is resolved/ merged 
#     F_a_space, #target
#     F_s_space, #source
#     weightfile = weightfile,
# )
function dummmy_remap!(target, source)  # delete when can use Tempest again
    parent(target) .= parent(source)
end

# init conservation info collector
CS = ConservationCheck([], [])

# coupling loop
for t in (tspan[1]:cpl_Δt:tspan[end])
    @show t

    ## Atmos
    # calculate surface fluxes for Atmos BCs
    coupler_atmos_boundary_flux(atmos_energy_bc, atmos_sim, slab_sim, F_a_space, F_s_space)
    coupler_atmos_boundary_flux(atmos_momentum_bc, atmos_sim, slab_sim, F_a_space, F_s_space)

    # run 
    step!(atmos_sim.integrator, t - atmos_sim.integrator.t, true) # NOTE: use (t - integ_atm.t) here instead of Δt_cpl to avoid accumulating roundoff error in our timestepping.

    ## Slab
    # pre: get accumulated flux from atmos
    F_S = ClimaCore.Fields.zeros(F_s_space)

    # ClimaCoreTempestRemap.remap!(F_S, atmoatmos_energy_bcs_bc.flux, R_atm2slab)
    dummmy_remap!(F_S, atmos_energy_bc.flux)

    # save the accumulated flux
    slab_F_sfc = slab_sim.integrator.p.F_sfc
    slab_F_sfc .= F_S .* cpl_Δt

    # run
    step!(slab_sim.integrator, t - slab_sim.integrator.t, true)

    # conservation info
    check_conservation(CS, atmos_sim, slab_sim)
end
