# coupeler_driver

# import packages
using Pkg
import SciMLBase: step!
using OrdinaryDiffEq
using OrdinaryDiffEq: ODEProblem, solve, SSPRK33, savevalues!, Euler
using LinearAlgebra
import Test: @test
import ClimaCore.Spaces as Spaces

# import coupler utils
include("coupler_utils/flux_calculator.jl")
include("coupler_utils/conservation_checker.jl")
include("coupler_utils/regridder.jl")

# initiate spatial and temporal info
t_end = 2592000 # 100e2
tspan = (0, t_end)
Δt_cpl = 2e2
saveat = Δt_cpl * 1000

# init MPI
include("mpi/mpi_init.jl")

# init model components
include("atmos/atmos_init.jl")
atmos_sim = atmos_init(FT, Y, spaces, integrator, params = params);

include("slab/slab_init.jl")
boundary_space = ClimaCore.Fields.level(atmos_sim.domain.center_space, 1)
slab_sim = slab_init(FT, tspan, dt = Δt_cpl, space = boundary_space, saveat = saveat);

# init boundary fields for regridding (TODO: technically this can be bypassed by directly rigridding on model grids)
T_S = ClimaCore.Fields.zeros(boundary_space)
F_S = ClimaCore.Fields.zeros(boundary_space)

# init conservation info collector
CS = ConservationCheck([], [])

# coupling loop
@show "Starting coupling loop"
walltime = @elapsed for t in (tspan[1]:Δt_cpl:tspan[end])
    @show t

    ## Atmos
    # calculate surface fluxes for Atmos BCs
    T_S .= FT(0)
    dummmy_remap!(T_S, slab_sim.integrator.u.T_sfc)

    #atmos_sim.integrator.p.ρ_dif_flux_h_tot .= ClimaCore.Geometry.WVector.(ClimaCore.Fields.zeros(axes(atmos_sim.integrator.p.ρ_dif_flux_h_tot)))
    parent(atmos_sim.integrator.p.ρ_dif_flux_h_tot) .= FT(0)
    calculate_surface_fluxes_atmos_grid!(atmos_sim.integrator, T_S)

    # run
    step!(atmos_sim.integrator, t - atmos_sim.integrator.t, true) # NOTE: use (t - integ_atm.t) here instead of Δt_cpl to avoid accumulating roundoff error in our timestepping.

    ## Slab
    # pre: get accumulated flux from atmos
    F_S .= ClimaCore.Fields.zeros(boundary_space)
    dummmy_remap!(F_S, atmos_sim.integrator.p.ρ_dif_flux_h_tot)

    # save the accumulated flux
    slab_F_sfc = slab_sim.integrator.p.F_sfc
    @. slab_F_sfc = -F_S

    # run
    step!(slab_sim.integrator, t - slab_sim.integrator.t, true)

    # conservation info "callback"
    if !is_distributed
        check_conservation_callback(CS, atmos_sim, slab_sim)
    end

end
@show walltime

@show "Postprocessing"
# collect solutions
sol_atm = atmos_sim.integrator.sol
sol_slab = slab_sim.integrator.sol

include("mpi/mpi_postprocess.jl")

# conservation  check
if !is_distributed || (is_distributed && ClimaComms.iamroot(comms_ctx))
    conservation_plot(atmos_sim, slab_sim, solu_atm, solu_slab, "conservation.png")
end
