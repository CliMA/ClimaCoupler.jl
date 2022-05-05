# coupeler_driver

# import packages
using Pkg
import SciMLBase: step!
using OrdinaryDiffEq
using OrdinaryDiffEq: ODEProblem, solve, SSPRK33, savevalues!, Euler
using LinearAlgebra
import Test: @test

# import coupler utils
include("coupler_utils/flux_calculator.jl")
include("coupler_utils/conservation_checker.jl")
include("coupler_utils/regridder.jl")

# initiate spatial and temporal info
t_end = 2592000 # 100e2 
tspan = (0, t_end) # 172800.0)
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

    #atmos_sim.integrator.p.dif_flux_energy .= ClimaCore.Geometry.WVector.(ClimaCore.Fields.zeros(axes(atmos_sim.integrator.p.dif_flux_energy)))
    parent(atmos_sim.integrator.p.dif_flux_energy) .= FT(0)
    calculate_surface_fluxes_atmos_grid!(atmos_sim.integrator, T_S)

    # run 
    step!(atmos_sim.integrator, t - atmos_sim.integrator.t, true) # NOTE: use (t - integ_atm.t) here instead of Δt_cpl to avoid accumulating roundoff error in our timestepping.

    ## Slab
    # pre: get accumulated flux from atmos
    F_S .= ClimaCore.Fields.zeros(boundary_space)
    dummmy_remap!(F_S, atmos_sim.integrator.p.dif_flux_energy)

    # save the accumulated flux
    slab_F_sfc = slab_sim.integrator.p.F_sfc
    @. slab_F_sfc = -F_S # Δt_cpl /  Δt_cpl

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

#=
# Strong scaling
# helem = 4, zelem = 10, nq = 5, timeend = 1e4 s, dt = 1e2 s
no_processors = [1,2,4,8]
walltime_list = [160, 118, 119, 118]
#walltime_list_weak = [170, 118, ]
no_processors_weak = [1,4,16]
walltime_list_weak = [170, 171, 170]
Plots.scatter(no_processors, walltime_list, label="strong_coupled")
Plots.scatter!(no_processors_weak, walltime_list_weak, label="weak_coupled", ylims = (110,190))
Plots.savefig("scaling.png")

#MacOS
#- atmos_omly: 90.98 s
#- coupled: 90.71 s
#- atmos_only_solve! 77.54 s

=#
