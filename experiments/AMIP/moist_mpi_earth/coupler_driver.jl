# coupeler_driver

# import packages
using Pkg
import SciMLBase: step!
using OrdinaryDiffEq
using OrdinaryDiffEq: ODEProblem, solve, SSPRK33, savevalues!, Euler
using LinearAlgebra
import Test: @test
using ClimaCore.Utilities: half, PlusHalf

# import coupler utils
include("coupler_utils/flux_calculator.jl")
include("coupler_utils/conservation_checker.jl")
include("coupler_utils/regridder.jl")
include("coupler_utils/ncfile_reader.jl")

# # initiate spatial and temporal info
t_end = 2592000 # 2592000 # 100e2 
tspan = (0, t_end) # 172800.0)
Δt_cpl = 2e2
saveat = Δt_cpl * 100

# init MPI
include("mpi/mpi_init.jl")

# init model components
include("atmos/atmos_init.jl")
atmos_sim = atmos_init(FT, Y, spaces, integrator, params = params);

include("slab/slab_init.jl")
boundary_space = ClimaCore.Fields.level(atmos_sim.domain.face_space, half)
slab_sim = slab_init(FT, tspan, dt = Δt_cpl, space = boundary_space, saveat = saveat, mask = mask);

include("slab_cold/slab_init.jl")
slab_cold_sim = slab_cold_init(FT, tspan, dt = Δt_cpl, space = boundary_space, saveat = saveat, mask = mask);

# init boundary fields for regridding (TODO: technically this can be bypassed by directly rigridding on model grids)
T_S = ClimaCore.Fields.zeros(boundary_space)
F_S = ClimaCore.Fields.zeros(boundary_space)
F_R = ClimaCore.Fields.zeros(boundary_space)

# init land-sea mask
infile = "data/seamask.nc"
R = atmos_sim.domain.face_space.horizontal_space.topology.mesh.domain.radius
h_elem = atmos_sim.domain.face_space.horizontal_space.topology.mesh.ne
Nq = Spaces.Quadratures.polynomial_degree(atmos_sim.domain.face_space.horizontal_space.quadrature_style) + 1
mask = LandSeaMask(FT, infile, "LSMASK", ne = h_elem, R = R, Nq = Nq)

# init conservation info collector
CS = ConservationCheck([], [])

combine_surface(mask, sfc_1, sfc_2) = (mask > FT(0.5) ?  sfc_1 : FT(0)) + (mask < 0.5 ?  sfc_2 : FT(0)) 

# coupling loop
@show "Starting coupling loop"
walltime = @elapsed for t in (tspan[1]:Δt_cpl:tspan[end])
    @show t

    ## Atmos
    # calculate surface fluxes for Atmos BCs
    T_S .= FT(0)

    T_sfc_all = zeros( axes(slab_sim.integrator.u.T_sfc) )    
    parent(T_sfc_all) .= combine_surface.(parent(mask), parent(slab_sim.integrator.u.T_sfc), parent(slab_cold_sim.integrator.u.T_sfc)) 
    
    dummmy_remap!(T_S, T_sfc_all)

    #atmos_sim.integrator.p.dif_flux_energy .= ClimaCore.Geometry.WVector.(ClimaCore.Fields.zeros(axes(atmos_sim.integrator.p.dif_flux_energy)))
    parent(atmos_sim.integrator.p.dif_flux_energy) .= FT(0)
    calculate_surface_fluxes_atmos_grid!(atmos_sim.integrator, T_S)

    # run 
    step!(atmos_sim.integrator, t - atmos_sim.integrator.t, true) # NOTE: use (t - integ_atm.t) here instead of Δt_cpl to avoid accumulating roundoff error in our timestepping.

    ## Slab
    # pre: get accumulated flux from atmos
    F_S .= ClimaCore.Fields.zeros(boundary_space)
    dummmy_remap!(F_S, atmos_sim.integrator.p.dif_flux_energy)
    F_R .= ClimaCore.Fields.zeros(boundary_space)
    dummmy_remap!(F_R, level(atmos_sim.integrator.p.ᶠradiation_flux, half)) # TODO: albedo hard coded...

    # save the accumulated flux
    slab_F_sfc = slab_sim.integrator.p.F_sfc
    @. slab_F_sfc = - F_S # Δt_cpl /  Δt_cpl
    slab_F_rad = slab_sim.integrator.p.F_rad
    @. slab_F_rad = - F_R # Δt_cpl /  Δt_cpl

    # run
    step!(slab_sim.integrator, t - slab_sim.integrator.t, true)

    ## Slab_cold
    # pre: get accumulated flux from atmos

    # save the accumulated flux
    slab_cold_F_sfc = slab_cold_sim.integrator.p.F_sfc 
    @. slab_cold_F_sfc = - F_S # Δt_cpl /  Δt_cpl
    slab_cold_F_rad = slab_cold_sim.integrator.p.F_rad
    @. slab_cold_F_rad = - F_R # Δt_cpl /  Δt_cpl

    # run
    step!(slab_cold_sim.integrator, t - slab_cold_sim.integrator.t, true)

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

# animations
using ClimaCorePlots

anim = Plots.@animate for u in sol_atm.u
    Plots.plot(Fields.level(Geometry.UVVector.(u.c.uₕ).components.data.:1,1))
end
Plots.mp4(anim, "anim_u.mp4", fps = 10)

anim = Plots.@animate for u in sol_atm.u
    Plots.plot(Fields.level(u.c.ρe,1))
end
Plots.mp4(anim, "anim_rhoe.mp4", fps = 10)

anim = Plots.@animate for u in sol_slab.u
    Plots.plot(u.T_sfc)
end
Plots.mp4(anim, "slab_T.mp4", fps = 10)

anim = Plots.@animate for u in sol_atm.u
    Plots.plot(Fields.level(u.c.ρe,1) .- Fields.level(sol_atm.u[1].c.ρe,1),  clims = (-1000, 7000) )
    println(parent(Fields.level(u.c.ρe,1) .- Fields.level(sol_atm.u[1].c.ρe,1))[1])
end
Plots.mp4(anim, "anim_rhoe_anom.mp4", fps = 10)