# coupler_driver

# import packages
using Pkg
import SciMLBase: step!
using OrdinaryDiffEq
using OrdinaryDiffEq: ODEProblem, solve, SSPRK33, savevalues!, Euler
using LinearAlgebra
import Test: @test
using ClimaCore.Utilities: half, PlusHalf
Pkg.add(PackageSpec(name = "ClimaCore", version = "0.10.3"))

# import coupler utils
include("coupler_utils/flux_calculator.jl")
include("coupler_utils/conservation_checker.jl")
include("coupler_utils/regridder.jl")
include("coupler_utils/masker.jl")
include("coupler_utils/general_helper.jl")

# # initiate spatial and temporal info
t_end = 2592000 * 3 # 100e2 # 2592000 # 100e2 #2592000 * 2 #
tspan = (0, t_end) # 172800.0)
Δt_cpl = 2e2
saveat = Δt_cpl * 100

# init MPI
include("mpi/mpi_init.jl")

# init atmos model component
include("atmos/atmos_init.jl")
atmos_sim = atmos_init(FT, Y, spaces, integrator, params = params);

# init land-sea mask
infile = "data/seamask.nc"
boundary_space = ClimaCore.Fields.level(atmos_sim.domain.face_space, half) # global surface grid
mask = LandSeaMask(FT, infile, "LSMASK", boundary_space)

# init surface (slab) model components
include("slab/slab_init.jl")
slab_sim = slab_init(FT, tspan, dt = Δt_cpl, space = boundary_space, saveat = saveat, mask = mask);

include("slab_ocean/slab_init.jl")
slab_ocean_sim = slab_ocean_init(FT, tspan, dt = Δt_cpl, space = boundary_space, saveat = saveat, mask = mask);

# init coupler's boundary fields for regridding (TODO: technically this can be bypassed by directly rigridding on model grids)
T_S = ClimaCore.Fields.zeros(boundary_space) # temperature
z0m_S = ClimaCore.Fields.zeros(boundary_space)
z0b_S = ClimaCore.Fields.zeros(boundary_space)

F_A = ClimaCore.Fields.zeros(boundary_space) # aerodynamic turbulent fluxes
F_R = ClimaCore.Fields.zeros(boundary_space) # radiative fluxes

# init conservation info collector
CS = ConservationCheck([], [])

# sample SST field
# SST = ncreader_rll_to_cgll_from_space(FT, "data/sst.nc",  "SST", boundary_space)    
# SST = swap_space!(SST,axes(mask)) .* (  abs.(mask .-1) ) .+ FT(273.15)

# coupling loop
@show "Starting coupling loop"
walltime = @elapsed for t in (tspan[1]:Δt_cpl:tspan[end])
    @show t

    ## Atmos
    ## Turbulent surface fluxes

    # coupler_get: T_sfc, z0m, z0b
    combined_field = zeros(boundary_space)
    parent(combined_field) .=
        combine_surface.(parent(mask), parent(slab_sim.integrator.u.T_sfc), parent(slab_ocean_sim.integrator.u.T_sfc))
    # parent(combined_field) .= combine_surface.(parent(mask), parent(slab_sim.integrator.u.T_sfc), parent(SST) ) # prescribed SSTs
    dummmy_remap!(T_S, combined_field)

    parent(combined_field) .=
        combine_surface.(
            parent(mask),
            parent(slab_sim.integrator.p.params.z0m .* mask),
            parent(slab_ocean_sim.integrator.p.params.z0m .* (abs.(mask .- 1))),
        )
    dummmy_remap!(z0m_S, combined_field)
    parent(combined_field) .=
        combine_surface.(
            parent(mask),
            parent(slab_sim.integrator.p.params.z0b .* mask),
            parent(slab_ocean_sim.integrator.p.params.z0b .* (abs.(mask .- 1))),
        )
    dummmy_remap!(z0b_S, combined_field)

    # calculate turbulent fluxes on atmos grid and save in atmos cache
    info_sfc = (; T_sfc = T_S, z0m = z0m_S, z0b = z0b_S)
    calculate_surface_fluxes_atmos_grid!(atmos_sim.integrator, info_sfc)

    # run 
    step!(atmos_sim.integrator, t - atmos_sim.integrator.t, true) # NOTE: instead of Δt_cpl, to avoid accumulating roundoff error

    #clip TODO: this is bad!! > limiters
    parent(atmos_sim.integrator.u.c.ρq_tot) .= heaviside.(parent(atmos_sim.integrator.u.c.ρq_tot)) # negligible for total energy cons

    # coupler_push!: get accumulated fluxes from atmos in the surface fields
    F_A .= ClimaCore.Fields.zeros(boundary_space)
    dummmy_remap!(F_A, atmos_sim.integrator.p.dif_flux_energy)
    F_R .= ClimaCore.Fields.zeros(boundary_space)
    parsed_args["rad"] == "gray" ? dummmy_remap!(F_R, level(atmos_sim.integrator.p.ᶠradiation_flux, half)) : nothing # TODO: albedo hard coded...

    ## Slab land
    # coupler_get: F_aero, F_rad
    slab_F_aero = slab_sim.integrator.p.F_aero
    @. slab_F_aero = -F_A
    slab_F_rad = slab_sim.integrator.p.F_rad
    @. slab_F_rad = -F_R

    # run
    step!(slab_sim.integrator, t - slab_sim.integrator.t, true)

    ## Slab ocean
    # coupler_get: F_aero, F_rad
    slab_ocean_F_aero = slab_ocean_sim.integrator.p.F_aero
    @. slab_ocean_F_aero = -F_A
    slab_ocean_F_rad = slab_ocean_sim.integrator.p.F_rad
    @. slab_ocean_F_rad = -F_R

    # run
    step!(slab_ocean_sim.integrator, t - slab_ocean_sim.integrator.t, true)

    # conservation info "callback" logging at every Δt_cpl
    if !is_distributed
        check_conservation_callback(CS, atmos_sim, slab_sim)
    end

end

@show walltime

@show "Postprocessing"
# collect solutions
sol_atm = atmos_sim.integrator.sol
sol_slab = slab_sim.integrator.sol
sol_slab_ocean = slab_ocean_sim.integrator.sol

include("mpi/mpi_postprocess.jl")

# conservation  check
if !is_distributed || (is_distributed && ClimaComms.iamroot(comms_ctx))
    conservation_plot(atmos_sim, slab_sim, solu_atm, solu_slab, "conservation.png")
end

# animations
plot_anim = nothing
if plot_anim !== nothing

    using ClimaCorePlots

    anim = Plots.@animate for u in sol_atm.u
        Plots.plot(Fields.level(Geometry.UVVector.(u.c.uₕ).components.data.:1, 1))
    end
    Plots.mp4(anim, "anim_u.mp4", fps = 10)

    anim = Plots.@animate for u in sol_atm.u
        Plots.plot(Fields.level(Geometry.UVVector.(u.c.uₕ).components.data.:1, 5))
    end
    Plots.mp4(anim, "anim_u_7km.mp4", fps = 10)

    anim = Plots.@animate for u in sol_atm.u
        Plots.plot(Fields.level(u.c.ρe, 1))
    end
    Plots.mp4(anim, "anim_rhoe.mp4", fps = 10)

    anim = Plots.@animate for u in sol_slab.u
        Plots.plot(u.T_sfc, clims = (240, 330))
    end
    Plots.mp4(anim, "slab_T.mp4", fps = 10)

    anim = Plots.@animate for u in sol_atm.u
        Plots.plot(Fields.level(u.c.ρe, 1) .- Fields.level(sol_atm.u[1].c.ρe, 1), clims = (-5000, 50000))
        println(parent(Fields.level(u.c.ρe, 1) .- Fields.level(sol_atm.u[1].c.ρe, 1))[1])
    end
    Plots.mp4(anim, "anim_rhoe_anom.mp4", fps = 10)

    anim = Plots.@animate for u in sol_atm.u
        Plots.plot(Fields.level(u.c.ρe, 7) .- Fields.level(sol_atm.u[1].c.ρe, 7), clims = (-1000, 3000))
    end
    Plots.mp4(anim, "anim_rhoe_anom_7km.mp4", fps = 10)

    anim = Plots.@animate for u in sol_atm.u
        Plots.plot(Fields.level(u.c.ρq_tot, 1))#.- Fields.level(sol_atm.u[1].c.ρt_tot,1),  clims = (-5000, 50000) )
    end
    Plots.mp4(anim, "anim_rhoqt.mp4", fps = 10)

    anim = Plots.@animate for u in sol_atm.u
        Plots.plot(Fields.level(u.c.ρq_tot, 2), clims = (0, 0.02))#.- Fields.level(sol_atm.u[1].c.ρt_tot,1),  clims = (-5000, 50000) )
    end
    Plots.mp4(anim, "anim_rhoqt_1km.mp4", fps = 10)

    # anim = Plots.@animate for u in sol_atm.u
    #     Plots.plot(Fields.level(Geometry.WVector.(u.f.w),half) )#.- Fields.level(sol_atm.u[1].c.ρt_tot,1),  clims = (-5000, 50000) )
    # end
    # Plots.mp4(anim, "anim_w.mp4", fps = 10)
end

# TODO:
# - update MPI, conservation plots 