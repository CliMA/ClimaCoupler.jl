# coupler_driver
# don't forget to run with threading: julia --project --threads 8 (MPI not that useful for debugging coarse runs)

# import packages

using Pkg
import SciMLBase: step!
using OrdinaryDiffEq
using OrdinaryDiffEq: ODEProblem, solve, SSPRK33, savevalues!, Euler
using LinearAlgebra
import Test: @test
using ClimaCore.Utilities: half, PlusHalf
Pkg.add(PackageSpec(name = "ClimaCore", version = "0.10.3"))
## When is FT set? It is already being used (as a global) in flux_calculator
# import coupler utils
include("coupler_utils/flux_calculator.jl")
include("coupler_utils/conservation_checker.jl")
include("coupler_utils/regridder.jl")
include("coupler_utils/masker.jl")
include("coupler_utils/general_helper.jl")

# # initiate spatial and temporal info
debug_mode = true
t_end = debug_mode ? 100e2 : 2592000 * 3
Δt_cpl = 2e2
saveat = debug_mode ? Δt_cpl * 1 : Δt_cpl * 100

tspan = (0, t_end)

# init MPI
include("mpi/mpi_init.jl")

# init atmos model component
include("atmos/atmos_init.jl")
atmos_sim = atmos_init(FT, Y, spaces, integrator, params = params);


@show debug_mode

# init a 2D bounary space at the surface, assuming the same instance (and MPI distribution if applicable) as the atmos domain above
boundary_space = ClimaCore.Fields.level(atmos_sim.domain.face_space, half) # global surface grid

# init land-sea mask
infile = "data/seamask.nc"
mask = LandSeaMask(FT, infile, "LSMASK", boundary_space) # TODO: split up the nc file to individual times for faster computation
mask .= FT(1.0)
# init surface (slab) model components
# ClimaLSM unregistered:
#Pkg.add( url = "https://github.com/CliMA/ClimaLSM.jl", rev = "move_coupled_types")
Pkg.develop(path="../../../../ClimaLSM.jl")

include("bucket/bucket_init.jl")
bucket_sim = bucket_init(FT, FT.(tspan); dt = FT(Δt_cpl), space = boundary_space, saveat = FT(saveat));

include("slab/slab_utils.jl")

include("slab_ocean/slab_init.jl")
prescribed_sst = false
if prescribed_sst == true
    SST = ncreader_rll_to_cgll_from_space(FT, "data/sst.nc", "SST", boundary_space)  # a sample SST field from https://gdex.ucar.edu/dataset/158_asphilli.html
    SST = swap_space!(SST, axes(mask)) .* (abs.(mask .- 1)) .+ FT(273.15) # TODO: avoids the "space not the same instance" error
    ocean_params = OceanSlabParameters(FT(20), FT(1500.0), FT(800.0), FT(280.0), FT(1e-3), FT(1e-5))
else
    slab_ocean_sim = slab_ocean_init(FT, tspan, dt = Δt_cpl, space = boundary_space, saveat = saveat, mask = mask)
end

include("slab_ice/slab_init.jl")
prescribed_sic = false
if prescribed_sic == true
    # sample SST field
    SIC = ncreader_rll_to_cgll_from_space(FT, "data/sic.nc", "SEAICE", boundary_space)
    SIC = swap_space!(SIC, axes(mask)) .* (abs.(mask .- 1))
    slab_ice_sim = slab_ice_init(FT, tspan, dt = Δt_cpl, space = boundary_space, saveat = saveat, prescribed_sic = SIC)
else
    slab_ice_sim = slab_ice_init(FT, tspan, dt = Δt_cpl, space = boundary_space, saveat = saveat)
end

# init coupler's boundary fields for regridding (TODO: technically this can be bypassed by directly rigridding on model grids)
T_S = ClimaCore.Fields.zeros(boundary_space) # temperature
z0m_S = ClimaCore.Fields.zeros(boundary_space)
z0b_S = ClimaCore.Fields.zeros(boundary_space)

F_A = ClimaCore.Fields.zeros(boundary_space) # aerodynamic turbulent fluxes
F_E = ClimaCore.Fields.zeros(boundary_space) # evaporation due to turbulent fluxes
F_R = ClimaCore.Fields.zeros(boundary_space) # radiative fluxes
dF_A = ClimaCore.Fields.zeros(boundary_space) # aerodynamic turbulent fluxes

# init conservation info collector
CS = ConservationCheck([], [],[])

# coupling loop
@show "Starting coupling loop"
walltime = @elapsed for t in (tspan[1]:Δt_cpl:tspan[end])
    @show t

    ## Atmos
    ## Turbulent surface fluxes

    # coupler_get: T_sfc, z_0m, z_0b
    combined_field = zeros(boundary_space)
    if prescribed_sst == true
        parent(combined_field) .=
            combine_surface.(
                parent(mask) .- parent(slab_ice_sim.integrator.p.ice_mask .* FT(2)),
                parent(bucket_sim.integrator.u.bucket.T_sfc),
                parent(SST),
                parent(slab_ice_sim.integrator.u.T_sfc),
            ) # prescribed SSTs
        dummmy_remap!(T_S, combined_field)
        parent(combined_field) .=
            combine_surface.(
                parent(mask),
                parent(bucket_sim.params.z_0m .* mask),
                parent(ocean_params.z0m .* (abs.(mask .- 1))),
            )
        dummmy_remap!(z0m_S, combined_field)
        parent(combined_field) .=
            combine_surface.(
                parent(mask),
                parent(bucket_sim.params.z_0b .* mask),
                parent(ocean_params.z0b .* (abs.(mask .- 1))),
            )
        dummmy_remap!(z0b_S, combined_field)
    else
        parent(combined_field) .=
            combine_surface.(
                parent(mask) .- parent(slab_ice_sim.integrator.p.ice_mask .* FT(2)),
                parent(bucket_sim.integrator.u.bucket.T_sfc),
                parent(slab_ocean_sim.integrator.u.T_sfc),
                parent(slab_ice_sim.integrator.u.T_sfc),
            )
        dummmy_remap!(T_S, combined_field)
        parent(combined_field) .=
            combine_surface.(
                parent(mask),
                parent(bucket_sim.params.z_0m .* mask),
                parent(slab_ocean_sim.integrator.p.params.z0m .* (abs.(mask .- 1))),
            )
        dummmy_remap!(z0m_S, combined_field)
        parent(combined_field) .=
            combine_surface.(
                parent(mask),
                parent(bucket_sim.params.z_0b .* mask),
                parent(slab_ocean_sim.integrator.p.params.z0b .* (abs.(mask .- 1))),
            )
        dummmy_remap!(z0b_S, combined_field)
    end

    # calculate turbulent fluxes on atmos grid and save in atmos cache
    info_sfc = (; T_sfc = T_S, z0m = z0m_S, z0b = z0b_S, ice_mask = slab_ice_sim.integrator.p.ice_mask)
    # This should also need q_sfc - currently it assumes q_sfc = q_sat
    calculate_surface_fluxes_atmos_grid!(atmos_sim.integrator, info_sfc)

    # run 
    step!(atmos_sim.integrator, t - atmos_sim.integrator.t, true) # NOTE: instead of Δt_cpl, to avoid accumulating roundoff error

    #clip TODO: this is bad!! > limiters
    parent(atmos_sim.integrator.u.c.ρq_tot) .= heaviside.(parent(atmos_sim.integrator.u.c.ρq_tot)) # negligible for total energy cons

    # coupler_push!: get accumulated fluxes from atmos in the surface fields
    F_A .= ClimaCore.Fields.zeros(boundary_space)
    dummmy_remap!(F_A, atmos_sim.integrator.p.dif_flux_energy)
    F_E .= ClimaCore.Fields.zeros(boundary_space)
    dummmy_remap!(F_E, atmos_sim.integrator.p.dif_flux_ρq_tot)
    F_R .= ClimaCore.Fields.zeros(boundary_space)
    parsed_args["rad"] == "gray" ? dummmy_remap!(F_R, level(atmos_sim.integrator.p.ᶠradiation_flux, half)) : nothing # TODO: albedo hard coded...
    dF_A .= ClimaCore.Fields.zeros(boundary_space)
    dummmy_remap!(dF_A, atmos_sim.integrator.p.∂F_aero∂T_sfc)

    ## Bucket Land
    # coupler_get: F_aero, F_rad
    # standin for ρsfc computation
    @. bucket_sim.integrator.p.bucket.ρ_sfc = FT(1.1)
    @. bucket_sim.integrator.p.bucket.SHF = F_A
    @. bucket_sim.integrator.p.bucket.LHF = FT(0.0)
    @. bucket_sim.integrator.p.bucket.E = F_E
    @. bucket_sim.integrator.p.bucket.R_n = F_R

    # run
    #step!(bucket_sim.integrator, FT(0.0), true) #???? do we need this each step? why?
    # Why is print t = t+dt rather than t?
    # Why is printed dY.T_sfc NE 0 but bucket_sim.u(t+dt) = bucket_sim.u(t)
    step!(bucket_sim.integrator, t - bucket_sim.integrator.t, true)

    ## Slab ocean
    # coupler_get: F_aero, F_rad
    if (prescribed_sst !== true) && (prescribed_sic == true)
        slab_ocean_F_aero = slab_ocean_sim.integrator.p.F_aero
        @. slab_ocean_F_aero = -F_A
        slab_ocean_F_rad = slab_ocean_sim.integrator.p.F_rad
        @. slab_ocean_F_rad = -F_R

        # run
        step!(slab_ocean_sim.integrator, t - slab_ocean_sim.integrator.t, true)
    end
    # conservation info "callback" logging at every Δt_cpl

    ## Slab ice
    # coupler_get: F_aero, F_rad
    slab_ice_F_aero = slab_ice_sim.integrator.p.F_aero
    @. slab_ice_F_aero = -F_A
    slab_ice_F_rad = slab_ice_sim.integrator.p.F_rad
    @. slab_ice_F_rad = -F_R
    slab_ice_∂F_aero∂T_sfc = slab_ice_sim.integrator.p.∂F_aero∂T_sfc
    @. slab_ice_∂F_aero∂T_sfc = dF_A

    # run
    step!(slab_ice_sim.integrator, t - slab_ice_sim.integrator.t, true)

    if !is_distributed
        check_conservation_callback(CS, atmos_sim, bucket_sim, F_A .+ F_R)
    end

end

@show walltime

@show "Postprocessing"
# collect solutions
sol_atm = atmos_sim.integrator.sol
sol_slab = bucket_sim.integrator.sol
sol_slab_ice = slab_ice_sim.integrator.sol
sol_slab_ocean = prescribed_sst !== true ? slab_ocean_sim.integrator.sol : nothing

include("mpi/mpi_postprocess.jl")

# conservation  check
if !is_distributed || (is_distributed && ClimaComms.iamroot(comms_ctx))
    conservation_plot(CS, "conservation.png")
end

# animations
include("coupler_utils/viz_explorer.jl")
plot_anim()

# TODO:
# - update MPI, conservation plots 
# - performance checks, incl threading (--threads=8)
# - add prescribable albedo to RRTMGP + ClimaAtmos
# - add in oupler specific abstractions
# - replace heavisides with smooth functions
# - test dynamical ea ice model 
# - test sea ice for conservation with slab ocean

# Next PR
# - keep adding coupler specific interface
# - fluxes: re-enable different ways to calculate / accumulate fluxes (at overy coupler timestep; at every atmos timestep via callback; via specification of an additional variable) 
# - formalize conservation/physical/performance tests: add error threshold and exception, interval, show option, and make a general interface for it
# - SurfaceFluxes: combine LHF and SHF into enthalpy flux formulation to avoid division by zero
# - Temporally varying SSTs/sea ice: https://pcmdi.llnl.gov/mips/amip/details/
