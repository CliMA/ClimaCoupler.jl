# coupler_driver
# don't forget to run with threading: julia --project --threads 8 (MPI not that useful for debugging coarse runs)

import SciMLBase: step!
using OrdinaryDiffEq
using OrdinaryDiffEq: ODEProblem, solve, SSPRK33, savevalues!, Euler
using LinearAlgebra
import Test: @test
using ClimaCore.Utilities: half, PlusHalf
using Dates

include("cli_options.jl")
(s, parsed_args) = parse_commandline()
# Read in some parsed args
prescribed_sst = false
energy_check = true
const FT = parsed_args["FLOAT_TYPE"] == "Float64" ? Float64 : Float32
land_sim = "bucket"
t_end = FT(86400*20)#time_to_seconds(parsed_args["t_end"]))
tspan = (0, t_end)
saveat = time_to_seconds(parsed_args["dt_save_to_sol"])
saveat = 3600
date0 = date = DateTime(1979, 01, 01)
date1 = Dates.firstdayofmonth(date0) # first date
# overwrite some parsed args :P
parsed_args["t_end"] = "20days"
parsed_args["coupled"] = true
Δt_cpl=FT(200)
parsed_args["dt"] = string(Δt_cpl) * "secs"
parsed_args["enable_threading"] = true
parsed_args["microphy"] = "0M"
parsed_args["forcing"] = nothing
parsed_args["idealized_h2o"] = false
parsed_args["vert_diff"] = true
parsed_args["rad"] = "gray"
parsed_args["hyperdiff"] = true
parsed_args["config"] = "sphere"
parsed_args["moist"] = "equil"
energy_check=true
prescribed_sst=false
anim = true
import ClimaCoupler
pkg_dir = pkgdir(ClimaCoupler)
coupler_output_dir = joinpath(pkg_dir, "experiments/AMIP/moist_mpi_earth")

# Get the paths to the necessary data files - land sea mask, sst map, sea ice concentration
include("artifacts.jl")

# import coupler utils
include("coupler_utils/flux_calculator.jl")
include("coupler_utils/conservation_checker.jl")
include("coupler_utils/regridder.jl")
include("coupler_utils/masker.jl")
include("coupler_utils/calendar_timer.jl")
include("coupler_utils/general_helper.jl")

# init MPI
include("mpi/mpi_init.jl")

# init atmos model component
include("atmos/atmos_init.jl")
atmos_sim = atmos_init(FT, Y, integrator, params = params);

# init a 2D bounary space at the surface, assuming the same instance (and MPI distribution if applicable) as the atmos domain above
boundary_space = ClimaCore.Fields.level(atmos_sim.domain.face_space, half) # global surface grid

# init land-sea mask

mask = LandSeaMask(FT, mask_data, "LSMASK", boundary_space)
mask = swap_space!(mask, boundary_space) # needed if we are reading from previous run

# init surface (slab) model components
# we need some types that are defined in these files
include("slab/slab_utils.jl")
include("bucket/bucket_init.jl")
include("slab/slab_init.jl")
include("slab_ocean/slab_init.jl")
include("slab_ice/slab_init.jl")

if land_sim == "slab"
    slab_sim = slab_init(FT; tspan, dt = Δt_cpl, space = boundary_space, saveat = saveat, mask = mask)
end
if land_sim == "bucket"
    slab_sim = bucket_init(FT, FT.(tspan); dt = FT(Δt_cpl), space = boundary_space, saveat = FT(saveat))
end

if prescribed_sst
    println("No ocean sim - do not expect energy conservation")
    weightfile, datafile_cgll, regrid_space =
        ncreader_rll_to_cgll_from_space(sst_data, "SST", boundary_space, outfile = "sst_cgll.nc")
    SST = ncreader_cgll_sparse_to_field(datafile_cgll, "SST", weightfile, (Int(1),), regrid_space)[1]
    SST = swap_space!(SST, axes(mask)) .* (abs.(mask .- 1)) .+ FT(273.15)

    ocean_params = OceanSlabParameters(FT(20), FT(1500.0), FT(800.0), FT(280.0), FT(1e-3), FT(1e-5), FT(0.06))
    slab_ocean_sim = nothing
else
    slab_ocean_sim =
        slab_ocean_init(FT; tspan = tspan, dt = Δt_cpl, space = boundary_space, saveat = saveat, mask = mask)
    SST = nothing
    ocean_params = nothing
end

# Currently, we only support a slab ice model with fixed area and depth.
weightfile, datafile_cgll, regrid_space =
    ncreader_rll_to_cgll_from_space(sic_data, "SEAICE", boundary_space, outfile = "sic_cgll.nc")
SIC = ncreader_cgll_sparse_to_field(datafile_cgll, "SEAICE", weightfile, (Int(1),), regrid_space)[1]
SIC = swap_space!(SIC, axes(mask)) .* (abs.(mask .- 1))
ice_mask = get_ice_mask.(SIC .- FT(25), FT) # here 25% and lower is considered ice free

slab_ice_sim =
    slab_ice_init(FT; tspan = tspan, dt = Δt_cpl, space = boundary_space, saveat = saveat, ice_mask = ice_mask)


# init coupler
coupler_sim = CouplerSimulation(FT(Δt_cpl), integrator.t, boundary_space, FT, mask)
include("./push_pull.jl")

# init conservation info collector
atmos_pull!(
    atmos_sim,
    slab_ice_sim,
    slab_sim,
    slab_ocean_sim,
    boundary_space,
    prescribed_sst,
    z0m_S,
    z0b_S,
    T_S,
    ρ_sfc,
    q_sfc,
    ocean_params,
    SST,
    mask,
)

atmos_push!(atmos_sim, boundary_space, F_A, F_R, F_E, P_liq, P_snow,parsed_args)
land_pull!(slab_sim, F_A, F_R, F_E, P_liq, P_snow,ρ_sfc)
ice_pull!(slab_ice_sim, F_A, F_R)
if !prescribed_sst
    ocean_pull!(slab_ocean_sim, F_A, F_R)
    reinit!(slab_ocean_sim.integrator)
end

reinit!(atmos_sim.integrator)
reinit!(slab_sim.integrator)
reinit!(slab_ice_sim.integrator)
if !is_distributed && energy_check && !prescribed_sst
    CS = OnlineConservationCheck([], [], [], [], [], [])
    check_conservation(CS, coupler_sim, atmos_sim, slab_sim, slab_ocean_sim, slab_ice_sim)
end

# coupling loop
@show "Starting coupling loop"

walltime = @elapsed for t in ((tspan[1] + Δt_cpl):Δt_cpl:tspan[end])

    date = current_date(t)

    @calendar_callback :(@show(date), date1 += Dates.Month(1)) date date1

    ## Atmos
    atmos_pull!(
        atmos_sim,
        slab_ice_sim,
        slab_sim,
        slab_ocean_sim,
        boundary_space,
        prescribed_sst,
        z0m_S,
        z0b_S,
        T_S,
        ρ_sfc,
        q_sfc,
        ocean_params,
        SST,
        mask,
    )
    step!(atmos_sim.integrator, t - atmos_sim.integrator.t, true) # NOTE: instead of Δt_cpl, to avoid accumulating roundoff error

    atmos_push!(atmos_sim, boundary_space, F_A, F_R, F_E, P_liq, P_snow, parsed_args)

    ## Slab land
    land_pull!(slab_sim, F_A, F_R, F_E, P_liq, P_snow, ρ_sfc)
    step!(slab_sim.integrator, t - slab_sim.integrator.t, true)

    ## Slab ocean
    if !prescribed_sst
        ocean_pull!(slab_ocean_sim, F_A, F_R)
        step!(slab_ocean_sim.integrator, t - slab_ocean_sim.integrator.t, true)
    end

    ## Slab ice
    ice_pull!(slab_ice_sim, F_A, F_R)
    step!(slab_ice_sim.integrator, t - slab_ice_sim.integrator.t, true)

    ## Compute energy
    if !is_distributed && energy_check && !prescribed_sst
        check_conservation(CS, coupler_sim, atmos_sim, slab_sim, slab_ocean_sim, slab_ice_sim)
    end
end

@show walltime

@show "Postprocessing"
if energy_check && !prescribed_sst
    plot_global_energy(CS, coupler_sim, "total_energy_bucket_snow.png", "total_energy_log_bucket_snow.png")
end

# # animations
if (land_sim == "bucket") && anim
    #make it so this works with slab land?
    include("coupler_utils/viz_explorer.jl")
    plot_anim(atmos_sim, slab_sim, slab_ocean_sim, slab_ice_sim, mask, prescribed_sst, SST)
end

# Cleanup temporary files
# TODO: Where should this live?
rm(REGRID_DIR; recursive = true, force = true)
