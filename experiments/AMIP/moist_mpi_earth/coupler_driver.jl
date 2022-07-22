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
mode_name = "aquaplanet"
energy_check = true
anim = true

energy_check = parsed_args["energy_check"]
const FT = parsed_args["FLOAT_TYPE"] == "Float64" ? Float64 : Float32
land_sim_name = "bucket"
t_end = FT(time_to_seconds(parsed_args["t_end"]))
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

if land_sim_name == "bucket"
    land_sim = bucket_init(FT, FT.(tspan); dt = FT(Δt_cpl), space = boundary_space, saveat = FT(saveat))
end

if mode_name == "amip"
    println("No ocean sim - do not expect energy conservation")

    # ocean
    weightfile, datafile_cgll, regrid_space =
        ncreader_rll_to_cgll_from_space(sst_data, "SST", boundary_space, outfile = "sst_cgll.nc")
    SST_init = ncreader_cgll_sparse_to_field(datafile_cgll, "SST", weightfile, (Int(1),), regrid_space)[1]
    SST_init = swap_space!(SST_init, axes(mask)) .* (abs.(mask .- 1)) .+ FT(273.15)

    ocean_params = OceanSlabParameters(FT(20), FT(1500.0), FT(800.0), FT(280.0), FT(1e-3), FT(1e-5), FT(0.06))
    ocean_sim = (; integrator = (; u = (; T_sfc = SST_init), p = (; params = ocean_params)))
    mode_specifics = (; name = mode_name)

    # sea ice
    # (Currently, we only support a slab ice model with fixed area and depth.)
    weightfile, datafile_cgll, regrid_space =
        ncreader_rll_to_cgll_from_space(sic_data, "SEAICE", boundary_space, outfile = "sic_cgll.nc")
    SIC = ncreader_cgll_sparse_to_field(datafile_cgll, "SEAICE", weightfile, (Int(1),), regrid_space)[1]
    SIC = swap_space!(SIC, axes(mask)) .* (abs.(mask .- 1))
    ice_mask = get_ice_mask.(SIC .- FT(25), FT) # here 25% and lower is considered ice free

    ice_sim = ice_init(FT; tspan = tspan, dt = Δt_cpl, space = boundary_space, saveat = saveat, ice_mask = ice_mask)

else
    # ocean
    ocean_sim = ocean_init(FT; tspan = tspan, dt = Δt_cpl, space = boundary_space, saveat = saveat, mask = mask)

    # sea ice
    ice_sim = (;
        integrator = (;
            u = (; T_sfc = ClimaCore.Fields.zeros(boundary_space)),
            p = (; params = ocean_sim.params, Ya = (; ice_mask = ClimaCore.Fields.zeros(boundary_space))),
        )
    )
    mode_specifics = (; name = mode_name)
end

# init coupler
coupler_field_names = (:T_S, :z0m_S, :z0b_S, :ρ_sfc, :q_sfc, :albedo, :F_A, :F_E, :F_R, :P_liq, :P_snow)
coupler_fields =
    NamedTuple{coupler_field_names}(ntuple(i -> ClimaCore.Fields.zeros(boundary_space), length(coupler_field_names)))
model_sims = (atm = atmos_sim, ice = ice_sim, lnd = land_sim, ocn = ocean_sim)
dates = (; date = [date], date0 = [date0], date1 = [date1])

coupler_sim = CouplerSimulation(
    FT(Δt_cpl),
    integrator.t,
    dates,
    boundary_space,
    FT,
    mask,
    coupler_fields,
    model_sims,
    mode_specifics,
    parsed_args,
);

include("./push_pull.jl")

# init conservation info collector
atmos_pull!(coupler_sim)
atmos_push!(coupler_sim)
land_pull!(coupler_sim)
reinit!(atmos_sim.integrator)
reinit!(land_sim.integrator)

mode_name == "amip" ? (ice_pull!(coupler_sim), reinit!(ice_sim.integrator)) :
(ocean_pull!(coupler_sim), reinit!(ocean_sim.integrator))

if !is_distributed && energy_check && mode_name == "aquaplanet"
    conservation_check = OnlineConservationCheck([], [], [], [], [], [])
    check_conservation(conservation_check, coupler_sim)
end

# coupling loop
@show "Starting coupling loop"

# step in time
walltime = @elapsed for t in ((tspan[1] + Δt_cpl):Δt_cpl:tspan[end])
    cs = coupler_sim

    date = current_date(t)

    @calendar_callback :(@show(date), date1 += Dates.Month(1)) date date1

    ## Atmos
    atmos_pull!(cs)
    step!(atmos_sim.integrator, t - atmos_sim.integrator.t, true) # NOTE: instead of Δt_cpl, to avoid accumulating roundoff error
    atmos_push!(cs)

    ## Slab land
    land_pull!(cs)
    step!(land_sim.integrator, t - land_sim.integrator.t, true)

    ## Slab ocean
    if cs.mode.name == "aquaplanet"
        ocean_pull!(cs)
        step!(ocean_sim.integrator, t - ocean_sim.integrator.t, true)
    end

    ## Slab ice
    if cs.mode.name == "amip"
        ice_pull!(cs)
        step!(ice_sim.integrator, t - ice_sim.integrator.t, true)
    end

    ## Compute energy
    if !is_distributed && energy_check && cs.mode.name == "aquaplanet"
        check_conservation(conservation_check, cs, atmos_sim, land_sim, ocean_sim, ice_sim)
    end

end

@show walltime

@show "Postprocessing"

if energy_check && coupler_sim.mode.name == "aquaplanet"
    plot_global_energy(
        conservation_check,
        coupler_sim,
        joinpath(coupler_output_dir, "total_energy_bucket.png"),
        joinpath(coupler_output_dir, "total_energy_log_bucket.png"),
    )
end

@show coupler_sim.fields.P_liq
# # animations
if (land_sim_name == "bucket") && parsed_args["anim"]
    #make it so this works with slab land?
    include("coupler_utils/viz_explorer.jl")
    plot_anim(
        coupler_sim.model_sims.atm,
        coupler_sim.model_sims.lnd,
        coupler_sim.model_sims.ocn,
        coupler_sim.model_sims.ice,
        coupler_sim.mask,
        coupler_sim.mode.name,
        coupler_sim.fields.T_S,
    )
end

# Cleanup temporary files
# TODO: Where should this live?
rm(REGRID_DIR; recursive = true, force = true)
