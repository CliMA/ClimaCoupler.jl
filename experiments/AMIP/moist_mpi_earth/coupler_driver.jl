# coupler_driver
# don't forget to run with threading: julia --project --threads 8 (MPI not that useful for debugging coarse runs)

import SciMLBase: step!
using OrdinaryDiffEq
using OrdinaryDiffEq: ODEProblem, solve, SSPRK33, savevalues!, Euler
using LinearAlgebra
import Test: @test
using ClimaCore.Utilities: half, PlusHalf
using Dates

using Pkg
Pkg.add(PackageSpec(name = "ClimaAtmos", rev = "cli_options")) # remove when ClimaAtmos@0.4.0 is released

include("cli_options.jl")
(s, parsed_args) = parse_commandline()

# Read in some parsed args
mode_name = parsed_args["mode_name"]
energy_check = parsed_args["energy_check"]
const FT = parsed_args["FLOAT_TYPE"] == "Float64" ? Float64 : Float32
land_sim_name = "bucket"
t_end = FT(time_to_seconds(parsed_args["t_end"]))
tspan = (0, t_end)
Δt_cpl = FT(parsed_args["dt_cpl"])
saveat = time_to_seconds(parsed_args["dt_save_to_sol"])
date0 = date = DateTime(1979, 01, 01)

# overwrite some parsed args :P
parsed_args["coupled"] = true
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
include("coupler_utils/bcfile_reader.jl")

# init MPI
include("mpi/mpi_init.jl")

# init atmos model component
include("atmos/atmos_init.jl")
atmos_sim = atmos_init(FT, Y, integrator, params = params);

# init a 2D bounary space at the surface, assuming the same instance (and MPI distribution if applicable) as the atmos domain above
boundary_space = ClimaCore.Fields.level(atmos_sim.domain.face_space, half) # global surface grid

# init land-sea mask
land_mask = LandSeaMask(FT, mask_data, "LSMASK", boundary_space)

# init surface (slab) model components
# we need some types that are defined in these files
include("slab/slab_utils.jl")
include("bucket/bucket_init.jl")
include("slab/slab_init.jl")
include("slab_ocean/slab_init.jl")
include("slab_ice/slab_init.jl")

land_sim = bucket_init(FT, FT.(tspan); dt = FT(Δt_cpl), space = boundary_space, saveat = FT(saveat))

if mode_name == "amip"
    println("No ocean sim - do not expect energy conservation")

    # ocean
    SST_info = bcfile_info_init(
        FT,
        sst_data,
        "SST",
        boundary_space,
        interpolate_daily = false,
        scaling_function = clean_sst,
        land_mask = land_mask,
        date0 = date0,
    )
    update_midmonth_data!(date0, SST_info)
    SST_init = deepcopy(SST_info.monthly_fields[1])
    ocean_params = OceanSlabParameters(FT(20), FT(1500.0), FT(800.0), FT(280.0), FT(1e-3), FT(1e-5), FT(0.06))
    ocean_sim = (; integrator = (; u = (; T_sfc = SST_init), p = (; params = ocean_params), SST_info = SST_info))

    # sea ice
    SIC_info = bcfile_info_init(
        FT,
        sic_data,
        "SEAICE",
        boundary_space,
        interpolate_daily = false,
        scaling_function = clean_sic,
        land_mask = land_mask,
        date0 = date0,
    )
    update_midmonth_data!(date0, SIC_info)
    SIC_init = deepcopy(SIC_info.monthly_fields[1])
    ice_mask = get_ice_mask.(SIC_init, FT) # here 50% and lower is considered ice free
    ice_sim = ice_init(FT; tspan = tspan, dt = Δt_cpl, space = boundary_space, saveat = saveat, ice_mask = ice_mask)
    mode_specifics = (; name = mode_name, SST_info = SST_info, SIC_info = SIC_info)

elseif mode_name == "slabplanet"
    # slab ocean for slabplanet runs
    ocean_sim =
        ocean_init(FT; tspan = tspan, dt = Δt_cpl, space = boundary_space, saveat = saveat, land_mask = land_mask)

    ice_sim = (;
        integrator = (;
            u = (; T_sfc = ClimaCore.Fields.zeros(boundary_space)),
            p = (; params = ocean_sim.params, Ya = (; ice_mask = ClimaCore.Fields.zeros(boundary_space))),
        )
    )
    mode_specifics = (; name = mode_name, SST_info = nothing, SIC_info = nothing)
end

# init coupler
coupler_field_names = (:T_S, :z0m_S, :z0b_S, :ρ_sfc, :q_sfc, :albedo, :F_A, :F_E, :F_R, :P_liq)
coupler_fields =
    NamedTuple{coupler_field_names}(ntuple(i -> ClimaCore.Fields.zeros(boundary_space), length(coupler_field_names)))
model_sims = (atmos_sim = atmos_sim, ice_sim = ice_sim, land_sim = land_sim, ocean_sim = ocean_sim)
dates = (; date = [date], date0 = [date0], date1 = [Dates.firstdayofmonth(date0)])

cs = CouplerSimulation(
    FT(Δt_cpl),
    integrator.t,
    tspan,
    dates,
    boundary_space,
    FT,
    land_mask,
    coupler_fields,
    model_sims,
    mode_specifics,
    parsed_args,
);

include("./push_pull.jl")

# init conservation info collector
atmos_pull!(cs)
atmos_push!(cs)
land_pull!(cs)
reinit!(atmos_sim.integrator)
reinit!(land_sim.integrator)

mode_name == "amip" ? (ice_pull!(cs), reinit!(ice_sim.integrator)) : nothing
mode_name == "slabplanet" ? (ocean_pull!(cs), reinit!(ocean_sim.integrator)) : nothing

if !is_distributed && energy_check && mode_name == "slabplanet"
    conservation_check = OnlineConservationCheck([], [], [], [], [], [])
    check_conservation(conservation_check, cs, atmos_sim, land_sim, ocean_sim, ice_sim)
end

# coupling loop
function solve_coupler!(cs, energy_check)
    @show "Starting coupling loop"

    @unpack model_sims, Δt_cpl, tspan = cs
    @unpack atmos_sim, land_sim, ocean_sim, ice_sim = model_sims

    # step in time
    walltime = @elapsed for t in ((tspan[1] + Δt_cpl):Δt_cpl:tspan[end])

        cs.dates.date[1] = current_date(t) # if not global, `date` is not updated.

        # print date on the first of month 
        @calendar_callback :(@show(cs.dates.date[1]), cs.dates.date1[1] += Dates.Month(1)) cs.dates.date[1] cs.dates.date1[1]

        # monthly read of boundary condition data
        if cs.mode.name == "amip"
            @calendar_callback :(update_midmonth_data!(cs.dates.date[1], cs.mode.SST_info)) cs.dates.date[1] next_date_in_file(
                cs.mode.SST_info,
            )
            SST = ocean_sim.integrator.u.T_sfc .= cs.mode.SST_info.monthly_fields[1]
            @calendar_callback :(update_midmonth_data!(cs.dates.date[1], cs.mode.SIC_info)) cs.dates.date[1] next_date_in_file(
                cs.mode.SIC_info,
            )
            SIC = cs.mode.SIC_info.monthly_fields[1]
            ice_mask = ice_sim.integrator.p.Ya.ice_mask .= get_ice_mask.(SIC, FT)

        end

        ## Atmos
        atmos_pull!(cs)
        step!(atmos_sim.integrator, t - atmos_sim.integrator.t, true) # NOTE: instead of Δt_cpl, to avoid accumulating roundoff error

        atmos_push!(cs)

        ## Slab land
        land_pull!(cs)
        step!(land_sim.integrator, t - land_sim.integrator.t, true)

        ## Slab ocean
        if cs.mode.name == "slabplanet"
            ocean_pull!(cs)
            step!(ocean_sim.integrator, t - ocean_sim.integrator.t, true)
        end

        ## Slab ice
        if cs.mode.name == "amip"
            ice_pull!(cs)
            step!(ice_sim.integrator, t - ice_sim.integrator.t, true)
        end

        ## Compute energy
        if !is_distributed && energy_check && cs.mode.name == "slabplanet"
            check_conservation(conservation_check, cs, atmos_sim, land_sim, ocean_sim, ice_sim)
        end
    end
    @show walltime
    return cs
end

solve_coupler!(cs, energy_check);

@show "Postprocessing"
if energy_check && cs.mode.name == "slabplanet"
    plot_global_energy(
        conservation_check,
        cs,
        joinpath(coupler_output_dir, "total_energy_bucket.png"),
        joinpath(coupler_output_dir, "total_energy_log_bucket.png"),
    )
end

@show cs.fields.P_liq
# # animations
if parsed_args["anim"]
    #make it so this works with slab land?
    include("coupler_utils/viz_explorer.jl")
    plot_anim(
        cs.model_sims.atmos_sim,
        cs.model_sims.land_sim,
        cs.model_sims.ocean_sim,
        cs.model_sims.ice_sim,
        cs.land_mask,
        cs.mode.name,
        cs.fields.T_S,
    )
end

# unit tests (to be moved to `test/`)
include("coupler_utils/unit_tester.jl")

# Cleanup temporary files
# TODO: Where should this live?
rm(REGRID_DIR; recursive = true, force = true)
# - cs needs to be global for the monthly macro - explote other solutions
# - SST_init is modified with SST_info even with deepcopy...
# - replace if statements with dipatches, write better abstractions
# - unit test for monthly file update 
