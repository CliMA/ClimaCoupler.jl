# coupler_driver
# don't forget to run with threading: julia --project --threads 8 (MPI not that useful for debugging coarse runs)

import SciMLBase: step!
using OrdinaryDiffEq
using OrdinaryDiffEq: ODEProblem, solve, SSPRK33, savevalues!, Euler
using LinearAlgebra
import Test: @test
using ClimaCore.Utilities: half, PlusHalf
using ClimaCore: InputOutput
using Dates
using UnPack

using Pkg

include("cli_options.jl")
(s, parsed_args) = parse_commandline()

# read in some parsed args
mode_name = parsed_args["mode_name"]
energy_check = parsed_args["energy_check"]
const FT = parsed_args["FLOAT_TYPE"] == "Float64" ? Float64 : Float32
land_sim_name = "bucket"
t_end = FT(time_to_seconds(parsed_args["t_end"]))
tspan = (0, t_end)
Δt_cpl = FT(parsed_args["dt_cpl"])
saveat = time_to_seconds(parsed_args["dt_save_to_sol"])
date0 = date = DateTime(1979, 01, 01)
mono_surface = parsed_args["mono_surface"]

# overwrite some parsed args
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
COUPLER_OUTPUT_DIR = joinpath(pkg_dir, "experiments/AMIP/moist_mpi_earth/output", mode_name)
isdir(COUPLER_OUTPUT_DIR) ? nothing : mkpath(COUPLER_OUTPUT_DIR)
REGRID_DIR = joinpath(COUPLER_OUTPUT_DIR, "regrid_tmp/")

# get the paths to the necessary data files - land sea mask, sst map, sea ice concentration
include("artifacts.jl")

# import coupler utils
include("coupler_utils/flux_calculator.jl")
include("coupler_utils/conservation_checker.jl")
include("coupler_utils/regridder.jl")
include("coupler_utils/masker.jl")
include("coupler_utils/calendar_timer.jl")
include("coupler_utils/general_helper.jl")
include("coupler_utils/bcfile_reader.jl")
include("coupler_utils/output_dumper.jl")

# init MPI
# include("mpi/mpi_init.jl") # rTODO: requires disabling MPI spec in ClimaAtmos 

# init atmos model component
include("atmos/atmos_init.jl")
atmos_sim = atmos_init(FT, Y, integrator, params = params);

# init a 2D bounary space at the surface, assuming the same instance (and MPI distribution if applicable) as the atmos domain above
boundary_space = atmos_sim.domain.face_space.horizontal_space # global surface grid

# init land-sea mask

land_mask = LandSeaMask(FT, mask_data, "LSMASK", boundary_space, mono = mono_surface)

# init surface (slab) model components
include("slab/slab_utils.jl")
include("bucket/bucket_init.jl")
include("slab/slab_init.jl")
include("slab_ocean/slab_init.jl")
include("slab_ice/slab_init.jl")

land_sim =
    bucket_init(FT, FT.(tspan), parsed_args["config"]; dt = FT(Δt_cpl), space = boundary_space, saveat = FT(saveat))

@info mode_name
if mode_name == "amip"
    println("No ocean sim - do not expect energy conservation")

    # ocean
    SST_info = bcfile_info_init(
        FT,
        sst_data,
        "SST",
        boundary_space,
        interpolate_daily = true,
        scaling_function = clean_sst, # convert to K
        land_mask = land_mask,
        date0 = date0,
        mono = mono_surface,
    )
    update_midmonth_data!(date0, SST_info)
    SST_init = interpolate_midmonth_to_daily(date0, SST_info)
    ocean_params = OceanSlabParameters(FT(20), FT(1500.0), FT(800.0), FT(280.0), FT(1e-3), FT(1e-5), FT(0.06))
    ocean_sim = (;
        integrator = (;
            u = (; T_sfc = SST_init),
            p = (; params = ocean_params, ocean_mask = (FT(1) .- land_mask)),
            SST_info = SST_info,
        )
    )

    # sea ice
    SIC_info = bcfile_info_init(
        FT,
        sic_data,
        "SEAICE",
        boundary_space,
        interpolate_daily = false,
        scaling_function = clean_sic, # convert to fractions
        land_mask = land_mask,
        date0 = date0,
        mono = mono_surface,
    )
    update_midmonth_data!(date0, SIC_info)
    SIC_init = interpolate_midmonth_to_daily(date0, SIC_info)
    ice_mask = get_ice_mask.(SIC_init, mono_surface)
    ice_sim = ice_init(FT; tspan = tspan, dt = Δt_cpl, space = boundary_space, saveat = saveat, ice_mask = ice_mask)
    mode_specifics = (; name = mode_name, SST_info = SST_info, SIC_info = SIC_info)

elseif mode_name == "slabplanet"
    # slab ocean for slabplanet runs
    # NB: ocean mask includes areas covered by sea ice
    ocean_sim = ocean_init(
        FT;
        tspan = tspan,
        dt = Δt_cpl,
        space = boundary_space,
        saveat = saveat,
        ocean_mask = (FT(1) .- land_mask),
    )

    ice_sim = (;
        integrator = (;
            u = (; T_sfc = ClimaCore.Fields.ones(boundary_space)),
            p = (; params = ocean_sim.params, ice_mask = ClimaCore.Fields.zeros(boundary_space)),
        )
    )
    mode_specifics = (; name = mode_name, SST_info = nothing, SIC_info = nothing)
end

# init coupler
coupler_field_names = (:T_S, :z0m_S, :z0b_S, :ρ_sfc, :q_sfc, :albedo, :F_A, :F_E, :F_R, :P_liq, :P_snow)
coupler_fields =
    NamedTuple{coupler_field_names}(ntuple(i -> ClimaCore.Fields.zeros(boundary_space), length(coupler_field_names)))
model_sims = (atmos_sim = atmos_sim, ice_sim = ice_sim, land_sim = land_sim, ocean_sim = ocean_sim);
dates = (; date = [date], date0 = [date0], date1 = [Dates.firstdayofmonth(date0)])

# online diagnostics
monthly_state_diags_names = (:ρ, :ρe_tot, :ρq_tot) # TODO: extend for :uₕ and cache variables
monthly_state_diags = (;
    fields = NamedTuple{monthly_state_diags_names}(
        ntuple(i -> ClimaCore.Fields.zeros(atmos_sim.domain.center_space), length(monthly_state_diags_names)),
    ),
    ct = [0],
)

cs = CouplerSimulation(
    FT(Δt_cpl),
    integrator.t,
    tspan,
    dates,
    boundary_space,
    FT,
    (; land = land_mask, ocean = zeros(boundary_space), ice = zeros(boundary_space)),
    coupler_fields,
    model_sims,
    mode_specifics,
    parsed_args,
    monthly_state_diags,
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

if !simulation.is_distributed && energy_check && mode_name == "slabplanet"
    conservation_check = OnlineConservationCheck([], [], [], [], [], [])
    check_conservation(conservation_check, cs)
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
        @calendar_callback :(@show(cs.dates.date[1])) cs.dates.date[1] cs.dates.date1[1]

        if cs.mode.name == "amip"

            # monthly read of boundary condition data for sea surface temperature (SST) and sea ice concentration (SIC)
            @calendar_callback :(update_midmonth_data!(cs.dates.date[1], cs.mode.SST_info)) cs.dates.date[1] next_date_in_file(
                cs.mode.SST_info,
            )
            SST = ocean_sim.integrator.u.T_sfc .= interpolate_midmonth_to_daily(cs.dates.date[1], cs.mode.SST_info)
            @calendar_callback :(update_midmonth_data!(cs.dates.date[1], cs.mode.SIC_info)) cs.dates.date[1] next_date_in_file(
                cs.mode.SIC_info,
            )
            SIC = interpolate_midmonth_to_daily(cs.dates.date[1], cs.mode.SIC_info)

            ice_mask = ice_sim.integrator.p.ice_mask .= get_ice_mask.(SIC_init, mono_surface)


            # accumulate data at each timestep
            accumulate_diags(atmos_sim.integrator.u.c, cs.monthly_state_diags)

            # save and reset monthly averages
            @calendar_callback :(
                map(x -> x ./= cs.monthly_state_diags.ct[1], cs.monthly_state_diags.fields),
                save_hdf5(cs.monthly_state_diags.fields, cs.dates.date[1], COUPLER_OUTPUT_DIR),
                map(x -> x .= FT(0), cs.monthly_state_diags.fields),
                cs.monthly_state_diags.ct .= FT(0),
            ) cs.dates.date[1] cs.dates.date1[1]

        end

        ## atmos
        atmos_pull!(cs)
        step!(atmos_sim.integrator, t - atmos_sim.integrator.t, true) # NOTE: instead of Δt_cpl, to avoid accumulating roundoff error

        atmos_push!(cs)

        ## land
        land_pull!(cs)
        step!(land_sim.integrator, t - land_sim.integrator.t, true)

        ## ocean
        if cs.mode.name == "slabplanet"
            ocean_pull!(cs)
            step!(ocean_sim.integrator, t - ocean_sim.integrator.t, true)
        end

        ## sea ice
        if cs.mode.name == "amip"
            ice_pull!(cs)
            step!(ice_sim.integrator, t - ice_sim.integrator.t, true)
        end

        ## compute energy
        if !simulation.is_distributed && energy_check && cs.mode.name == "slabplanet"
            check_conservation(conservation_check, cs)
        end

        ## step to next calendar month
        @calendar_callback :(cs.dates.date1[1] += Dates.Month(1)) cs.dates.date[1] cs.dates.date1[1]

    end
    @show walltime

    return cs
end

solve_coupler!(cs, energy_check);

@info "Postprocessing"

isdir(COUPLER_OUTPUT_DIR * "_artifacts") ? nothing : mkpath(COUPLER_OUTPUT_DIR * "_artifacts")

if energy_check && cs.mode.name == "slabplanet"
    @info "Energy Check"
    plot_global_energy(
        conservation_check,
        cs,
        joinpath(COUPLER_OUTPUT_DIR * "_artifacts", "total_energy_bucket.png"),
        joinpath(COUPLER_OUTPUT_DIR * "_artifacts", "total_energy_log_bucket.png"),
    )
end

# # animations
if parsed_args["anim"]
    @info "Animations"
    include("coupler_utils/viz_explorer.jl")
    plot_anim(cs, COUPLER_OUTPUT_DIR * "_artifacts")
end

# unit tests (to be moved to `test/`)
if cs.mode.name == "amip"
    @info "Unit Tests"
    include("coupler_utils/unit_tester.jl")
end

# Cleanup temporary files
# TODO: Where should this live?
rm(COUPLER_OUTPUT_DIR; recursive = true, force = true)
# - cs needs to be global for the monthly macro - explote other solutions
# - SST_init is modified with SST_info even with deepcopy...
# - replace if statements with dipatches, write better abstractions
# - add spin up for AMIP
# - do nothing for the parts of the domain that are masked out - aware of MPI and fracitonal surfaces  
