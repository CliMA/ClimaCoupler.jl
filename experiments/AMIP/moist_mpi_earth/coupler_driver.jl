# coupler_driver
# don't forget to run with threading: julia --project --threads 8 (MPI not that useful for debugging coarse runs)

# setup MPI context if distributed
include("mpi/mpi_init.jl")

import SciMLBase: step!
using OrdinaryDiffEq
using OrdinaryDiffEq: ODEProblem, solve, SSPRK33, savevalues!, Euler
using LinearAlgebra
import Test: @test
using ClimaCore.Utilities: half, PlusHalf
using ClimaCore: InputOutput
using Dates
using UnPack

include("cli_options.jl")
(s, parsed_args) = parse_commandline()

# read in some parsed args
mode_name = parsed_args["mode_name"]
run_name = parsed_args["run_name"]
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
COUPLER_OUTPUT_DIR = joinpath(pkg_dir, "experiments/AMIP/moist_mpi_earth/output", joinpath(mode_name, run_name))
!isdir(COUPLER_OUTPUT_DIR) && mkpath(COUPLER_OUTPUT_DIR)
REGRID_DIR = joinpath(COUPLER_OUTPUT_DIR, "regrid_tmp/")
mkpath(REGRID_DIR)

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
include("coupler_utils/variable_definer.jl")
include("coupler_utils/diagnostics_gatherer.jl")
include("coupler_utils/offline_postprocessor.jl")

# init atmos model component
include("atmos/atmos_init.jl")
atmos_sim = atmos_init(FT, Y, integrator, params = params);

# init a 2D bounary space at the surface, assuming the same instance (and MPI distribution if applicable) as the atmos domain above
boundary_space = atmos_sim.domain.face_space.horizontal_space # global surface grid

# init land-sea mask
land_mask = LandSeaMask(FT, comms_ctx, mask_data, "LSMASK", boundary_space, mono = mono_surface)

# init surface (slab) model components
include("slab/slab_utils.jl")
include("bucket/bucket_init.jl")
include("slab/slab_init.jl")
include("slab_ocean/slab_init.jl")
include("slab_ice/slab_init.jl")

# 1. land
land_sim =
    bucket_init(FT, FT.(tspan), parsed_args["config"]; dt = FT(Δt_cpl), space = boundary_space, saveat = FT(saveat))

@info mode_name
if mode_name == "amip"
    println("AMIP boundary conditions - do not expect energy conservation")

    # 2. ocean
    SST_info = bcfile_info_init(
        FT,
        comms_ctx,
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
    # 3. sea ice
    SIC_info = bcfile_info_init(
        FT,
        comms_ctx,
        sic_data,
        "SEAICE",
        boundary_space,
        interpolate_daily = true,
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
    # 2. ocean
    ocean_sim = ocean_init(
        FT;
        tspan = tspan,
        dt = Δt_cpl,
        space = boundary_space,
        saveat = saveat,
        ocean_mask = (FT(1) .- land_mask), # NB: ocean mask includes areas covered by sea ice
    )

    # 3. sea ice
    ice_sim = (;
        integrator = (;
            u = (; T_sfc = ClimaCore.Fields.ones(boundary_space)),
            p = (; params = ocean_sim.params, ice_mask = ClimaCore.Fields.zeros(boundary_space)),
        )
    )
    mode_specifics = (; name = mode_name, SST_info = nothing, SIC_info = nothing)
end

# init coupler 

# 1. coupler fields
coupler_field_names = (:T_S, :z0m_S, :z0b_S, :ρ_sfc, :q_sfc, :albedo, :F_A, :F_E, :F_R, :P_liq, :P_snow)
coupler_fields =
    NamedTuple{coupler_field_names}(ntuple(i -> ClimaCore.Fields.zeros(boundary_space), length(coupler_field_names)))

# 2. model simulations
model_sims = (atmos_sim = atmos_sim, ice_sim = ice_sim, land_sim = land_sim, ocean_sim = ocean_sim);

# 3. dates
dates = (; date = [date], date0 = [date0], date1 = [Dates.firstdayofmonth(date0)])

# 4. online diagnostics
monthly_3d_diags_names = (:T, :u, :q_tot) # need to be specified in coupler_utils/variable_definer.jl
monthly_3d_diags = (;
    fields = NamedTuple{monthly_3d_diags_names}(
        ntuple(i -> ClimaCore.Fields.zeros(atmos_sim.domain.center_space), length(monthly_3d_diags_names)),
    ),
    ct = [0],
)
monthly_2d_diags_names = (:precipitation, :toa, :T_sfc) # need to be specified in coupler_utils/variable_definer.jl
monthly_2d_diags = (;
    fields = NamedTuple{monthly_2d_diags_names}(
        ntuple(i -> ClimaCore.Fields.zeros(boundary_space), length(monthly_2d_diags_names)),
    ),
    ct = [0],
)

# 5. coupler simulation
cs = CouplerSimulation(
    comms_ctx,
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
    monthly_3d_diags,
    monthly_2d_diags,
);

# share states between models
include("./push_pull.jl")
atmos_pull!(cs)
if parsed_args["ode_algo"] == "ARS343"
    step!(atmos_sim.integrator, Δt_cpl, true) # this is necessary to set values to the unitialized cache. In `ODE.jl` this is done as part of `reinit!``
end
atmos_push!(cs)
land_pull!(cs)

# reinitialize (TODO: avoid with interfaces)
reinit!(atmos_sim.integrator)
reinit!(land_sim.integrator)
mode_name == "amip" ? (ice_pull!(cs), reinit!(ice_sim.integrator)) : nothing
mode_name == "slabplanet" ? (ocean_pull!(cs), reinit!(ocean_sim.integrator)) : nothing

# init conservation info collector
if !is_distributed && energy_check && mode_name == "slabplanet"
    conservation_check = OnlineConservationCheck([], [], [], [], [], [], [])
    check_conservation(conservation_check, cs)
end

# coupling loop
function solve_coupler!(cs, energy_check)
    @info "Starting coupling loop"

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

            # accumulate diagnostics at each timestep
            accumulate_diags(collect_diags(cs, propertynames(cs.monthly_3d_diags.fields)), cs.monthly_3d_diags)
            accumulate_diags(collect_diags(cs, propertynames(cs.monthly_2d_diags.fields)), cs.monthly_2d_diags)

            # save and reset monthly averages
            @calendar_callback :(
                map(x -> x ./= cs.monthly_3d_diags.ct[1], cs.monthly_3d_diags.fields),
                save_hdf5(
                    cs.comms_ctx,
                    cs.monthly_3d_diags.fields,
                    cs.dates.date[1],
                    COUPLER_OUTPUT_DIR,
                    name_tag = "3d_",
                ),
                map(x -> x .= FT(0), cs.monthly_3d_diags.fields),
                cs.monthly_3d_diags.ct .= FT(0),
            ) cs.dates.date[1] cs.dates.date1[1]
            @calendar_callback :(
                map(x -> x ./= cs.monthly_2d_diags.ct[1], cs.monthly_2d_diags.fields),
                save_hdf5(
                    cs.comms_ctx,
                    cs.monthly_2d_diags.fields,
                    cs.dates.date[1],
                    COUPLER_OUTPUT_DIR,
                    name_tag = "2d_",
                ),
                map(x -> x .= FT(0), cs.monthly_2d_diags.fields),
                cs.monthly_2d_diags.ct .= FT(0),
            ) cs.dates.date[1] cs.dates.date1[1]

        end

        # run component models sequentially for one coupling timestep (Δt_cpl)
        # 1. atmos
        ClimaComms.barrier(comms_ctx)

        atmos_pull!(cs)
        step!(atmos_sim.integrator, t - atmos_sim.integrator.t, true) # NOTE: instead of Δt_cpl, to avoid accumulating roundoff error
        atmos_push!(cs)

        # 2. land
        land_pull!(cs)
        step!(land_sim.integrator, t - land_sim.integrator.t, true)

        # 3. ocean
        if cs.mode.name == "slabplanet"
            ocean_pull!(cs)
            step!(ocean_sim.integrator, t - ocean_sim.integrator.t, true)
        end

        # 4. sea ice
        if cs.mode.name == "amip"
            ice_pull!(cs)
            step!(ice_sim.integrator, t - ice_sim.integrator.t, true)
        end

        # compute global energy
        if !simulation.is_distributed && energy_check && cs.mode.name == "slabplanet"
            check_conservation(conservation_check, cs)
        end

        # step to the next calendar month
        @calendar_callback :(cs.dates.date1[1] += Dates.Month(1)) cs.dates.date[1] cs.dates.date1[1]

    end
    @show walltime

    return cs
end

# run the coupled simulation
solve_coupler!(cs, energy_check);

# postprocessing
if ClimaComms.iamroot(comms_ctx)
    isdir(COUPLER_OUTPUT_DIR * "_artifacts") ? nothing : mkpath(COUPLER_OUTPUT_DIR * "_artifacts")

    # 1. energy check plots
    if !is_distributed && energy_check && cs.mode.name == "slabplanet"
        @info "Energy Check"
        plot_global_energy(
            conservation_check,
            cs,
            joinpath(COUPLER_OUTPUT_DIR * "_artifacts", "total_energy_bucket.png"),
            joinpath(COUPLER_OUTPUT_DIR * "_artifacts", "total_energy_log_bucket.png"),
        )
    end

    # 2. animations
    if !is_distributed && parsed_args["anim"]
        @info "Animations"
        include("coupler_utils/viz_explorer.jl")
        plot_anim(cs, COUPLER_OUTPUT_DIR * "_artifacts")
    end

    # 3. plotting AMIP results
    if cs.mode.name == "amip"
        @info "AMIP plots"

        include("coupler_utils/plotter.jl")

        # ClimaESM
        include("coupler_utils/amip_visualizer.jl")
        post_spec = (;
            T = (:regridded_3d, :zonal_mean),
            u = (:regridded_3d, :zonal_mean),
            q_tot = (:regridded_3d, :zonal_mean),
            toa = (:regridded_2d, :horizontal_2d),
            precipitation = (:regridded_2d, :horizontal_2d),
            T_sfc = (:regridded_2d, :horizontal_2d),
        )

        plot_spec = (; # optional plotting args
            T = (; clims = (190, 320), units = "K"),
            u = (; clims = (-50, 50), units = "m/s"),
            q_tot = (; clims = (0, 50), units = "g/kg"),
            toa = (; clims = (-250, 210), units = "W/m^2"),
            precipitation = (clims = (0, 1e-6), units = "kg/m^2/s"),
            T_sfc = (clims = (225, 310), units = "K"),
        )
        amip_paperplots(
            post_spec,
            plot_spec,
            COUPLER_OUTPUT_DIR,
            files_root = ".monthly",
            output_dir = COUPLER_OUTPUT_DIR * "_artifacts",
        )

        # NCEP reanalysis
        @info "NCEP plots"
        include("coupler_utils/ncep_visualizer.jl")
        ncep_post_spec = (;
            T = (:zonal_mean,),
            u = (:zonal_mean,),
            q_tot = (:zonal_mean,),
            toa = (:horizontal_2d,),
            precipitation = (:horizontal_2d,),
            T_sfc = (:horizontal_2d,),
        )
        ncep_plot_spec = plot_spec
        ncep_paperplots(
            ncep_post_spec,
            ncep_plot_spec,
            COUPLER_OUTPUT_DIR,
            output_dir = COUPLER_OUTPUT_DIR * "_artifacts",
        ) # plot data that correspond to the model's last save_hdf5 call (i.e., last month)
    end

    rm(COUPLER_OUTPUT_DIR; recursive = true, force = true) # TODO: Where should this live?
end

# unit tests (to be moved to `test/`)
if !is_distributed && cs.mode.name == "amip"
    @info "Unit Tests"
    include("coupler_utils/unit_tester.jl")
end

# Cleanup temporary files
# TODO:
# - cs needs to be global for the monthly macro - explore other solutions
# - SST_init is modified with SST_info even with deepcopy...
# - replace if statements with dipatches, write better abstractions
# - add spin up for AMIP
# - do nothing for the parts of the domain that are masked out - aware of MPI and fracitonal surfaces  
# - add an energy check for MPI runs
