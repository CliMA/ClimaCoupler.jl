redirect_stderr(IOContext(stderr, :stacktrace_types_limited => Ref(false)))

import Dates
import DelimitedFiles
import YAML

import ClimaAtmos as CA
import ClimaComms
import ClimaCore as CC
import ClimaCalibrate

import ClimaCoupler
import ClimaCoupler:
    ConservationChecker, Checkpointer, FieldExchanger, FluxCalculator, Interfacer, TimeManager, Utilities
import ClimaCoupler.Interfacer:
    AbstractSlabplanetSimulationMode,
    AMIPMode,
    SlabplanetAquaMode,
    SlabplanetEisenmanMode,
    SlabplanetMode,
    SlabplanetTerraMode

import ClimaUtilities.SpaceVaryingInputs: SpaceVaryingInput
import ClimaUtilities.TimeVaryingInputs: TimeVaryingInput, evaluate!
import ClimaUtilities.Utils: period_to_seconds_float
import ClimaUtilities.ClimaArtifacts: @clima_artifact
import Interpolations # triggers InterpolationsExt in ClimaUtilities
# Random is used by RRMTGP for some cloud properties
import Random

# TODO: Move to ClimaUtilities once we move the Schedules to ClimaUtilities
import ClimaDiagnostics.Schedules: EveryCalendarDtSchedule

pkg_dir = pkgdir(ClimaCoupler)

## helpers for component models
include(joinpath(pkg_dir, "experiments/ClimaEarth/components/atmosphere/climaatmos.jl"))
include(joinpath(pkg_dir, "experiments/ClimaEarth/components/land/climaland_bucket.jl"))
include(joinpath(pkg_dir, "experiments/ClimaEarth/components/ocean/slab_ocean.jl"))
include(joinpath(pkg_dir, "experiments/ClimaEarth/components/ocean/prescr_seaice.jl"))
include(joinpath(pkg_dir, "experiments/ClimaEarth/components/ocean/eisenman_seaice.jl"))
include(joinpath(pkg_dir, "experiments/ClimaEarth/cli_options.jl"))
include(joinpath(pkg_dir, "experiments/ClimaEarth/user_io/arg_parsing.jl"))

function ClimaCalibrate.forward_model(iter, member)
    default_coupler_dict = parse_commandline(argparse_settings())
    config_file = joinpath(pkg_dir, "config/nightly_configs/amip_coarse_random.yml")
    config_dict = merge(default_coupler_dict, YAML.load_file(config_file))

    # Select the correct timestep for each component model based on which are available
    parse_component_dts!(config_dict)
    # Add extra diagnostics if specified
    add_extra_diagnostics!(config_dict)

    (;
        job_id,
        sim_mode,
        random_seed,
        FT,
        comms_ctx,
        t_end,
        t_start,
        date0,
        date,
        Δt_cpl,
        component_dt_dict,
        saveat,
        checkpoint_dt,
        restart_dir,
        restart_t,
        use_coupler_diagnostics,
        use_land_diagnostics,
        calendar_dt,
        evolving_ocean,
        mono_surface,
        turb_flux_partition,
        land_domain_type,
        land_albedo_type,
        land_initial_condition,
        land_temperature_anomaly,
        energy_check,
        conservation_softfail,
        output_dir_root,
        plot_diagnostics,
    ) = get_coupler_args(config_dict)
        config_dict["calibration_toml"] = ClimaCalibrate.parameter_path(output_dir_root, iter, member)
        random_seed = 1234

    @show output_dir_root
    member_output_dir = ClimaCalibrate.path_to_ensemble_member(output_dir_root, iter, member)
    dir_paths = Utilities.setup_output_dirs(output_dir = member_output_dir, comms_ctx = comms_ctx)
    @info "Coupler output directory $(dir_paths.output)"
    @info "Coupler artifacts directory $(dir_paths.artifacts)"
    @info "Coupler checkpoint directory $(dir_paths.checkpoints)"

    atmos_output_dir = joinpath(dir_paths.output, "clima_atmos")
    isdir(atmos_output_dir) || mkpath(atmos_output_dir)
    land_output_dir = joinpath(dir_paths.output, "clima_land")
    isdir(land_output_dir) || mkpath(land_output_dir)


    ## get component model dictionaries (if applicable)
    ## Note this step must come after parsing the coupler config dictionary, since
    ##  some parameters are passed from the coupler config to the component model configs
    job_id = "AMIP_calibration"
    atmos_config_dict = get_atmos_config_dict(config_dict, job_id, atmos_output_dir)
    (; dt_rad, output_default_diagnostics) = get_atmos_args(atmos_config_dict)

    ## set unique random seed if desired, otherwise use default
    Random.seed!(random_seed)
    @info "Random seed set to $(random_seed)"

    tspan = (t_start, t_end)

    #=
    ## Data File Paths
    =#
    sst_data, sic_data = try
        joinpath(@clima_artifact("historical_sst_sic", comms_ctx), "MODEL.SST.HAD187001-198110.OI198111-202206.nc"),
        joinpath(@clima_artifact("historical_sst_sic", comms_ctx), "MODEL.ICE.HAD187001-198110.OI198111-202206.nc")
    catch error
        @warn "Using lowres sst sic. If you want the higher resolution version, you have to obtain it from ClimaArtifacts"
        joinpath(
            @clima_artifact("historical_sst_sic_lowres", comms_ctx),
            "MODEL.SST.HAD187001-198110.OI198111-202206_lowres.nc",
        ),
        joinpath(
            @clima_artifact("historical_sst_sic_lowres", comms_ctx),
            "MODEL.ICE.HAD187001-198110.OI198111-202206_lowres.nc",
        )
    end
    co2_data = joinpath(@clima_artifact("co2_dataset", comms_ctx), "co2_mm_mlo.txt")
    land_mask_data = joinpath(@clima_artifact("landsea_mask_60arcseconds", comms_ctx), "landsea_mask.nc")

    Utilities.show_memory_usage()

    ## init atmos model component
    @info "TOML" atmos_config_dict["toml"]
    atmos_sim = atmos_init(CA.AtmosConfig(atmos_config_dict));
    # Get surface elevation from `atmos` coordinate field
    surface_elevation = CC.Fields.level(CC.Fields.coordinate_field(atmos_sim.integrator.u.f).z, CC.Utilities.half)
    Utilities.show_memory_usage()

    thermo_params = get_thermo_params(atmos_sim) # TODO: this should be shared by all models #342

    boundary_space = CC.Spaces.horizontal_space(atmos_sim.domain.face_space)

# Preprocess the file to be 1s and 0s before remapping into onto the grid
land_area_fraction = SpaceVaryingInput(land_mask_data, "landsea", boundary_space)
if !mono_surface
    land_area_fraction = Utilities.binary_mask.(land_area_fraction)
end
Utilities.show_memory_usage()


@info(sim_mode)
if sim_mode <: AMIPMode
    @info("AMIP boundary conditions - do not expect energy conservation")

    ## land model
    land_sim = bucket_init(
        FT,
        tspan,
        land_domain_type,
        land_albedo_type,
        land_initial_condition,
        land_temperature_anomaly,
        land_output_dir;
        dt = component_dt_dict["dt_land"],
        space = boundary_space,
        saveat = saveat,
        area_fraction = land_area_fraction,
        date_ref = date0,
        t_start = t_start,
        energy_check = energy_check,
        surface_elevation,
        use_land_diagnostics,
    )

    ## ocean stub
    SST_timevaryinginput = TimeVaryingInput(
        sst_data,
        "SST",
        boundary_space,
        reference_date = date0,
        file_reader_kwargs = (; preprocess_func = (data) -> data + FT(273.15),), ## convert to Kelvin
    )

    SST_init = zeros(boundary_space)
    evaluate!(SST_init, SST_timevaryinginput, t_start)

    ocean_sim = Interfacer.SurfaceStub((;
        T_sfc = SST_init,
        ρ_sfc = zeros(boundary_space),
        # ocean roughness follows GFDL model
        # (https://github.com/NOAA-GFDL/ice_param/blob/main/ocean_rough.F90#L47)
        z0m = FT(5.8e-5),
        z0b = FT(5.8e-5),
        beta = FT(1),
        α_direct = ones(boundary_space) .* FT(0.06),
        α_diffuse = ones(boundary_space) .* FT(0.06),
        area_fraction = (FT(1) .- land_area_fraction),
        phase = TD.Liquid(),
        thermo_params = thermo_params,
    ))

    ## sea ice model
    SIC_timevaryinginput = TimeVaryingInput(
        sic_data,
        "SEAICE",
        boundary_space,
        reference_date = date0,
        file_reader_kwargs = (; preprocess_func = (data) -> data / 100,), ## convert to fraction
    )

    SIC_init = zeros(boundary_space)
    evaluate!(SIC_init, SIC_timevaryinginput, t_start)

    ice_fraction = get_ice_fraction.(SIC_init, mono_surface)
    ice_sim = ice_init(
        FT;
        tspan = tspan,
        dt = component_dt_dict["dt_seaice"],
        space = boundary_space,
        saveat = saveat,
        area_fraction = ice_fraction,
        thermo_params = thermo_params,
    )

    ## CO2 concentration from temporally varying file
    CO2_text = DelimitedFiles.readdlm(co2_data, Float64; comments = true)
    # The text file only has month and year, so we set the day to 15th of the month
    years = CO2_text[:, 1]
    months = CO2_text[:, 2]
    CO2_dates = Dates.DateTime.(years, months) + Dates.Day(14)
    CO2_times = period_to_seconds_float.(CO2_dates .- date0)
    # convert from ppm to fraction, data is in fourth column of the text file
    CO2_vals = CO2_text[:, 4] .* 10^(-6)
    CO2_timevaryinginput = TimeVaryingInput(CO2_times, CO2_vals;)

    CO2_init = zeros(boundary_space)
    evaluate!(CO2_init, CO2_timevaryinginput, t_start)
    CO2_field = Interfacer.update_field!(atmos_sim, Val(:co2), CO2_init)

    mode_specifics = (;
        type = sim_mode,
        SST_timevaryinginput = SST_timevaryinginput,
        SIC_timevaryinginput = SIC_timevaryinginput,
        CO2_timevaryinginput = CO2_timevaryinginput,
    )
    Utilities.show_memory_usage()
end

## coupler exchange fields
coupler_field_names = (
    :T_S,
    :z0m_S,
    :z0b_S,
    :ρ_sfc,
    :q_sfc,
    :surface_direct_albedo,
    :surface_diffuse_albedo,
    :beta,
    :F_turb_energy,
    :F_turb_moisture,
    :F_turb_ρτxz,
    :F_turb_ρτyz,
    :F_radiative,
    :P_liq,
    :P_snow,
    :radiative_energy_flux_toa,
    :P_net,
    :temp1,
    :temp2,
)
coupler_fields = NamedTuple{coupler_field_names}(ntuple(i -> zeros(boundary_space), length(coupler_field_names)))
Utilities.show_memory_usage()

## model simulations
model_sims = (atmos_sim = atmos_sim, ice_sim = ice_sim, land_sim = land_sim, ocean_sim = ocean_sim);

## dates
dates = (; date = [date], date0 = [date0])

conservation_checks = nothing
schedule_checkpoint = EveryCalendarDtSchedule(TimeManager.time_to_period(checkpoint_dt); start_date = date0)
checkpoint_cb = TimeManager.TimeManager.Callback(schedule_checkpoint, Checkpointer.checkpoint_sims)

schedule_albedo = EveryCalendarDtSchedule(TimeManager.time_to_period(dt_rad); start_date = date0)

albedo_cb = TimeManager.Callback(schedule_albedo, FluxCalculator.water_albedo_from_atmosphere!)

callbacks = (; checkpoint = checkpoint_cb, water_albedo = albedo_cb)

turbulent_fluxes = nothing
if turb_flux_partition == "PartitionedStateFluxes"
    turbulent_fluxes = FluxCalculator.PartitionedStateFluxes()
elseif turb_flux_partition == "CombinedStateFluxesMOST"
    turbulent_fluxes = FluxCalculator.CombinedStateFluxesMOST()
else
    error("turb_flux_partition must be either PartitionedStateFluxes or CombinedStateFluxesMOST")
end

#= Set up default AMIP diagnostics
Use ClimaDiagnostics for default AMIP diagnostics, which currently include turbulent energy fluxes.
=#
# if sim_mode <: AMIPMode && use_coupler_diagnostics
#     include("user_io/amip_diagnostics.jl")
#     coupler_diags_path = joinpath(dir_paths.output, "coupler")
#     isdir(coupler_diags_path) || mkpath(coupler_diags_path)
#     amip_diags_handler =
#         amip_diagnostics_setup(coupler_fields, coupler_diags_path, dates.date0[1], tspan[1], calendar_dt)
# else
#     amip_diags_handler = nothing
# end
amip_diags_handler = nothing

cs = Interfacer.CoupledSimulation{FT}(
    comms_ctx,
    dates,
    boundary_space,
    coupler_fields,
    conservation_checks,
    [tspan[1], tspan[2]],
    Δt_cpl,
    model_sims,
    mode_specifics,
    callbacks,
    dir_paths,
    turbulent_fluxes,
    thermo_params,
    amip_diags_handler,
);
Utilities.show_memory_usage()

if !isnothing(restart_dir)
    for sim in cs.model_sims
        if Checkpointer.get_model_prog_state(sim) !== nothing
            Checkpointer.restart_model_state!(sim, comms_ctx, restart_t; input_dir = restart_dir)
        end
    end
end

FieldExchanger.update_surface_fractions!(cs)

# 2.surface density (`ρ_sfc`): calculated by the coupler by adiabatically extrapolating atmospheric thermal state to the surface.
# For this, we need to import surface and atmospheric fields. The model sims are then updated with the new surface density.
FieldExchanger.import_combined_surface_fields!(cs.fields, cs.model_sims, cs.turbulent_fluxes)
FieldExchanger.import_atmos_fields!(cs.fields, cs.model_sims, cs.boundary_space, cs.turbulent_fluxes)
FieldExchanger.update_model_sims!(cs.model_sims, cs.fields, cs.turbulent_fluxes)

# 3.surface vapor specific humidity (`q_sfc`): step surface models with the new surface density to calculate their respective `q_sfc` internally
## TODO: the q_sfc calculation follows the design of the bucket q_sfc, but it would be neater to abstract this from step! (#331)
Interfacer.step!(land_sim, Δt_cpl)
Interfacer.step!(ocean_sim, Δt_cpl)
Interfacer.step!(ice_sim, Δt_cpl)

# 4.turbulent fluxes: now we have all information needed for calculating the initial turbulent
# surface fluxes using either the combined state or the partitioned state method
if cs.turbulent_fluxes isa FluxCalculator.CombinedStateFluxesMOST
    ## import the new surface properties into the coupler (note the atmos state was also imported in step 3.)
    FieldExchanger.import_combined_surface_fields!(cs.fields, cs.model_sims, cs.turbulent_fluxes) # i.e. T_sfc, albedo, z0, beta, q_sfc
    ## calculate turbulent fluxes inside the atmos cache based on the combined surface state in each grid box
    FluxCalculator.combined_turbulent_fluxes!(cs.model_sims, cs.fields, cs.turbulent_fluxes) # this updates the atmos thermo state, sfc_ts
elseif cs.turbulent_fluxes isa FluxCalculator.PartitionedStateFluxes
    ## calculate turbulent fluxes in surface models and save the weighted average in coupler fields
    FluxCalculator.partitioned_turbulent_fluxes!(
        cs.model_sims,
        cs.fields,
        cs.boundary_space,
        FluxCalculator.MoninObukhovScheme(),
        cs.thermo_params,
    )

    ## update atmos sfc_conditions for surface temperature
    ## TODO: this is hard coded and needs to be simplified (req. CA modification) (#479)
    new_p = get_new_cache(atmos_sim, cs.fields)
    CA.SurfaceConditions.update_surface_conditions!(atmos_sim.integrator.u, new_p, atmos_sim.integrator.t) ## sets T_sfc (but SF calculation not necessary - requires split functionality in CA)
    atmos_sim.integrator.p.precomputed.sfc_conditions .= new_p.precomputed.sfc_conditions
end

# 5.reinitialize models + radiative flux: prognostic states and time are set to their initial conditions. For atmos, this also triggers the callbacks and sets a nonzero radiation flux (given the new sfc_conditions)
FieldExchanger.reinit_model_sims!(cs.model_sims)

# 6.update all fluxes: coupler re-imports updated atmos fluxes (radiative fluxes for both `turbulent_fluxes` types
# and also turbulent fluxes if `turbulent_fluxes isa CombinedStateFluxesMOST`,
# and sends them to the surface component models. If `turbulent_fluxes isa PartitionedStateFluxes`
# atmos receives the turbulent fluxes from the coupler.
FieldExchanger.import_atmos_fields!(cs.fields, cs.model_sims, cs.boundary_space, cs.turbulent_fluxes)
FieldExchanger.update_model_sims!(cs.model_sims, cs.fields, cs.turbulent_fluxes)

function solve_coupler!(cs)
    (; model_sims, Δt_cpl, tspan, comms_ctx) = cs
    (; atmos_sim, land_sim, ocean_sim, ice_sim) = model_sims

    @info("Starting coupling loop")
    ## step in time
    for t in ((tspan[begin] + Δt_cpl):Δt_cpl:tspan[end])
        # Update date
        cs.dates.date[] = TimeManager.current_date(cs, t)

        if cs.mode.type <: AMIPMode
            evaluate!(Interfacer.get_field(ocean_sim, Val(:surface_temperature)), cs.mode.SST_timevaryinginput, t)
            evaluate!(Interfacer.get_field(ice_sim, Val(:area_fraction)), cs.mode.SIC_timevaryinginput, t)

            # TODO: get_field with :co2 is not implemented, so this is a little awkward
            current_CO2 = zeros(boundary_space)
            evaluate!(current_CO2, cs.mode.CO2_timevaryinginput, t)
            Interfacer.update_field!(atmos_sim, Val(:co2), current_CO2)
        end

        ## compute global energy and water conservation checks
        ## (only for slabplanet if tracking conservation is enabled)
        !isnothing(cs.conservation_checks) && ConservationChecker.check_conservation!(cs)
        ClimaComms.barrier(comms_ctx)

        ## update water albedo from wind at dt_water_albedo
        ## (this will be extended to a radiation callback from the coupler)
        TimeManager.maybe_trigger_callback(cs.callbacks.water_albedo, cs, t)

        ## update the surface fractions for surface models,
        ## and update all component model simulations with the current fluxes stored in the coupler
        FieldExchanger.update_surface_fractions!(cs)
        FieldExchanger.update_model_sims!(cs.model_sims, cs.fields, cs.turbulent_fluxes)

        ## step component model simulations sequentially for one coupling timestep (Δt_cpl)
        FieldExchanger.step_model_sims!(cs.model_sims, t)

        ## update the coupler with the new surface properties and calculate the turbulent fluxes
        FieldExchanger.import_combined_surface_fields!(cs.fields, cs.model_sims, cs.turbulent_fluxes) # i.e. T_sfc, surface_albedo, z0, beta
        if cs.turbulent_fluxes isa FluxCalculator.CombinedStateFluxesMOST
            FluxCalculator.combined_turbulent_fluxes!(cs.model_sims, cs.fields, cs.turbulent_fluxes) # this updates the surface thermo state, sfc_ts, in ClimaAtmos (but also unnecessarily calculates fluxes)
        elseif cs.turbulent_fluxes isa FluxCalculator.PartitionedStateFluxes
            ## calculate turbulent fluxes in surfaces and save the weighted average in coupler fields
            FluxCalculator.partitioned_turbulent_fluxes!(
                cs.model_sims,
                cs.fields,
                cs.boundary_space,
                FluxCalculator.MoninObukhovScheme(),
                cs.thermo_params,
            )

            ## update atmos sfc_conditions for surface temperature - TODO: this needs to be simplified (need CA modification)
            new_p = get_new_cache(atmos_sim, cs.fields)
            CA.SurfaceConditions.update_surface_conditions!(atmos_sim.integrator.u, new_p, atmos_sim.integrator.t) # to set T_sfc (but SF calculation not necessary - CA modification)
            atmos_sim.integrator.p.precomputed.sfc_conditions .= new_p.precomputed.sfc_conditions
        end

        ## update the coupler with the new atmospheric properties
        FieldExchanger.import_atmos_fields!(cs.fields, cs.model_sims, cs.boundary_space, cs.turbulent_fluxes) # radiative and/or turbulent

        ## callback to checkpoint model state
        TimeManager.maybe_trigger_callback(cs.callbacks.checkpoint, cs, t)

        ## compute/output AMIP diagnostics if scheduled for this timestep
        ## wrap the current CoupledSimulation fields and time in a NamedTuple to match the ClimaDiagnostics interface
        cs_nt = (; u = cs.fields, p = nothing, t = t, step = round(t / Δt_cpl))
        (cs.mode.type <: AMIPMode && !isnothing(cs.amip_diags_handler)) &&
            CD.orchestrate_diagnostics(cs_nt, cs.amip_diags_handler)
    end
    return nothing
end

walltime = ClimaComms.@elapsed comms_ctx.device begin
    s = CA.@timed_str begin
        solve_coupler!(cs)
    end
end
@info(walltime)

end
