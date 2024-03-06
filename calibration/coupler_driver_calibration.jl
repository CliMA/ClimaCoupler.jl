redirect_stderr(IOContext(stderr, :stacktrace_types_limited => Ref(false)))

using ClimaComms
comms_ctx = ClimaComms.context()
const pid, nprocs = ClimaComms.init(comms_ctx)


import SciMLBase: ODEProblem, solve, step!, init, reinit!
using LinearAlgebra
import Test: @test
using Dates
using Plots
using Statistics: mean
import ClimaAtmos as CA
import YAML

using ClimaCore.Utilities: half, PlusHalf
using ClimaCore: InputOutput, Fields
import ClimaCore.Spaces as Spaces

## coupler specific imports
import ClimaCoupler
import ClimaCoupler.Regridder
import ClimaCoupler.Regridder:
    update_surface_fractions!, combine_surfaces!, combine_surfaces_from_sol!, dummmy_remap!, binary_mask
import ClimaCoupler.ConservationChecker:
    EnergyConservationCheck, WaterConservationCheck, check_conservation!, plot_global_conservation
import ClimaCoupler.Utilities: swap_space!
import ClimaCoupler.BCReader:
    bcfile_info_init, float_type_bcf, update_midmonth_data!, next_date_in_file, interpolate_midmonth_to_daily
import ClimaCoupler.TimeManager:
    current_date,
    datetime_to_strdate,
    trigger_callback,
    Monthly,
    EveryTimestep,
    HourlyCallback,
    MonthlyCallback,
    update_firstdayofmonth!,
    trigger_callback!
import ClimaCoupler.Diagnostics: get_var, init_diagnostics, accumulate_diagnostics!, save_diagnostics, TimeMean
import ClimaCoupler.PostProcessor: postprocess

import ClimaCoupler.Interfacer:
    CoupledSimulation,
    float_type,
    AtmosModelSimulation,
    SurfaceModelSimulation,
    SurfaceStub,
    SeaIceModelSimulation,
    LandModelSimulation,
    OceanModelSimulation,
    get_field,
    update_field!
import ClimaCoupler.FluxCalculator:
    PartitionedStateFluxes,
    CombinedStateFluxes,
    combined_turbulent_fluxes!,
    MoninObukhovScheme,
    partitioned_turbulent_fluxes!
import ClimaCoupler.FieldExchanger:
    import_atmos_fields!,
    import_combined_surface_fields!,
    update_sim!,
    update_model_sims!,
    reinit_model_sims!,
    step_model_sims!
import ClimaCoupler.Checkpointer: checkpoint_model_state, get_model_state_vector, restart_model_state!

## helpers for component models
include("../experiments/AMIP/components/atmosphere/climaatmos_init.jl")
include("../experiments/AMIP/components/land/bucket_init.jl")
include("../experiments/AMIP/components/land/bucket_utils.jl")
include("../experiments/AMIP/components/ocean/slab_ocean_init.jl")
include("../experiments/AMIP/components/ocean/prescr_seaice_init.jl")
include("../experiments/AMIP/user_io/user_diagnostics.jl")
include("../experiments/AMIP/user_io/user_logging.jl")

## coupler defaults
# get component model dictionaries
include("../experiments/AMIP/cli_options.jl")
parsed_args = parse_commandline(argparse_settings())
config_dict = YAML.load_file("./experiments/amip_coupled/coupler_config.yml")
config_dict = YAML.load_file(joinpath(experiment_dir, "coupler_config.yml"));
config_dict["t_end"] = "150secs";
config_dict["output_dir"] = output_dir;
config_dict = merge(parsed_args, config_dict)
config_dict_atmos = get_atmos_config(config_dict)

# merge dictionaries of command line arguments, coupler dictionary and component model dictionaries
# (if there are common keys, the last dictorionary in the `merge` arguments takes precedence)
config_dict = merge(config_dict_atmos, config_dict)


## read in some parsed command line arguments
mode_name = config_dict["mode_name"]
run_name = config_dict["run_name"]
energy_check = config_dict["energy_check"]
FT = config_dict["FLOAT_TYPE"] == "Float64" ? Float64 : Float32
land_sim_name = "bucket"
t_end = Float64(time_to_seconds(config_dict["t_end"]))
t_start = 0.0
tspan = (t_start, t_end)
Δt_cpl = Float64(config_dict["dt_cpl"])
saveat = Float64(time_to_seconds(config_dict["dt_save_to_sol"]))
date0 = date = DateTime(config_dict["start_date"], dateformat"yyyymmdd")
mono_surface = config_dict["mono_surface"]
hourly_checkpoint = config_dict["hourly_checkpoint"]
restart_dir = config_dict["restart_dir"]
restart_t = Int(config_dict["restart_t"])
evolving_ocean = config_dict["evolving_ocean"]
config_dict["print_config_dict"] = false

## I/O directory setup
COUPLER_OUTPUT_DIR = "/Users/akshaysridhar/Research/Codes/ClimaCoupler.jl/calibration/output/amip/"
mkpath(COUPLER_OUTPUT_DIR)

REGRID_DIR = joinpath(COUPLER_OUTPUT_DIR, "regrid_tmp/")
mkpath(REGRID_DIR)

COUPLER_ARTIFACTS_DIR = COUPLER_OUTPUT_DIR * "_artifacts"
isdir(COUPLER_ARTIFACTS_DIR) ? nothing : mkpath(COUPLER_ARTIFACTS_DIR)

config_dict["print_config_dict"] ? @info(config_dict) : nothing

# get the paths to the necessary data files: land-sea mask, sst map, sea ice concentration
include(joinpath(pkgdir(ClimaCoupler), "artifacts", "artifact_funcs.jl"))
sst_data = joinpath(sst_dataset_path(), "sst.nc")
sic_data = joinpath(sic_dataset_path(), "sic.nc")
co2_data = joinpath(co2_dataset_path(), "mauna_loa_co2.nc")
land_mask_data = joinpath(mask_dataset_path(), "seamask.nc")

config_dict_atmos["output_dir"] = COUPLER_OUTPUT_DIR
atmos_sim = atmos_init(FT, config_dict_atmos);
thermo_params = get_thermo_params(atmos_sim) # TODO: this should be shared by all models

#=
We use a common `Space` for all global surfaces. This enables the MPI processes to operate on the same columns in both
the atmospheric and surface components, so exchanges are parallelized. Note this is only possible when the
atmosphere and surface are of the same horizontal resolution.
=#
## init a 2D boundary space at the surface
boundary_space = Spaces.horizontal_space(atmos_sim.domain.face_space)

# init land-sea fraction
land_fraction =
    FT.(
        Regridder.land_fraction(
            FT,
            REGRID_DIR,
            comms_ctx,
            land_mask_data,
            "LSMASK",
            boundary_space,
            mono = mono_surface,
        )
    )

@info mode_name
if mode_name == "amip"
    @info "AMIP boundary conditions - do not expect energy conservation"

    ## land
    land_sim = bucket_init(
        FT,
        tspan,
        config_dict["land_domain_type"],
        config_dict["land_albedo_type"],
        config_dict["land_temperature_anomaly"],
        comms_ctx,
        REGRID_DIR;
        dt = Δt_cpl,
        space = boundary_space,
        saveat = saveat,
        area_fraction = land_fraction,
        date_ref = date0,
        t_start = t_start,
    )

    ## ocean
    SST_info = bcfile_info_init(
        FT,
        REGRID_DIR,
        sst_data,
        "SST",
        boundary_space,
        comms_ctx,
        interpolate_daily = true,
        scaling_function = clean_sst, ## convert to Kelvin
        land_fraction = land_fraction,
        date0 = date0,
        mono = mono_surface,
    )

    update_midmonth_data!(date0, SST_info)
    SST_init = interpolate_midmonth_to_daily(date0, SST_info)
    ocean_sim = SurfaceStub((;
        T_sfc = SST_init,
        ρ_sfc = ClimaCore.Fields.zeros(boundary_space),
        z0m = FT(1e-3),
        z0b = FT(1e-3),
        beta = FT(1),
        α = FT(0.06),
        area_fraction = (FT(1) .- land_fraction),
        phase = TD.Liquid(),
        thermo_params = thermo_params,
    ))

    ## sea ice
    SIC_info = bcfile_info_init(
        FT,
        REGRID_DIR,
        sic_data,
        "SEAICE",
        boundary_space,
        comms_ctx,
        interpolate_daily = true,
        scaling_function = clean_sic, ## convert to fraction
        land_fraction = land_fraction,
        date0 = date0,
        mono = mono_surface,
    )
    update_midmonth_data!(date0, SIC_info)
    SIC_init = interpolate_midmonth_to_daily(date0, SIC_info)
    ice_fraction = get_ice_fraction.(SIC_init, mono_surface)
    ice_sim = ice_init(
        FT;
        tspan = tspan,
        dt = Δt_cpl,
        space = boundary_space,
        saveat = saveat,
        area_fraction = ice_fraction,
        thermo_params = thermo_params,
    )

    ## CO2 concentration
    CO2_info = bcfile_info_init(
        FT,
        REGRID_DIR,
        co2_data,
        "co2",
        boundary_space,
        comms_ctx,
        interpolate_daily = true,
        land_fraction = ones(boundary_space),
        date0 = date0,
        mono = mono_surface,
    )

    update_midmonth_data!(date0, CO2_info)
    CO2_init = interpolate_midmonth_to_daily(date0, CO2_info)
    update_field!(atmos_sim, Val(:co2_gm), CO2_init)

    mode_specifics = (; name = mode_name, SST_info = SST_info, SIC_info = SIC_info, CO2_info = CO2_info)

elseif mode_name in ("slabplanet", "slabplanet_aqua", "slabplanet_terra")

    land_fraction = mode_name == "slabplanet_aqua" ? land_fraction .* 0 : land_fraction
    land_fraction = mode_name == "slabplanet_terra" ? land_fraction .* 0 .+ 1 : land_fraction

    ## land
    land_sim = bucket_init(
        FT,
        tspan,
        config_dict["land_domain_type"],
        config_dict["land_albedo_type"],
        config_dict["land_temperature_anomaly"],
        comms_ctx,
        REGRID_DIR;
        dt = Δt_cpl,
        space = boundary_space,
        saveat = saveat,
        area_fraction = land_fraction,
        date_ref = date0,
        t_start = t_start,
    )

    ## ocean
    ocean_sim = ocean_init(
        FT;
        tspan = tspan,
        dt = Δt_cpl,
        space = boundary_space,
        saveat = saveat,
        area_fraction = (FT(1) .- land_fraction), ## NB: this ocean fraction includes areas covered by sea ice (unlike the one contained in the cs)
        thermo_params = thermo_params,
        evolving = evolving_ocean,
    )

    ## sea ice (here set to zero area coverage)
    ice_sim = SurfaceStub((;
        T_sfc = ClimaCore.Fields.ones(boundary_space),
        ρ_sfc = ClimaCore.Fields.zeros(boundary_space),
        z0m = FT(0),
        z0b = FT(0),
        beta = FT(1),
        α = FT(1),
        area_fraction = ClimaCore.Fields.zeros(boundary_space),
        phase = TD.Ice(),
        thermo_params = thermo_params,
    ))

    mode_specifics = (; name = mode_name, SST_info = nothing, SIC_info = nothing)
end

#=
## Coupler Initialization
The coupler needs to contain exchange information, manage the calendar and be able to access all component models. It can also optionally
save online diagnostics. These are all initialized here and saved in a global `CouplerSimulation` struct, `cs`.
=#

## coupler exchange fields
coupler_field_names = (
    :T_S,
    :z0m_S,
    :z0b_S,
    :ρ_sfc,
    :q_sfc,
    :albedo,
    :beta,
    :F_turb_energy,
    :F_turb_moisture,
    :F_turb_ρτxz,
    :F_turb_ρτyz,
    :F_radiative,
    :P_liq,
    :P_snow,
    :F_radiative_TOA,
    :P_net,
)
coupler_fields =
    NamedTuple{coupler_field_names}(ntuple(i -> ClimaCore.Fields.zeros(boundary_space), length(coupler_field_names)))

## model simulations
model_sims = (atmos_sim = atmos_sim, ice_sim = ice_sim, land_sim = land_sim, ocean_sim = ocean_sim);

## dates
dates = (; date = [date], date0 = [date0], date1 = [Dates.firstdayofmonth(date0)], new_month = [false])

#=
### Online Diagnostics
User can write custom diagnostics in the `user_diagnostics.jl`.
=#
monthly_3d_diags = init_diagnostics(
    (:T, :u, :q_tot, :q_liq_ice),
    atmos_sim.domain.center_space;
    save = Monthly(),
    operations = (; accumulate = TimeMean([Int(0)])),
    output_dir = COUPLER_OUTPUT_DIR,
    name_tag = "monthly_mean_3d_",
)

monthly_2d_diags = init_diagnostics(
    (:precipitation_rate, :toa_fluxes, :T_sfc, :tubulent_energy_fluxes),
    boundary_space;
    save = Monthly(),
    operations = (; accumulate = TimeMean([Int(0)])),
    output_dir = COUPLER_OUTPUT_DIR,
    name_tag = "monthly_mean_2d_",
)

diagnostics = (monthly_3d_diags, monthly_2d_diags)

#=
## Initialize Conservation Checks
=#
## init conservation info collector
conservation_checks = nothing
if energy_check
    @assert(
        mode_name[1:10] == "slabplanet" && !CA.is_distributed(ClimaComms.context(boundary_space)),
        "Only non-distributed slabplanet allowable for energy_check"
    )
    conservation_checks = (; energy = EnergyConservationCheck(model_sims), water = WaterConservationCheck(model_sims))
end

dir_paths = (; output = COUPLER_OUTPUT_DIR, artifacts = COUPLER_ARTIFACTS_DIR)
checkpoint_cb =
    HourlyCallback(dt = FT(480), func = checkpoint_sims, ref_date = [dates.date[1]], active = hourly_checkpoint) # 20 days
update_firstdayofmonth!_cb =
    MonthlyCallback(dt = FT(1), func = update_firstdayofmonth!, ref_date = [dates.date1[1]], active = true) # for BCReader
callbacks = (; checkpoint = checkpoint_cb, update_firstdayofmonth! = update_firstdayofmonth!_cb)

## coupler simulation
cs = CoupledSimulation{FT}(
    comms_ctx,
    dates,
    boundary_space,
    coupler_fields,
    config_dict,
    conservation_checks,
    [tspan[1], tspan[2]],
    atmos_sim.integrator.t,
    Δt_cpl,
    (; land = land_fraction, ocean = zeros(boundary_space), ice = zeros(boundary_space)),
    model_sims,
    mode_specifics,
    diagnostics,
    callbacks,
    dir_paths,
);

#=
## Restart component model states if specified
=#
#if restart_dir !== "unspecified"
#    for sim in cs.model_sims
#        if get_model_state_vector(sim) !== nothing
#            @skipping restart
#            restart_model_state!(sim, comms_ctx, restart_t; input_dir = restart_dir)
#        end
#    end
#end

#=
## Initialize Component Model Exchange
=#
turbulent_fluxes = nothing
if config_dict["turb_flux_partition"] == "PartitionedStateFluxes"
    turbulent_fluxes = PartitionedStateFluxes()
elseif config_dict["turb_flux_partition"] == "CombinedStateFluxes"
    turbulent_fluxes = CombinedStateFluxes()
else
    error("turb_flux_partition must be either PartitionedStateFluxes or CombinedStateFluxes")
end

# 1) coupler combines surface states and calculates rho_sfc using surface and atmos variables
update_surface_fractions!(cs)
import_combined_surface_fields!(cs.fields, cs.model_sims, cs.boundary_space, turbulent_fluxes)
import_atmos_fields!(cs.fields, cs.model_sims, cs.boundary_space, turbulent_fluxes)
update_model_sims!(cs.model_sims, cs.fields, turbulent_fluxes)

# 2) each surface component model calculates its own vapor specific humidity (q_sfc)
# TODO: the q_sfc calculation follows the design of the bucket q_sfc, but it would be neater to abstract this from step!
step!(land_sim, Δt_cpl)
step!(ocean_sim, Δt_cpl)
step!(ice_sim, Δt_cpl)

# 3) coupler re-imports updated surface fields and calculates turbulent fluxes, while updating atmos sfc_conditions
if turbulent_fluxes isa CombinedStateFluxes
    # calculate fluxes using combined surface states on the atmos grid
    import_combined_surface_fields!(cs.fields, cs.model_sims, cs.boundary_space, turbulent_fluxes) # i.e. T_sfc, albedo, z0, beta, q_sfc
    combined_turbulent_fluxes!(cs.model_sims, cs.fields, turbulent_fluxes) # this updates the atmos thermo state, sfc_ts
elseif turbulent_fluxes isa PartitionedStateFluxes
    # calculate turbulent fluxes in surface models and save the weighted average in coupler fields
    partitioned_turbulent_fluxes!(cs.model_sims, cs.fields, cs.boundary_space, MoninObukhovScheme(), thermo_params)

    # update atmos sfc_conditions for surface temperature
    # TODO: this is hard coded and needs to be simplified (need CA modification)
    new_p = get_new_cache(atmos_sim, cs.fields)
    CA.SurfaceConditions.update_surface_conditions!(atmos_sim.integrator.u, new_p, atmos_sim.integrator.t) # sets T_sfc (but SF calculation not necessary - CA)
    atmos_sim.integrator.p.precomputed.sfc_conditions .= new_p.precomputed.sfc_conditions
end

# 4) given the new sfc_conditions, atmos calls the radiative flux callback
reinit_model_sims!(cs.model_sims) # NB: for atmos this sets a nonzero radiation flux

# 5) coupler re-imports updated atmos fluxes (radiative fluxes for both `turbulent_fluxes` types
# and also turbulent fluxes if `turbulent_fluxes isa CombinedStateFluxes`,
# and sends them to the surface component models. If `turbulent_fluxes isa PartitionedStateFluxes`
# atmos receives the turbulent fluxes from the coupler.
import_atmos_fields!(cs.fields, cs.model_sims, cs.boundary_space, turbulent_fluxes)
update_model_sims!(cs.model_sims, cs.fields, turbulent_fluxes)

#=
## Coupling Loop
=#
function solve_coupler!(cs::ClimaCoupler.Interfacer.CoupledSimulation)
    @info "Starting coupling loop"

    (; model_sims, Δt_cpl, tspan) = cs
    (; atmos_sim, land_sim, ocean_sim, ice_sim) = model_sims

    ## step in time
    walltime = @elapsed for t in ((tspan[1] + Δt_cpl):Δt_cpl:tspan[end])

        cs.dates.date[1] = current_date(cs, t) # if not global, `date` is not updated.

        ## print date on the first of month
        if cs.dates.date[1] >= cs.dates.date1[1]
            ClimaComms.iamroot(comms_ctx) ? @show(cs.dates.date[1]) : nothing
        end

        if cs.mode.name == "amip"

            ## monthly read of boundary condition data for SST and SIC and CO2
            if cs.dates.date[1] >= next_date_in_file(cs.mode.SST_info)
                update_midmonth_data!(cs.dates.date[1], cs.mode.SST_info)
            end
            SST_current = interpolate_midmonth_to_daily(cs.dates.date[1], cs.mode.SST_info)
            update_field!(ocean_sim, Val(:surface_temperature), SST_current)

            if cs.dates.date[1] >= next_date_in_file(cs.mode.SIC_info)
                update_midmonth_data!(cs.dates.date[1], cs.mode.SIC_info)
            end
            SIC_current =
                get_ice_fraction.(interpolate_midmonth_to_daily(cs.dates.date[1], cs.mode.SIC_info), mono_surface)
            update_field!(ice_sim, Val(:area_fraction), SIC_current)

            if cs.dates.date[1] >= next_date_in_file(cs.mode.CO2_info)
                update_midmonth_data!(cs.dates.date[1], cs.mode.CO2_info)
            end
            CO2_current = interpolate_midmonth_to_daily(cs.dates.date[1], cs.mode.CO2_info)
            update_field!(atmos_sim, Val(:co2_gm), CO2_current)

            ## calculate and accumulate diagnostics at each timestep
            ClimaComms.barrier(comms_ctx)
            accumulate_diagnostics!(cs)

            ## save and reset monthly averages
            save_diagnostics(cs)

        end

        ## compute global energy
        !isnothing(cs.conservation_checks) ? check_conservation!(cs) : nothing

        ## run component models sequentially for one coupling timestep (Δt_cpl)
        ClimaComms.barrier(comms_ctx)
        update_surface_fractions!(cs)
        update_model_sims!(cs.model_sims, cs.fields, turbulent_fluxes)

        ## step sims
        step_model_sims!(cs.model_sims, t)

        ## exchange combined fields and (if specified) calculate fluxes using combined states
        import_combined_surface_fields!(cs.fields, cs.model_sims, cs.boundary_space, turbulent_fluxes) # i.e. T_sfc, albedo, z0, beta
        if turbulent_fluxes isa CombinedStateFluxes
            combined_turbulent_fluxes!(cs.model_sims, cs.fields, turbulent_fluxes) # this updates the surface thermo state, sfc_ts, in ClimaAtmos (but also unnecessarily calculates fluxes)
        elseif turbulent_fluxes isa PartitionedStateFluxes
            # calculate turbulent fluxes in surfaces and save the weighted average in coupler fields
            partitioned_turbulent_fluxes!(cs.model_sims, cs.fields, cs.boundary_space, MoninObukhovScheme(), thermo_params)

            # update atmos sfc_conditions for surface temperature - TODO: this needs to be simplified (need CA modification)
            new_p = get_new_cache(atmos_sim, cs.fields)
            CA.SurfaceConditions.update_surface_conditions!(atmos_sim.integrator.u, new_p, atmos_sim.integrator.t) # to set T_sfc (but SF calculation not necessary - CA modification)
            atmos_sim.integrator.p.precomputed.sfc_conditions .= new_p.precomputed.sfc_conditions
        end

        import_atmos_fields!(cs.fields, cs.model_sims, cs.boundary_space, turbulent_fluxes) # radiative and/or turbulent

        ## callback to update the fist day of month if needed (for BCReader)
        trigger_callback!(cs, cs.callbacks.update_firstdayofmonth!)

        ## callback to checkpoint model state
        trigger_callback!(cs, cs.callbacks.checkpoint)

    end
    @show walltime

    return cs
end
