# # Cloudy Slabplanet
## This follows the setup of the Cloudy Aquaplanet experiment, but with a slab ocean and
## bucket land model for the surface, partitioned using a land-sea mask.

redirect_stderr(IOContext(stderr, :stacktrace_types_limited => Ref(false)))

#=
## Configuration initialization
=#

#=
### Package Import
=#

## standard packages
import Dates

# ## ClimaESM packages
import ClimaComms
@static pkgversion(ClimaComms) >= v"0.6" && ClimaComms.@import_required_backends
import ClimaAtmos as CA
import ClimaCore as CC

# ## Coupler specific imports
import ClimaCoupler
import ClimaCoupler:
    ConservationChecker, Checkpointer, FieldExchanger, FluxCalculator, Interfacer, TimeManager, Utilities

import ClimaUtilities.ClimaArtifacts: @clima_artifact
import ClimaUtilities.SpaceVaryingInputs: SpaceVaryingInput
import Interpolations # triggers InterpolationsExt in ClimaUtilities

# TODO: Move to ClimaUtilities once we move the Schedules to ClimaUtilities
import ClimaDiagnostics.Schedules: EveryCalendarDtSchedule
import ClimaUtilities.TimeManager: ITime

pkg_dir = pkgdir(ClimaCoupler)

#=
## Setup Communication Context
=#

comms_ctx = Utilities.get_comms_context(Dict("device" => "auto"))

#=
### Helper Functions
=#

## helpers for component models
include("components/atmosphere/climaatmos.jl")
include("components/ocean/slab_ocean.jl")
include("components/land/climaland_bucket.jl")

#=
### Setup simulation parameters
=#

## run names
job_id = "cloudy_slabplanet"
coupler_output_dir = "$job_id"
const FT = Float64
restart_dir = nothing
restart_t = Int(0)

## coupler simulation specific configuration
Δt_cpl = Float64(100)
t_end = "1000days"
tspan = (Float64(0.0), Float64(Utilities.time_to_seconds(t_end)))
start_date = "19790321"
use_itime = true
if use_itime
    tspan = (
        ITime(0.0, epoch = Dates.DateTime(1979, 3, 1)),
        ITime(Utilities.time_to_seconds(t_end), epoch = Dates.DateTime(1979, 3, 1)),
    )
end
checkpoint_dt = "20days"
dt_rad = "6hours"

#=
### I/O Directory Setup
=#

dir_paths = Utilities.setup_output_dirs(output_dir = coupler_output_dir, comms_ctx = ClimaComms.context())


## namelist
config_dict = Dict(
    # general
    "FLOAT_TYPE" => string(FT),
    # file paths
    "atmos_config_file" => nothing,
    "coupler_toml" => [],
    "coupler_output_dir" => coupler_output_dir,
    "mode_name" => "",
    "job_id" => job_id,
    "atmos_config_repo" => "ClimaAtmos",
    # timestepping
    "dt" => "$(Δt_cpl)secs",
    "dt_save_to_sol" => "1days",
    "t_end" => t_end,
    "start_date" => "19790301",
    # domain
    "h_elem" => 4,
    "z_elem" => 10,
    "z_max" => 30000.0, # semi-high top
    "dz_bottom" => 300.0,
    "nh_poly" => 4,
    # output
    "dt_save_to_sol" => "1days",
    "checkpoint_dt" => "1days",
    # numerics
    "apply_limiter" => false,
    "viscous_sponge" => false,
    "rayleigh_sponge" => false,
    # "vert_diff" => "true", #required
    "hyperdiff" => "CAM_SE",
    "ode_algo" => "ARS343",
    # run
    "surface_setup" => "PrescribedSurface",
    # diagnostic (nested with period and short_name)
    "output_default_diagnostics" => false,
    "extra_atmos_diagnostics" => [
        Dict(
            "short_name" =>
                ["mse", "lr", "mass_strf", "stab", "vt", "egr", "ua", "va", "wa", "ta", "rhoa", "pfull"],
            "period" => "6hours",
            "reduction" => "inst",
        ),
    ],
    # cloudy aquaplanet specific
    "precip_model" => "0M",
    "moist" => "equil",
    "prognostic_surface" => "PrescribedSurfaceTemperature",
    "turb_flux_partition" => "CombinedStateFluxesMOST",
    "rad" => "allskywithclear",
    "idealized_insolation" => true, # perpetual equinox with no diurnal cycle
    "dt_rad" => dt_rad,
    "turbconv" => "diagnostic_edmfx",
    "dt_cloud_fraction" => "1hours",
    "implicit_diffusion" => true,
    "approximate_linear_solve_iters" => 2,
    "prognostic_tke" => true,
    "edmfx_upwinding" => "first_order",
    "edmfx_entr_model" => "Generalized",
    "edmfx_detr_model" => "Generalized",
    "edmfx_nh_pressure" => true,
    "edmfx_sgs_mass_flux" => true,
    "edmfx_sgs_diffusive_flux" => true,
    "override_precip_timescale" => false,
    "use_itime" => use_itime,
)

atmos_output_dir = joinpath(dir_paths.output, "clima_atmos")
land_output_dir = joinpath(dir_paths.output, "clima_land")

## merge dictionaries of command line arguments, coupler dictionary and component model dictionaries
atmos_config_dict = get_atmos_config_dict(config_dict, job_id, atmos_output_dir)
atmos_config_object = CA.AtmosConfig(atmos_config_dict)

# override default toml parameters
atmos_config_object.toml_dict["precipitation_timescale"]["value"] = 600
atmos_config_object.toml_dict["entr_inv_tau"]["value"] = 0.002
atmos_config_object.toml_dict["entr_coeff"]["value"] = 0
atmos_config_object.toml_dict["detr_inv_tau"]["value"] = 0
atmos_config_object.toml_dict["detr_vertdiv_coeff"]["value"] = 0.6
atmos_config_object.toml_dict["detr_buoy_coeff"]["value"] = 0.12
atmos_config_object.toml_dict["min_area_limiter_scale"]["value"] = 0
atmos_config_object.toml_dict["max_area_limiter_scale"]["value"] = 0

#=
## Component Model Initialization
=#

## start date
start_date = Dates.DateTime(start_date, Dates.dateformat"yyyymmdd")

#=
### Atmosphere
This uses the `ClimaAtmos.jl` model, with parameterization options specified in the `config_dict_atmos` dictionary.
=#

## init atmos model component
atmos_sim = ClimaAtmosSimulation(atmos_config_object);
surface_elevation = CC.Fields.level(CC.Fields.coordinate_field(atmos_sim.integrator.u.f).z, CC.Utilities.half)
thermo_params = get_thermo_params(atmos_sim)

#=
### Boundary Space
=#

## init a 2D boundary space at the surface
boundary_space = CC.Spaces.horizontal_space(atmos_sim.domain.face_space) # TODO: specify this in the coupler and pass it to all component models #665

# Land initial condition
# Use the default land initial condition (not reading from a file)
land_initial_condition = ""

#=
### Land-sea Fraction
This is a static field that contains the area fraction of land and sea, ranging from 0 to 1. If applicable, sea ice is included in the sea fraction. at this stage.
Note that land-sea area fraction is different to the land-sea mask, which is a binary field (masks are used internally by the coupler to indicate passive cells that are not populated by a given component model).
=#
land_mask_data = joinpath(@clima_artifact("landsea_mask_60arcseconds", comms_ctx), "landsea_mask.nc")
land_area_fraction = SpaceVaryingInput(land_mask_data, "landsea", boundary_space)
#=
### Surface Model: Bucket Land and Slab Ocean
=#

saveat = Float64(Utilities.time_to_seconds(config_dict["dt_save_to_sol"]))
if use_itime
    saveat = ITime(saveat)
    Δt_cpl, t0, tf = promote(ITime(Δt_cpl), tspan[1], tspan[2])
    tspan = (t0, tf)
end
saveat = [tspan[1]:saveat:tspan[1]..., tspan[2]]

## land model
land_sim = BucketSimulation(
    FT;
    dt = Δt_cpl,
    tspan,
    start_date,
    output_dir = land_output_dir,
    boundary_space,
    area_fraction = land_area_fraction,
    saveat,
    surface_elevation,
    land_temperature_anomaly = "aquaplanet",
    land_initial_condition,
)

ocean_sim = SlabOceanSimulation(
    FT;
    tspan = tspan,
    dt = Δt_cpl,
    space = boundary_space,
    saveat = saveat,
    area_fraction = ones(boundary_space),
    thermo_params = thermo_params,
    evolving = true,
)

#=
## Coupler Initialization
=#

## coupler exchange fields
coupler_field_names = [
    :T_sfc,
    :z0m_sfc,
    :z0b_sfc,
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
]
coupler_fields = Interfacer.init_coupler_fields(FT, coupler_field_names, boundary_space)
Utilities.show_memory_usage()

## model simulations
model_sims = (atmos_sim = atmos_sim, ocean_sim = ocean_sim);

#=
## Initialize Callbacks
=#
schedule_checkpoint = EveryCalendarDtSchedule(TimeManager.time_to_period(checkpoint_dt); start_date)
checkpoint_cb = TimeManager.Callback(schedule_checkpoint, Checkpointer.checkpoint_sims)

schedule_albedo = EveryCalendarDtSchedule(TimeManager.time_to_period(dt_rad); start_date)
albedo_cb = TimeManager.Callback(schedule_albedo, FluxCalculator.water_albedo_from_atmosphere!)

callbacks = (; checkpoint = checkpoint_cb, water_albedo = albedo_cb)

#=
## Initialize turbulent fluxes
=#
turbulent_fluxes = nothing
if config_dict["turb_flux_partition"] == "CombinedStateFluxesMOST"
    turbulent_fluxes = FluxCalculator.CombinedStateFluxesMOST()
else
    error("turb_flux_partition must be CombinedStateFluxesMOST")
end

#=
## Initialize Coupled Simulation
=#

cs = Interfacer.CoupledSimulation{FT}(
    comms_ctx,
    Ref(start_date),
    boundary_space,
    coupler_fields,
    nothing, # conservation checks
    [tspan[1], tspan[2]],
    Δt_cpl,
    Ref(tspan[1]),
    model_sims,
    callbacks,
    dir_paths,
    turbulent_fluxes,
    thermo_params,
    nothing, # diags_handler
);

#=
## Restart component model states if specified in the config_dict
=#

if !isnothing(restart_dir)
    for sim in cs.model_sims
        if Checkpointer.get_model_prog_state(sim) !== nothing
            Checkpointer.restart_model_state!(sim, comms_ctx, restart_t; input_dir = restart_dir)
        end
    end
end

#=
## Initialize Component Model Exchange
=#

# 1.surface density (`ρ_sfc`): calculated by the coupler by adiabatically extrapolating atmospheric thermal state to the surface.
# For this, we need to import surface and atmospheric fields. The model sims are then updated with the new surface density.
FieldExchanger.import_combined_surface_fields!(cs.fields, cs.model_sims, cs.turbulent_fluxes)
FieldExchanger.import_atmos_fields!(cs.fields, cs.model_sims, cs.boundary_space, cs.turbulent_fluxes)
FieldExchanger.update_model_sims!(cs.model_sims, cs.fields, cs.turbulent_fluxes)

# 2.surface vapor specific humidity (`q_sfc`): step surface models with the new surface density to calculate their respective `q_sfc` internally
Interfacer.step!(ocean_sim, Δt_cpl)

# 3.turbulent fluxes
## import the new surface properties into the coupler (note the atmos state was also imported in step 3.)
FieldExchanger.import_combined_surface_fields!(cs.fields, cs.model_sims, cs.turbulent_fluxes) # i.e. T_sfc, albedo, z0, beta, q_sfc
## calculate turbulent fluxes inside the atmos cache based on the combined surface state in each grid box
FluxCalculator.combined_turbulent_fluxes!(cs.model_sims, cs.fields, cs.turbulent_fluxes) # this updates the atmos thermo state, sfc_ts

# 4.reinitialize models + radiative flux: prognostic states and time are set to their initial conditions.
FieldExchanger.reinit_model_sims!(cs.model_sims)

# 5.update all fluxes: coupler re-imports updated atmos fluxes
FieldExchanger.import_atmos_fields!(cs.fields, cs.model_sims, cs.boundary_space, cs.turbulent_fluxes)
FieldExchanger.update_model_sims!(cs.model_sims, cs.fields, cs.turbulent_fluxes)

#=
## Coupling Loop
=#

function solve_coupler!(cs)
    (; Δt_cpl, tspan, comms_ctx) = cs

    @info("Starting coupling loop")
    ## step in time
    for t in ((tspan[begin] + Δt_cpl):Δt_cpl:tspan[end])
        # Update current time
        cs.t[] = t

        ClimaComms.barrier(comms_ctx)

        ## update water albedo from wind at dt_water_albedo (this will be extended to a radiation callback from the coupler)
        TimeManager.maybe_trigger_callback(cs.callbacks.water_albedo, cs)

        ## run component models sequentially for one coupling timestep (Δt_cpl)
        FieldExchanger.update_model_sims!(cs.model_sims, cs.fields, cs.turbulent_fluxes)

        ## step sims
        FieldExchanger.step_model_sims!(cs.model_sims, t)

        ## exchange combined fields and (if specified) calculate fluxes using combined states
        FieldExchanger.import_combined_surface_fields!(cs.fields, cs.model_sims, cs.turbulent_fluxes) # i.e. T_sfc, surface_albedo, z0, beta
        FluxCalculator.combined_turbulent_fluxes!(cs.model_sims, cs.fields, cs.turbulent_fluxes)

        FieldExchanger.import_atmos_fields!(cs.fields, cs.model_sims, cs.boundary_space, cs.turbulent_fluxes) # radiative and/or turbulent

        ## callback to checkpoint model state
        TimeManager.maybe_trigger_callback(cs.callbacks.checkpoint, cs)
    end

    return nothing
end

solve_coupler!(cs)
