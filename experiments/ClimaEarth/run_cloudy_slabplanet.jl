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
import YAML

# ## ClimaESM packages
import ClimaComms
@static pkgversion(ClimaComms) >= v"0.6" && ClimaComms.@import_required_backends
import ClimaAtmos as CA
import ClimaCore as CC

# ## Coupler specific imports
import ClimaCoupler
import ClimaCoupler:
    ConservationChecker,
    Checkpointer,
    Diagnostics,
    FieldExchanger,
    FluxCalculator,
    Interfacer,
    Regridder,
    TimeManager,
    Utilities

pkg_dir = pkgdir(ClimaCoupler)

#=
### Helper Functions
=#

## helpers for component models
include("components/atmosphere/climaatmos.jl")
include("components/ocean/slab_ocean.jl")
include("components/land/climaland_bucket.jl")

## helpers for user-specified IO
include("user_io/user_diagnostics.jl")
include("user_io/user_logging.jl")
include("user_io/io_helpers.jl")

#=
### Setup simulation parameters
=#

## run names
job_id = "cloudy_slabplanet"
coupler_output_dir = "$job_id"
const FT = Float64
restart_dir = "unspecified"
restart_t = Int(0)

## coupler simulation specific configuration
Δt_cpl = Float64(100)
t_end = "1000days"
tspan = (Float64(0.0), Float64(time_to_seconds(t_end)))
start_date = "19790321"
hourly_checkpoint = true
dt_rad = "6hours"

## namelist
config_dict = Dict(
    # general
    "FLOAT_TYPE" => string(FT),
    # file paths
    "atmos_config_file" => nothing,
    "coupler_toml_file" => nothing,
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
    "diagnostics" => [
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
    "override_τ_precip" => false,
)

## merge dictionaries of command line arguments, coupler dictionary and component model dictionaries
atmos_config_dict, config_dict = get_atmos_config_dict(config_dict, job_id)
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
## Setup Communication Context
=#

comms_ctx = Utilities.get_comms_context(Dict("device" => "auto"))

#=
### I/O Directory Setup
=#

dir_paths = setup_output_dirs(output_dir = coupler_output_dir, comms_ctx = comms_ctx)
@info(config_dict)

#=
## Data File Paths
=#
include(joinpath(pkgdir(ClimaCoupler), "artifacts", "artifact_funcs.jl"))
land_mask_data = artifact_data(mask_dataset_path(), "seamask")

#=
## Component Model Initialization
=#

## dates
date0 = date = Dates.DateTime(start_date, Dates.dateformat"yyyymmdd")
dates = (; date = [date], date0 = [date0], date1 = [Dates.firstdayofmonth(date0)], new_month = [false])


#=
### Atmosphere
This uses the `ClimaAtmos.jl` model, with parameterization options specified in the `config_dict_atmos` dictionary.
=#

## init atmos model component
atmos_sim = atmos_init(atmos_config_object);
surface_elevation = CC.Fields.level(CC.Fields.coordinate_field(atmos_sim.integrator.u.f).z, CC.Utilities.half)
thermo_params = get_thermo_params(atmos_sim)

#=
### Boundary Space
=#

## init a 2D boundary space at the surface
boundary_space = CC.Spaces.horizontal_space(atmos_sim.domain.face_space) # TODO: specify this in the coupler and pass it to all component models #665

#=
### Land-sea Fraction
This is a static field that contains the area fraction of land and sea, ranging from 0 to 1. If applicable, sea ice is included in the sea fraction. at this stage.
Note that land-sea area fraction is different to the land-sea mask, which is a binary field (masks are used internally by the coupler to indicate passive cells that are not populated by a given component model).
=#

land_area_fraction =
    FT.(
        Regridder.land_fraction(
            FT,
            dir_paths.regrid,
            comms_ctx,
            land_mask_data,
            "LSMASK",
            boundary_space,
            mono = false,
        )
    )

#=
### Surface Model: Bucket Land and Slab Ocean
=#

saveat = Float64(time_to_seconds(config_dict["dt_save_to_sol"]))

## land model
land_sim = bucket_init(
    FT,
    tspan,
    "sphere",
    "map_static",
    "aquaplanet",
    dir_paths.regrid;
    dt = Δt_cpl,
    space = boundary_space,
    saveat = saveat,
    area_fraction = land_area_fraction,
    date_ref = date0,
    t_start = tspan[1],
    energy_check = false,
    surface_elevation,
)

ocean_sim = ocean_init(
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
coupler_fields =
    NamedTuple{coupler_field_names}(ntuple(i -> CC.Fields.zeros(boundary_space), length(coupler_field_names)))
Utilities.show_memory_usage()

## model simulations
model_sims = (atmos_sim = atmos_sim, ocean_sim = ocean_sim);

#=
## Initialize Callbacks
=#

checkpoint_cb = TimeManager.HourlyCallback(
    dt = FT(480),
    func = checkpoint_sims,
    ref_date = [dates.date[1]],
    active = hourly_checkpoint,
) # 20 days
update_firstdayofmonth!_cb = TimeManager.MonthlyCallback(
    dt = FT(1),
    func = TimeManager.update_firstdayofmonth!,
    ref_date = [dates.date1[1]],
    active = true,
)
dt_water_albedo = parse(FT, filter(x -> !occursin(x, "hours"), dt_rad))
albedo_cb = TimeManager.HourlyCallback(
    dt = dt_water_albedo,
    func = FluxCalculator.water_albedo_from_atmosphere!,
    ref_date = [dates.date[1]],
    active = true,
)
callbacks =
    (; checkpoint = checkpoint_cb, update_firstdayofmonth! = update_firstdayofmonth!_cb, water_albedo = albedo_cb)

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
    dates,
    boundary_space,
    coupler_fields,
    config_dict,
    nothing, # conservation checks
    [tspan[1], tspan[2]],
    atmos_sim.integrator.t,
    Δt_cpl,
    (; land = zeros(boundary_space), ocean = ones(boundary_space), ice = zeros(boundary_space)),
    model_sims,
    (;), # mode_specifics
    (), # coupler diagnostics
    callbacks,
    dir_paths,
    turbulent_fluxes,
    thermo_params,
);

#=
## Restart component model states if specified in the config_dict
=#

if restart_dir !== "unspecified"
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

        cs.dates.date[1] = TimeManager.current_date(cs, t)

        ## print date on the first of month
        cs.dates.date[1] >= cs.dates.date1[1] && @info(cs.dates.date[1])

        ClimaComms.barrier(comms_ctx)

        ## update water albedo from wind at dt_water_albedo (this will be extended to a radiation callback from the coupler)
        TimeManager.trigger_callback!(cs, cs.callbacks.water_albedo)

        ## run component models sequentially for one coupling timestep (Δt_cpl)
        FieldExchanger.update_model_sims!(cs.model_sims, cs.fields, cs.turbulent_fluxes)

        ## step sims
        FieldExchanger.step_model_sims!(cs.model_sims, t)

        ## exchange combined fields and (if specified) calculate fluxes using combined states
        FieldExchanger.import_combined_surface_fields!(cs.fields, cs.model_sims, cs.turbulent_fluxes) # i.e. T_sfc, surface_albedo, z0, beta
        FluxCalculator.combined_turbulent_fluxes!(cs.model_sims, cs.fields, cs.turbulent_fluxes)

        FieldExchanger.import_atmos_fields!(cs.fields, cs.model_sims, cs.boundary_space, cs.turbulent_fluxes) # radiative and/or turbulent

        ## callback to update the fist day of month if needed
        TimeManager.trigger_callback!(cs, cs.callbacks.update_firstdayofmonth!)

        ## callback to checkpoint model state
        TimeManager.trigger_callback!(cs, cs.callbacks.checkpoint)

    end

    return nothing
end

solve_coupler!(cs)
