# # Dry Held-Suarez
## This script runs an idealized global circulation model, as in Held and Suarez (1994).
## There is no moisture, the radiation is approximated by a Newtonian cooling scheme,
## and dissipation is applied via a Rayleigh damping scheme.
## The numerics follow ClimaCore's dry Held-Suarez `held_suarez_rhoe` example.

redirect_stderr(IOContext(stderr, :stacktrace_types_limited => Ref(false)))

#=
## Configuration initialization
=#

#=
### Package Import
=#

## standard packages
using Dates

## ClimaESM packages
using ClimaComms
@static pkgversion(ClimaComms) >= v"0.6" && ClimaComms.@import_required_backends
import ClimaAtmos as CA
import ClimaCore

## Coupler specific imports
import ClimaCoupler
import ClimaCoupler: Checkpointer, FieldExchanger, Interfacer, TimeManager, Utilities

pkg_dir = pkgdir(ClimaCoupler)

#=
### Helper Functions
=#

## helpers for component models
include("components/atmosphere/climaatmos.jl")

#=
### Setup simulation parameters
Here we follow ClimaCore's dry Held-Suarez `held_suarez_rhoe` example.
=#

## run names
job_id = "dry_held_suarez"
coupler_output_dir = "$job_id"
const FT = Float64
restart_dir = nothing
restart_t = Int(0)

## coupler simulation specific configuration
Δt_cpl = Float64(400)
t_end = "1000days"
tspan = (Float64(0.0), Float64(Utilities.time_to_seconds(t_end)))
start_date = "19790301"
hourly_checkpoint = true

#=
### I/O Directory Setup
=#

dir_paths = Utilities.setup_output_dirs(output_dir = coupler_output_dir, comms_ctx = ClimaComms.context())
@info(config_dict)

## namelist
config_dict = Dict(
    # general
    "FLOAT_TYPE" => string(FT),
    # file paths
    "atmos_config_file" => nothing,
    "coupler_toml_file" => nothing,
    "coupler_output_dir" => coupler_output_dir,
    "mode_name" => "",
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
    "vert_diff" => "false",
    "hyperdiff" => "CAM_SE",
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
    # held-suarez specific
    "forcing" => "held_suarez",
)

## merge dictionaries of command line arguments, coupler dictionary and component model dictionaries
atmos_output_dir = joinpath(dir_paths.output, "clima_atmos")
atmos_config_dict = get_atmos_config_dict(config_dict, job_id, atmos_output_dir)
atmos_config_object = CA.AtmosConfig(atmos_config_dict)

#=
## Setup Communication Context
=#

comms_ctx = Utilities.get_comms_context(Dict("device" => "auto"))

#=
## Component Model Initialization
=#

## init atmos model component
atmos_sim = atmos_init(atmos_config_object);
thermo_params = get_thermo_params(atmos_sim)

#=
### Boundary Space
=#

## init a 2D boundary space at the surface
boundary_space = ClimaCore.Spaces.horizontal_space(atmos_sim.domain.face_space)

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
    NamedTuple{coupler_field_names}(ntuple(i -> ClimaCore.Fields.zeros(boundary_space), length(coupler_field_names)))

## model simulations
model_sims = (atmos_sim = atmos_sim,);

## dates
date0 = date = Dates.DateTime(start_date, Dates.dateformat"yyyymmdd")
dates = (; date = [date], date0 = [date0], date1 = [Dates.firstdayofmonth(date0)], new_month = [false])

#=
## Initialize Callbacks
=#
checkpoint_cb = TimeManager.HourlyCallback(
    dt = FT(480),
    func = Checkpointer.checkpoint_sims,
    ref_date = [dates.date[1]],
    active = hourly_checkpoint,
) # 20 days TODO: not GPU friendly
update_firstdayofmonth!_cb = TimeManager.MonthlyCallback(
    dt = FT(1),
    func = TimeManager.update_firstdayofmonth!,
    ref_date = [dates.date1[1]],
    active = true,
)
callbacks = (; checkpoint = checkpoint_cb, update_firstdayofmonth! = update_firstdayofmonth!_cb)

cs = Interfacer.CoupledSimulation{FT}(
    comms_ctx,
    dates,
    boundary_space,
    coupler_fields,
    nothing, # conservation checks
    [tspan[1], tspan[2]],
    Δt_cpl,
    model_sims,
    (;), # mode_specifics
    callbacks,
    dir_paths,
    nothing, # turbulent_fluxes
    thermo_params,
    nothing, # amip_diags_handler
);

#=
## Restart component model states if specified in the config_dict
=#

if !isnothing(restart_dir)
    for sim in cs.model_sims
        if get_model_prog_state(sim) !== nothing
            Checkpointer.restart_model_state!(sim, comms_ctx, restart_t; input_dir = restart_dir)
        end
    end
end

#=
## Coupling Loop
=#

function solve_coupler!(cs)
    (; model_sims, Δt_cpl, tspan, comms_ctx) = cs
    (; atmos_sim) = model_sims

    @info("Starting coupling loop")

    ## step in time
    walltime = @elapsed for t in ((tspan[begin] + Δt_cpl):Δt_cpl:tspan[end])

        cs.dates.date[1] = TimeManager.current_date(cs, t)

        ## print date on the first of month
        cs.dates.date[1] >= cs.dates.date1[1] && @info(cs.dates.date[1])

        ## step sims
        FieldExchanger.step_model_sims!(cs.model_sims, t)

        FieldExchanger.import_atmos_fields!(cs.fields, cs.model_sims, cs.boundary_space, cs.turbulent_fluxes) # radiative and/or turbulent

        ## callback to update the fist day of month if needed
        TimeManager.trigger_callback!(cs, cs.callbacks.update_firstdayofmonth!)

        ## callback to checkpoint model state
        TimeManager.trigger_callback!(cs, cs.callbacks.checkpoint)

    end
    @info(walltime)

    return cs
end

## run the coupled simulation
solve_coupler!(cs);
