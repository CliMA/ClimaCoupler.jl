# # Dry Held-Suarez

redirect_stderr(IOContext(stderr, :stacktrace_types_limited => Ref(false)))

#=
### Package Import
=#

## standard packages
using Dates
import YAML

## ClimaESM packages
import ClimaAtmos as CA
using ClimaCore

## Coupler specific imports
using ClimaCoupler
using ClimaCoupler.Checkpointer: restart_model_state!
using ClimaCoupler.FieldExchanger: import_atmos_fields!, step_model_sims!
using ClimaCoupler.Interfacer: CoupledSimulation
using ClimaCoupler.TimeManager:
    current_date, HourlyCallback, MonthlyCallback, update_firstdayofmonth!, trigger_callback!
import ClimaCoupler.Utilities: get_comms_context

pkg_dir = pkgdir(ClimaCoupler)

#=
### Helper Functions
=#

## helpers for component models
include("components/atmosphere/climaatmos.jl")

## helpers for user-specified IO
include("user_io/user_diagnostics.jl")
include("user_io/user_logging.jl")

include("user_io/io_helpers.jl")

#=
### Setup simulation parameters
Here we follow ClimaCore's dry Held-Suarez `held_suarez_rhoe` example.
=#

## run names
run_name = "dry_held_suarez"
coupler_output_dir = "$run_name"
const FT = Float64
restart_dir = "unspecified"
restart_t = Int(0)

## coupler simulation specific configuration
Δt_cpl = Float64(400)
t_end = "1000days"
tspan = (Float64(0.0), Float64(time_to_seconds(t_end)))
start_date = "19790301"
hourly_checkpoint = true

## namelist
config_dict = Dict(
    # file paths
    "atmos_config_file" => nothing,
    "coupler_toml_file" => nothing,
    "coupler_output_dir" => coupler_output_dir,
    "mode_name" => "",
    "run_name" => run_name,
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
    "hyperdiff" => "ClimaHyperdiffusion",
    # run
    "job_id" => run_name,
    "surface_setup" => "PrescribedSurface",
    # diagnostic (nested with period and short_name)
    "output_default_diagnostics" => false,
    "diagnostics" => [
        Dict(
            "short_name" => [
                "mse",
                "lr",
                "mass_strf",
                "stab",
                "vt",
                "egr",
                "ua",
                "va",
                "wa",
                "ta",
                "rhoa",
                "pfull",
                "stab",
            ],
            "period" => "6hours",
            "reduction" => "inst",
        ),
    ],
    # held-suarez specific
    "forcing" => "held_suarez",
)

## merge dictionaries of command line arguments, coupler dictionary and component model dictionaries
config_dict_atmos = get_atmos_config_dict(config_dict)
config_dict = merge(config_dict_atmos, config_dict) # config_dict overrides config_dict_atmos
atmos_config_object = CA.AtmosConfig(config_dict_atmos)

#=
## Setup Communication Context
=#

using ClimaComms
comms_ctx = get_comms_context(Dict("device" => "auto"))
ClimaComms.init(comms_ctx)

#=
### I/O Directory Setup
=#

dir_paths = setup_output_dirs(output_dir = coupler_output_dir, comms_ctx = comms_ctx)
ClimaComms.iamroot(comms_ctx) ? @info(config_dict) : nothing

#=
## Component Model Initialization
=#

## init atmos model component
atmos_sim = atmos_init(FT, atmos_config_object);
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
date0 = date = DateTime(start_date, dateformat"yyyymmdd")
dates = (; date = [date], date0 = [date0], date1 = [Dates.firstdayofmonth(date0)], new_month = [false])

#=
## Initialize Callbacks
=#
checkpoint_cb =
    HourlyCallback(dt = FT(480), func = checkpoint_sims, ref_date = [dates.date[1]], active = hourly_checkpoint) # 20 days TODO: not GPU friendly
update_firstdayofmonth!_cb =
    MonthlyCallback(dt = FT(1), func = update_firstdayofmonth!, ref_date = [dates.date1[1]], active = true)
callbacks = (; checkpoint = checkpoint_cb, update_firstdayofmonth! = update_firstdayofmonth!_cb)

cs = CoupledSimulation{FT}(
    comms_ctx,
    dates,
    boundary_space,
    coupler_fields,
    config_dict,
    nothing, # conservation checks
    [tspan[1], tspan[2]],
    atmos_sim.integrator.t,
    Δt_cpl,
    (;),
    model_sims,
    (;), # mode_specifics
    (), # coupler diagnostics
    callbacks,
    dir_paths,
    nothing, # turbulent_fluxes
    thermo_params,
);

#=
## Restart component model states if specified in the config_dict
=#

if restart_dir !== "unspecified"
    for sim in cs.model_sims
        if get_model_prog_state(sim) !== nothing
            restart_model_state!(sim, comms_ctx, restart_t; input_dir = restart_dir)
        end
    end
end

#=
## Coupling Loop
=#

function solve_coupler!(cs)
    (; model_sims, Δt_cpl, tspan, comms_ctx) = cs
    (; atmos_sim) = model_sims

    ClimaComms.iamroot(comms_ctx) ? @info("Starting coupling loop") : nothing
    ## step in time
    walltime = @elapsed for t in ((tspan[begin] + Δt_cpl):Δt_cpl:tspan[end])

        cs.dates.date[1] = current_date(cs, t)

        ## print date on the first of month
        if cs.dates.date[1] >= cs.dates.date1[1]
            ClimaComms.iamroot(comms_ctx) ? @show(cs.dates.date[1]) : nothing
        end

        ## step sims
        step_model_sims!(cs.model_sims, t)

        import_atmos_fields!(cs.fields, cs.model_sims, cs.boundary_space, cs.turbulent_fluxes) # radiative and/or turbulent

        ## callback to update the fist day of month if needed (for BCReader)
        trigger_callback!(cs, cs.callbacks.update_firstdayofmonth!)

        ## callback to checkpoint model state
        trigger_callback!(cs, cs.callbacks.checkpoint)

    end
    ClimaComms.iamroot(comms_ctx) ? @show(walltime) : nothing

    return cs
end

## run the coupled simulation
solve_coupler!(cs);
