import EnsembleKalmanProcesses as EKP
import ClimaCoupler as CCo
import YAML
import CalibrateAtmos: get_forward_model, AbstractPhysicalModel, get_config

struct CoupledModel <: AbstractPhysicalModel end

function get_forward_model(
    experiment_id::Val{:amip_coupled}
)
    return CoupledModel()
end

function get_config(
    model::CoupledModel,
    member,
    iteration, 
    experiment_id::AbstractString
)
    config_dict = YAML.load_file("experiments/$experiment_id/model_config.yml")
    return get_config(model, member, iteration, config_dict)
end

function get_config(
    ::CoupledModel,
    member,
    iteration,
    config_dict::AbstractDict,
)
    # Specify member path for output_dir
    # Set TOML to use EKP parameter(s)
    config_dict = YAML.load_file("./experiments/amip_coupled/model_config.yml")
    output_dir = "output"
    member_path =
        EKP.TOMLInterface.path_to_ensemble_member(output_dir, iteration, member)
    config_dict = merge(parsed_args, config_dict)

    ## get component model dictionaries (if applicable)
    config_dict_atmos = get_atmos_config(config_dict)
   
    ## merge dictionaries of command line arguments, coupler dictionary and component model dictionaries
    ## (if there are common keys, the last dictorionary in the `merge` arguments takes precedence)
    config_dict = merge(config_dict_atmos, config_dict)
    # COPY Coupler Driver
    config_dict["output_dir"] = member_path
    include("coupler_driver_calibration.jl")
    config_dict["output_dir"] = member_path
    # END Coupler Driver
    parameter_path = joinpath(member_path, "parameters.toml")
    if haskey(config_dict, "toml")
        push!(config_dict["toml"], parameter_path)
    else
        config_dict["toml"] = [parameter_path]
    end
    # Turn off default diagnostics
    config_dict["output_default_diagnostics"] = false
    return (;config_dict=config_dict)
end

function run_forward_model(
    ::CoupledModel,
    config;
    lk = nothing,
)   
    cs = get_simulation(config);
    sol_res = solve_coupler!(cs);
    return sol_res
end
