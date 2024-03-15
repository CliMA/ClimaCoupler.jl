import EnsembleKalmanProcesses as EKP
import ClimaCoupler as CCo
import YAML

"""
    get_coupler_sim(member, iteration, experiment_id::AbstractString)

Returns a CouplerSimulation object for the given member and iteration. 
If given an experiment id string, it will load the config from the corresponding YAML file.
Turns off default diagnostics and sets the TOML parameter file to the member's path.
This assumes that the config dictionary has `output_dir` and `restart_file` keys.
"""
function get_coupler_sim(member, iteration, experiment_id::AbstractString)
    # Specify member path for output_dir
    # Set TOML to use EKP parameter(s)
    config_dict = YAML.load_file("./experiments/amip_coupled/coupler_config.yml")
    output_dir = "output"
    member_path =
        EKP.TOMLInterface.path_to_ensemble_member(output_dir, iteration, member)
    config_dict["output_dir"] = member_path
    # COPY Coupler Driver
    include("coupler_driver_calibration.jl")
    # END Coupler Driver

    parameter_path = joinpath(member_path, "parameters.toml")
    if haskey(config_dict, "toml")
        push!(config_dict["toml"], parameter_path)
    else
        config_dict["toml"] = [parameter_path]
    end

    # Turn off default diagnostics
    config_dict["output_default_diagnostics"] = false

    # Set restart file for initial equilibrium state
    return cs
end

"""
    run_forward_model(coupled_sim::CCo.CoupledSimulation)

Runs the coupled model with the given a CoupledSimulation object.
Note that running an AtmosModel can be considered a special case 
of running a CoupledSimulation. 
Currently only has basic error handling.
"""
function run_forward_model(coupled_sim::CCo.Interfacer.CoupledSimulation)
    sol_res = solve_coupler!(coupled_sim)
    # Is there a return code in coupler sims? 
    # Can we get the atmos version instead ?
end
