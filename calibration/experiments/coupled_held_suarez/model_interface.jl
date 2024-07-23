import ClimaCalibrate: set_up_forward_model, run_forward_model, path_to_ensemble_member, ExperimentConfig
import ClimaComms
@static pkgversion(ClimaComms) >= v"0.6" && ClimaComms.@import_required_backends
using ClimaUtilities.ClimaArtifacts
import YAML

using ClimaCoupler

"""
    set_up_forward_model(member, iteration, experiment_dir::AbstractString)

Return CoupledSimulation object for the given member and iteration. 

Turns off default diagnostics and sets the TOML parameter file to the member's path.
"""
function set_up_forward_model(member, iteration, experiment_dir::AbstractString)
    include(joinpath(experiment_dir,"forward_model.jl"));
    output_dir = cs.parsed_args["coupler_output_dir"]
    cs.parsed_args["coupler_output_dir"] = joinpath(experiment_dir,"output", output_dir);
    return cs
end

"""
    run_forward_model(coupled_simulation)

Run the coupled model with the given a coupled simulation object.
Currently only has basic error handling.
"""
function run_forward_model(cs)
    sol_res = solve_coupler!(cs);
    return sol_res
end
