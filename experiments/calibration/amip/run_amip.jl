# TOML front end for the three-input AMIP pipeline.
#
# Usage:
#   julia --project=experiments/AMIP experiments/calibration/amip/run_amip.jl \
#       my_experiment.toml [--no-submit]
#
# See example_experiment.toml for the file format and README.md for the
# workflow. The TOML carries the three inputs (coupler YAML, parameters,
# observations); everything else has validated defaults.

import Dates
import TOML

include(joinpath(@__DIR__, "pipeline.jl"))

function _parse_dates(entries)
    return [
        begin
            d = Dates.DateTime(String(e))
            (d, d)
        end for e in entries
    ]
end

function _parse_experiment(path)
    spec = TOML.parsefile(path)
    haskey(spec, "config_file") || error("experiment TOML needs config_file")
    haskey(spec, "output_dir") || error("experiment TOML needs output_dir")
    haskey(spec, "observations") || error("experiment TOML needs observations")

    haskey(spec, "priors") || error("experiment TOML needs [[priors]] tables")
    parameters = [
        PD.constrained_gaussian(
            p["name"],
            Float64(p["mean"]),
            Float64(p["sigma"]),
            Float64(p["lower"]),
            Float64(p["upper"]),
        ) for p in spec["priors"]
    ]

    kwargs = Dict{Symbol, Any}(
        :parameters => parameters,
        :observations => spec["observations"],
        :output_dir => spec["output_dir"],
    )
    if haskey(spec, "reduction")
        r = spec["reduction"]
        kwargs[:reduction] = r == "zonal" ? :zonal : Int(r)
    end
    haskey(spec, "short_names") && (kwargs[:short_names] = spec["short_names"])
    haskey(spec, "n_iterations") && (kwargs[:n_iterations] = spec["n_iterations"])
    haskey(spec, "rng_seed") && (kwargs[:rng_seed] = spec["rng_seed"])
    haskey(spec, "sample_dates") &&
        (kwargs[:sample_date_ranges] = _parse_dates(spec["sample_dates"]))
    haskey(spec, "covariance_dates") &&
        (kwargs[:covariance_date_ranges] = _parse_dates(spec["covariance_dates"]))
    haskey(spec, "spinup_days") && (kwargs[:spinup] = Dates.Day(spec["spinup_days"]))
    haskey(spec, "extend_months") && (kwargs[:extend] = Dates.Month(spec["extend_months"]))
    haskey(spec, "workers_per_node") &&
        (kwargs[:workers_per_node] = spec["workers_per_node"])
    if haskey(spec, "noise_groups")
        kwargs[:noise_groups] = [
            (
                short_names = collect(String, g["short_names"]),
                model_error_scale = Float64(g["model_error_scale"]),
            ) for g in spec["noise_groups"]
        ]
    end
    return spec["config_file"], kwargs
end

if abspath(PROGRAM_FILE) == @__FILE__
    isempty(ARGS) && error("Usage: run_amip.jl <experiment.toml> [--no-submit]")
    experiment_path = abspath(ARGS[1])
    submit = !("--no-submit" in ARGS)
    empty!(ARGS)

    config_file, kwargs = _parse_experiment(experiment_path)
    result = amip_pipeline(config_file; submit, kwargs...)
    @info "Pipeline finished" result.config_path result.job_id
end
