# Parameter sweep utility for the coupled AMIP model.
#
# Two inputs, given in a small Julia experiment file: a coupler YAML and the
# parameters to vary. The sweep is a single ClimaCalibrate iteration with no
# ensemble update: ClimaCalibrate runs one forward model per member with the
# same `CouplerModelInterface` used by the AMIP calibration. Each member writes
# its own model output; comparing the members is left to the user.
#
# The experiment file mixes exact values and priors freely. Per parameter:
#   scalar -> pinned at that value in every member,
#   vector -> swept over exactly those values,
#   prior (PD.constrained_gaussian) -> swept over evenly spaced quantiles of
#       the prior (default 3, set prior_points to change).
# The members are every combination of the swept values. See example_sweep.jl.
#
# Run (from the repository root):
#
#   julia --project=experiments/AMIP \
#       experiments/calibration/param_sweep/run_sweep.jl my_sweep.jl
#
# ClimaCalibrate submits GPU workers and runs every member on the worker pool.
# Keep the driver process alive while the members run. Pass --no-submit to write
# the member inputs without running them. The experiment file is the only record
# of the design, so keep it next to the sweep.

import Dates
import TOML
import ClimaCalibrate
import ClimaCoupler
import ClimaCoupler: CalibrationTools
import ClimaParams as CP
import EnsembleKalmanProcesses.ParameterDistributions as PD
import EnsembleKalmanProcesses.TOMLInterface as TI
using Distributed

const MODEL_INTERFACE_FILEPATH = joinpath(
    pkgdir(ClimaCoupler),
    "experiments",
    "calibration",
    "amip",
    "model_interface.jl",
)

# The forward model and the parameter file layout both come from the AMIP
# calibration interface. The driver needs it to build the interface; the workers
# get it again once they have joined, in `run_sweep`.
include(MODEL_INTERFACE_FILEPATH)

# A sweep is one iteration that is never updated.
const SWEEP_ITERATION = 0

# One worker per member, each holding one Derecho GPU. Set `worker_options` in
# the experiment file to change these; anything `add_workers` does not consume
# goes to the scheduler.
const DEFAULT_WORKER_OPTIONS = (; device = :gpu, time = 240)

# Guard against a value list that accidentally asks for hundreds of 1.2-hour GPU
# runs. Raise it with `max_members` in the experiment file if a large sweep is
# what you want.
const DEFAULT_MAX_MEMBERS = 64

# ---------------------------------------------------------------------------
# Member expansion
# ---------------------------------------------------------------------------

"""
    expand_parameters(parameters; prior_points = 3, max_members = 64)

Turn the per-parameter specification into the member list.

A scalar pins the parameter, a vector sweeps exactly those values, and a
ParameterDistribution sweeps `prior_points` evenly spaced quantiles of the
prior: quantile `i / (prior_points + 1)` for `i in 1:prior_points`, so the
default of 3 points gives the 25th, 50th, and 75th percentiles. The endpoints
are deliberately left out, since a prior's 0th and 100th percentiles are often
unphysical.

Members are every combination of the swept values, so the count is the product
of the value-list lengths. `max_members` caps that count.
"""
function expand_parameters(parameters; prior_points = 3, max_members = DEFAULT_MAX_MEMBERS)
    names = sort!(collect(String, keys(parameters)))
    value_axes = map(names) do n
        spec = parameters[n]
        if spec isa Number
            [Float64(spec)]
        elseif spec isa AbstractVector
            Float64.(collect(spec))
        elseif spec isa PD.ParameterDistribution
            unconstrained = only(values(PD.get_distribution(spec)))
            map((i / (prior_points + 1) for i in 1:prior_points)) do q
                z = PD.quantile(unconstrained, q)
                only(PD.transform_unconstrained_to_constrained(spec, [z]))
            end
        else
            error(
                "parameter $n must be a number, a vector of values, or a " *
                "ParameterDistribution, got $(typeof(spec))",
            )
        end
    end
    n_members = prod(length.(value_axes))
    n_members <= max_members || error(
        "This sweep has $n_members members, which is more than " *
        "max_members = $max_members. Every member is a separate ~1.2 hour " *
        "GPU run, so shrink the value lists, pin more parameters, or raise " *
        "max_members in the experiment file.",
    )
    members = [
        Dict{String, Float64}(zip(names, Float64.(collect(combo)))) for
        combo in Iterators.product(value_axes...)
    ]
    return vec(members)
end

"""
    check_parameter_names(names, config_file)

Every parameter name must exist in the ClimaParams registry or the coupler
TOML files. This catches typos in seconds. A registered parameter the loaded
model ignores is not detectable this cheaply.
"""
function check_parameter_names(names, config_file)
    registry = CP.create_toml_dict(Float32)
    config_dict = ClimaCoupler.Input.get_coupler_config_dict(config_file)
    toml_names = Set{String}()
    for f in config_dict["coupler_toml"]
        path = isfile(f) ? f : joinpath(pkgdir(ClimaCoupler), f)
        isfile(path) && union!(toml_names, keys(TOML.parsefile(path)))
    end
    unknown = filter(n -> !(haskey(registry.data, n) || n in toml_names), names)
    isempty(unknown) || error(
        "Unknown parameter names: $(unknown). Not in the ClimaParams " *
        "registry or the coupler TOML files, so they are almost surely typos.",
    )
    return nothing
end

"""
    write_parameter_ensemble(members, output_dir)

Write one parameter file per member into the ClimaCalibrate iteration layout.

`save_parameter_ensemble` writes the member parameter files from a distribution
and a value array. Sweep values are already in the constrained space, so the
distribution only has to carry the parameter names.
"""
function write_parameter_ensemble(members, output_dir)
    names = sort!(collect(String, keys(first(members))))
    prior = mapreduce(
        name -> PD.ParameterDistribution(
            Dict(
                "distribution" => PD.Parameterized(PD.Normal(0, 1)),
                "constraint" => PD.no_constraint(),
                "name" => name,
            ),
        ),
        (x, y) -> PD.combine_distributions([x, y]),
        names,
    )
    values = [member[name] for name in names, member in members]
    TI.save_parameter_ensemble(
        values,
        prior,
        ClimaCalibrate.get_param_dict(prior),
        output_dir,
        "parameters.toml",
        SWEEP_ITERATION;
        apply_constraints = false,
    )
    return nothing
end

# ---------------------------------------------------------------------------
# Design and run: expand the experiment file, write the ClimaCalibrate inputs,
# run the ensemble.
# ---------------------------------------------------------------------------

"""
    design(experiment_path)

Expand the experiment file into the members, the model interface, and the worker
settings of a single iteration. This writes nothing.
"""
function design(experiment_path)
    Base.include(Main, experiment_path)
    for key in (:config_file, :output_dir, :parameters)
        isdefined(Main, key) || error("experiment file must define $key")
    end
    _get(key, default) = isdefined(Main, key) ? getfield(Main, key) : default

    config_file = abspath(Main.config_file)
    isfile(config_file) || error("config_file $config_file does not exist")
    output_dir = Main.output_dir
    isabspath(output_dir) ||
        error("output_dir must be an absolute path on scratch, got $output_dir")
    occursin("CHANGE_ME", output_dir) && error(
        "output_dir is still the placeholder from example_sweep.jl. Point it " *
        "at your own scratch directory: $output_dir",
    )

    members = expand_parameters(
        Main.parameters;
        prior_points = _get(:prior_points, 3),
        max_members = _get(:max_members, DEFAULT_MAX_MEMBERS),
    )
    param_names = sort!(collect(keys(first(members))))
    check_parameter_names(param_names, config_file)
    @info "Sweep design" n_members = length(members) param_names

    # One week of spin-up, then the sampled month.
    sample_date = Dates.DateTime(_get(:sample_date, "2010-10-01"))
    config = CalibrationTools.CalibrateConfig(;
        config_file,
        # `CalibrateConfig` requires `short_names` because a calibration needs
        # observations to compare against, but these are never used in the sweep.
        short_names = ["ta"],
        minibatch_size = 1,
        n_iterations = 1,
        sample_date_ranges = [(sample_date, sample_date)],
        extend = Dates.Month(1),
        spinup = Dates.Day(_get(:spinup_days, 7)),
        output_dir,
        rng_seed = _get(:rng_seed, 42),
    )
    return (;
        interface = CouplerModelInterface(config),
        members,
        n_workers = _get(:n_workers, length(members)),
        worker_options = _get(:worker_options, DEFAULT_WORKER_OPTIONS),
    )
end

"""
    run_sweep(experiment_path; submit)

Run every member of the sweep described by `experiment_path`. Members that
already completed are not rerun, and a member that stopped partway resumes from
its last checkpoint.
"""
function run_sweep(experiment_path; submit)
    sweep = design(experiment_path)
    (; interface, members, n_workers, worker_options) = sweep
    output_dir = interface.config.output_dir
    write_parameter_ensemble(members, output_dir)
    if !submit
        @info "Wrote the member inputs but did not run them" output_dir
        return nothing
    end

    # Each worker is its own scheduler allocation holding one GPU. Members are
    # handed to workers as they free up, so fewer workers than members just
    # means the members run in waves.
    @info "Submitting workers" n_workers n_members = length(members) worker_options
    ClimaCalibrate.add_workers(n_workers; worker_options...)
    @everywhere include($MODEL_INTERFACE_FILEPATH)

    ClimaCalibrate.run_iteration(
        ClimaCalibrate.WorkerBackend(),
        interface,
        SWEEP_ITERATION,
        length(members),
        output_dir,
    )
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    isempty(ARGS) && error("Usage: run_sweep.jl <sweep.jl> [--no-submit]")
    experiment_path = abspath(ARGS[1])
    submit = !("--no-submit" in ARGS)
    empty!(ARGS)
    run_sweep(experiment_path; submit)
end
