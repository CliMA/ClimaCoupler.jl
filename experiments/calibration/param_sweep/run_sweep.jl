# Parameter sweep utility for the coupled AMIP model.
#
# Two inputs, given in a small Julia experiment file: a coupler YAML and the
# parameters to vary. The sweep is a single ClimaCalibrate iteration with no
# ensemble update: ClimaCalibrate runs one forward model per member with the
# same `CouplerModelInterface` used by the AMIP calibration. The analysis
# reports, per diagnostic variable, how much the output responds to each
# parameter compared with the model's internal variability (measured by
# near-identical replicate runs). There is no comparison to observations.
#
# The experiment file mixes exact values and priors freely. Per parameter:
#   scalar -> pinned at that value in every member,
#   vector -> swept over exactly those values,
#   prior (PD.constrained_gaussian) -> swept over deterministic quantiles of
#       the prior (default 3, set prior_points to change).
# The ensemble is the full factorial across the swept parameters. See
# example_sweep.jl.
#
# Run (from the repository root):
#
#   julia --project=experiments/AMIP \
#       experiments/calibration/param_sweep/run_sweep.jl my_sweep.jl
#
# ClimaCalibrate submits one GPU worker per member, runs every member on the
# worker pool, then logs the analysis. Keep the driver process alive while the
# members run. Pass --no-submit to write the member inputs without running
# them, or --analyze to only analyze the output of a finished sweep. Both modes
# take the same experiment file, which is the only record of the design.

import Dates
import Random
import Statistics
import TOML
import ClimaAnalysis
import ClimaCalibrate
import ClimaCoupler
import ClimaCoupler: CalibrationTools
import ClimaParams as CP
import EnsembleKalmanProcesses.ParameterDistributions as PD
import EnsembleKalmanProcesses.TOMLInterface as TI

const MODEL_INTERFACE_FILEPATH = joinpath(
    pkgdir(ClimaCoupler),
    "experiments",
    "calibration",
    "amip",
    "model_interface.jl",
)

# The forward model, the parameter file layout, and the analysis paths all come
# from the AMIP calibration interface. `@worker_setup` loads it here and replays
# it on every worker as it joins, unlike `@everywhere`.
ClimaCalibrate.@worker_setup include($MODEL_INTERFACE_FILEPATH)

# A sweep is one iteration that is never updated.
const SWEEP_ITERATION = 0

# One worker per member, each holding one Derecho GPU. Set `worker_options` in
# the experiment file to change these; anything `add_workers` does not consume
# goes to the scheduler, for example `q = "preempt"`.
const DEFAULT_WORKER_OPTIONS = (; device = :gpu, time = 240)

# ---------------------------------------------------------------------------
# Member expansion
# ---------------------------------------------------------------------------

"""
    expand_parameters(parameters; prior_points = 3)

Turn the per-parameter specification into the member list. A scalar pins the
parameter, a vector sweeps exactly those values, and a ParameterDistribution
sweeps `prior_points` deterministic quantiles of the prior. Members are the
full factorial across the swept parameters.
"""
function expand_parameters(parameters; prior_points = 3)
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
    n_members <= 64 || error(
        "The factorial sweep has $n_members members. Shrink the value lists " *
        "or pin more parameters (limit 64).",
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
model ignores is not detectable this cheaply; the analysis flags it
afterwards through members with no response.
"""
function check_parameter_names(names, config_file)
    registry = CP.create_toml_dict(Float32)
    config_dict = ClimaCoupler.Input.get_coupler_config_dict(config_file)
    toml_names = Set{String}()
    for f in config_dict["coupler_toml"]
        path = isfile(f) ? f : joinpath(pkgdir(ClimaCoupler), f)
        isfile(path) && union!(toml_names, keys(TOML.parsefile(path)))
    end
    unknown =
        filter(n -> !(haskey(registry.data, n) || n in toml_names), names)
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
# Analysis: response across members against the replicate noise floor.
# Member 1 is the base point. The response of member m for a variable is the
# latitude-weighted RMS difference of its last monthly mean from member 1's.
# Replicates (0.1% jitter off member 1) measure how large that difference is
# from internal variability alone.
# ---------------------------------------------------------------------------

"""
    weighted_rms_diff(var_a, var_b)

Latitude-weighted RMS difference between two `OutputVar`s.

NaNs are ignored. Dimensions other than longitude and latitude are averaged
with equal weight.
"""
function weighted_rms_diff(var_a, var_b)
    squared_diff = (var_a - var_b) * (var_a - var_b)
    mse = ClimaAnalysis.weighted_average_lonlat(squared_diff)
    finite_mse = filter(isfinite, mse.data)
    return isempty(finite_mse) ? NaN : sqrt(Statistics.mean(finite_mse))
end

"Last monthly mean of `short_name`, or `nothing` when unavailable."
function load_last_month(simdir_path, short_name)
    isdir(simdir_path) || return nothing
    var = try
        get(
            ClimaAnalysis.SimDir(simdir_path);
            short_name,
            reduction = "average",
            period = "1M",
        )
    catch
        return nothing
    end
    return ClimaAnalysis.slice(var; time = last(ClimaAnalysis.times(var)))
end

"Path to a member's diagnostics folder."
function member_simdir_path(config, member)
    member_path = ClimaCalibrate.path_to_ensemble_member(
        config.output_dir,
        SWEEP_ITERATION,
        member,
    )
    return joinpath(member_path, get_job_id(config), "output_active")
end

"""
    analyze(sweep)

Report the response of every member of `sweep` against the replicate noise
floor, and flag parameters that the loaded model configuration ignores.
"""
function analyze(sweep)
    (; interface, n_sweep, n_replicates) = sweep
    (; config) = interface
    (; short_names, output_dir) = config
    n_members = n_sweep + n_replicates

    members = map(1:n_members) do m
        toml = TOML.parsefile(
            ClimaCalibrate.parameter_path(output_dir, SWEEP_ITERATION, m),
        )
        Dict{String, Float64}(k => v["value"] for (k, v) in toml)
    end
    param_names = sort!(collect(keys(first(members))))

    available = filter(m -> isdir(member_simdir_path(config, m)), 1:n_members)
    missing_members = setdiff(1:n_members, available)
    isempty(missing_members) ||
        @warn "Members without output are skipped" missing_members
    1 in available || error("Member 1 (the base member) has no output")

    fields = Dict(
        m => Dict(
            n => load_last_month(member_simdir_path(config, m), n) for
            n in short_names
        ) for m in available
    )

    responses = Dict{String, Dict{Int, Float64}}()
    for name in short_names
        ref = fields[1][name]
        if isnothing(ref)
            @warn "Variable $name missing in member 1, skipped"
            continue
        end
        response = Dict(
            m => weighted_rms_diff(fields[m][name], ref) for
            m in available if !isnothing(fields[m][name])
        )
        responses[name] = response

        sweep_r = [r for (m, r) in response if m <= n_sweep]
        rep_r = [r for (m, r) in response if m > n_sweep]
        signal = isempty(sweep_r) ? NaN : maximum(sweep_r)
        noise = isempty(rep_r) ? NaN : Statistics.mean(rep_r)
        @info "  $name" max_response = signal replicate_noise = noise signal_to_noise =
            signal / noise
    end

    # Member 1 is the base corner of the factorial, so every swept parameter has
    # a member that differs from member 1 in that parameter alone. A member with
    # no response at all reproduced the base output bit for bit.
    for m in available
        (m != 1 && m <= n_sweep && !isempty(responses)) || continue
        all(r -> get(r, m, NaN) == 0, values(responses)) || continue
        differing = [n for n in param_names if members[m][n] != members[1][n]]
        @error "Member $m produced bit-identical output to member 1 although " *
               "they differ in $(differing). Those parameters are very " *
               "likely NOT wired into the loaded model configuration."
    end
    return nothing
end

# ---------------------------------------------------------------------------
# Design and run: expand the experiment file, write the ClimaCalibrate inputs,
# run the ensemble, analyze.
# ---------------------------------------------------------------------------

"""
    design(experiment_path)

Expand the experiment file into the members, the model interface, and the worker
settings of a single iteration. This writes nothing, so the analysis can rebuild
the design of a finished sweep from the same experiment file.
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

    sweep_members = expand_parameters(
        Main.parameters;
        prior_points = _get(:prior_points, 3),
    )
    param_names = sort!(collect(keys(first(sweep_members))))
    check_parameter_names(param_names, config_file)

    # Replicates jitter the first member by 0.1% so their output spread
    # measures the model's internal variability, the floor any real
    # parameter response must beat.
    n_replicates = _get(:replicates, 3)
    rng = Random.MersenneTwister(_get(:rng_seed, 42))
    base = first(sweep_members)
    replicate_members = [
        Dict(n => v * (1 + 1e-3 * randn(rng)) for (n, v) in base) for
        _ in 1:n_replicates
    ]
    members = vcat(sweep_members, replicate_members)
    @info "Sweep design" n_sweep = length(sweep_members) n_replicates param_names

    # One week of spin-up, then the sampled month.
    sample_date = Dates.DateTime(_get(:sample_date, "2010-10-01"))
    config = CalibrationTools.CalibrateConfig(;
        config_file,
        short_names = collect(
            String,
            _get(:short_names, ["lwp", "pr", "rsut", "rlut", "tas"]),
        ),
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
        n_sweep = length(sweep_members),
        n_replicates,
        n_workers = _get(:n_workers, length(members)),
        worker_options = _get(:worker_options, DEFAULT_WORKER_OPTIONS),
    )
end

"""
    run_sweep(experiment_path; submit)

Run every member of the sweep described by `experiment_path`, then analyze the
output. Members that already completed are not rerun.
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

    # Workers are submitted as individual allocations and join the pool as they
    # start, so the members begin before the whole pool is up. Fewer workers
    # than members just means the members run in waves.
    @info "Submitting workers" n_workers n_members = length(members) worker_options
    ClimaCalibrate.add_workers(n_workers; worker_options...)

    ClimaCalibrate.Calibration.run_iteration(
        ClimaCalibrate.WorkerBackend(),
        interface,
        SWEEP_ITERATION,
        length(members),
        output_dir,
    )
    analyze(sweep)
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    isempty(ARGS) &&
        error("Usage: run_sweep.jl <sweep.jl> [--no-submit | --analyze]")
    experiment_path = abspath(ARGS[1])
    analyze_only = "--analyze" in ARGS
    submit = !("--no-submit" in ARGS)
    empty!(ARGS)
    if analyze_only
        analyze(design(experiment_path))
    else
        run_sweep(experiment_path; submit)
    end
end
