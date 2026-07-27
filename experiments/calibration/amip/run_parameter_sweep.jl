# Parameter sweep / sensitivity harness for the AMIP calibration.
#
# This runs an ARBITRARY ensemble (any set of parameter values, any size) through
# the exact same ClimaCalibrate machinery the calibration uses -- it writes one
# parameter TOML per member, submits the forward models via the backend, and runs
# the observation map to get the forward-model output matrix G (N_obs x N_ens).
#
# Its purpose is to answer the question the flat-error calibration could not:
# does the parameter actually move the observables ABOVE the model's internal
# (chaotic) variability on a single monthly sample? We do this by including, in a
# single ensemble:
#   * a WIDE sweep of the parameter across its prior range (the "signal"), and
#   * a tight cluster of near-default replicates (the "noise floor" -- tiny
#     parameter jitter whose output spread reflects chaotic divergence over the
#     1-month window, i.e. the variability that competes with the signal).
# If the sweep response is not clearly larger than the replicate scatter, no
# amount of retuning noise/priors will make this parameter identifiable on this
# observation with a single-month sample.
#
# Usage (same as run_calibration.jl):
#   julia --project=experiments/AMIP/ experiments/calibration/amip/run_parameter_sweep.jl
#
# Edit the SWEEP SPECIFICATION block below to sweep a different parameter / grid.

using Dates
import Random
import Statistics
import ClimaCalibrate
import ClimaAnalysis
import ClimaCoupler
import ClimaCoupler: CalibrationTools
import EnsembleKalmanProcesses as EKP
import EnsembleKalmanProcesses.ParameterDistributions as PD
import JLD2
import CairoMakie

# Brings in CouplerModelInterface, forward_model, observation_map, and (via the
# config include below) the priors + calibration config.
include(joinpath(@__DIR__, "model_interface.jl"))
# Config-selectable like run_calibration.jl. Defaults to the lwp+pr config so the
# sweep probes the warm-rain parameters against {lwp, pr}. Override with
# ENV["CALIBRATION_CONFIG"].
_sweep_config_path = get(
    ENV,
    "CALIBRATION_CONFIG",
    joinpath(@__DIR__, "config", "lwp_pr.jl"),
)
@info "Sweep using calibration configuration in: $_sweep_config_path"
include(_sweep_config_path)

# ---------------------------------------------------------------------------
# Analysis: per-variable RMSE vs parameter, with the replicate noise floor.
# ---------------------------------------------------------------------------
"""
    analyze_sweep(obs_series, g_ens, samples, param_names, is_replicate, iteration, output_dir)

Report, per observed variable, how much the forward-model output responds to the
swept parameters (the "signal", = spread of RMSE-to-obs across the wide sweep)
relative to the near-default replicate scatter (the "noise floor", = std of
RMSE-to-obs across replicates), and save a figure of RMSE-to-obs vs each swept
parameter (one row per variable, one column per parameter).

`samples` is the `N_param x N_ens` matrix of physical parameter values, ordered
by `param_names`.

A signal-to-noise ratio well above 1 means the parameters are identifiable on
this observation; a ratio near or below 1 means their effect is buried in the
model's internal variability and calibration cannot succeed as configured. The
per-parameter columns show WHICH parameter drives the response for each variable.
"""
function analyze_sweep(
    obs_series,
    g_ens,
    samples,
    param_names,
    is_replicate,
    iteration,
    output_dir,
)
    # Observation target and per-variable index ranges for this iteration's
    # minibatch, reconstructed from the EKP observation metadata.
    mb_obs = ClimaCalibrate.get_observations_for_nth_iteration(obs_series, iteration)
    y = reduce(vcat, EKP.get_obs.(mb_obs))
    obs_vars =
        mapreduce(ClimaCalibrate.ObservationRecipe.reconstruct_vars, vcat, mb_obs)
    var_names = ClimaAnalysis.short_name.(obs_vars)
    var_lengths = [length(ClimaAnalysis.flatten(v).data) for v in obs_vars]
    @assert sum(var_lengths) == length(y) == size(g_ens, 1) "obs/G length mismatch"

    # NaN-aware RMSE between a model column and the observation over a slice.
    function slice_rmse(gcol, lo, hi)
        r = gcol[lo:hi] .- y[lo:hi]
        finite = filter(isfinite, r)
        isempty(finite) ? NaN : sqrt(Statistics.mean(abs2, finite))
    end

    sweep_idx = findall(!, is_replicate)
    rep_idx = findall(is_replicate)
    n_param = length(param_names)

    @info "=== Parameter sweep sensitivity ==="
    @info "Per-variable RMSE-to-obs. signal = RMSE spread across the wide sweep; " *
          "noise floor = RMSE std across near-default replicates. " *
          "signal/noise >> 1 means the parameters are identifiable."

    fig = CairoMakie.Figure(size = (500 * n_param, 420 * length(obs_vars)))
    offset = 0
    for (vi, (name, len)) in enumerate(zip(var_names, var_lengths))
        lo, hi = offset + 1, offset + len
        offset += len
        rmse = [slice_rmse(view(g_ens, :, m), lo, hi) for m in 1:size(g_ens, 2)]

        sweep_rmse = rmse[sweep_idx]
        rep_rmse = rmse[rep_idx]
        signal = maximum(sweep_rmse) - minimum(sweep_rmse)
        noise = length(rep_rmse) > 1 ? Statistics.std(rep_rmse) : NaN
        ratio = signal / noise
        @info "  $name" signal_rmse_range = signal noise_floor_std = noise signal_to_noise =
            ratio

        # Reachability: does the observed mean fall inside the model's attainable
        # range across the sweep? If not, the bias is structural for this
        # parameter — no value of it can reach the observation.
        blockmean(v) = Statistics.mean(filter(isfinite, v))
        obs_mean = blockmean(view(y, lo:hi))
        model_means = [blockmean(view(g_ens, lo:hi, m)) for m in sweep_idx]
        m_lo, m_hi = extrema(model_means)
        @info "  $name reachability" obs_mean model_range = (m_lo, m_hi) reachable =
            (m_lo <= obs_mean <= m_hi)

        rep_mean = length(rep_rmse) > 1 ? Statistics.mean(rep_rmse) : NaN
        rep_std = length(rep_rmse) > 1 ? Statistics.std(rep_rmse) : NaN
        for p in 1:n_param
            ax = CairoMakie.Axis(
                fig[vi, p];
                title = p == 1 ? "$name  (overall S/N ≈ $(round(ratio; digits = 2)))" :
                        name,
                xlabel = String(param_names[p]),
                ylabel = "RMSE to obs",
            )
            CairoMakie.scatter!(
                ax,
                samples[p, sweep_idx],
                sweep_rmse;
                label = "sweep",
            )
            # Replicate noise floor: mean +/- std band.
            if length(rep_rmse) > 1
                CairoMakie.hlines!(ax, [rep_mean]; color = :red, label = "replicate mean")
                CairoMakie.hspan!(ax, rep_mean - rep_std, rep_mean + rep_std; color = (:red, 0.15))
            end
            CairoMakie.axislegend(ax; position = :rt)
        end
    end
    figpath = joinpath(output_dir, "parameter_sweep_sensitivity.png")
    CairoMakie.save(figpath, fig)
    @info "Saved sensitivity figure" figpath
    return nothing
end

# ---------------------------------------------------------------------------
# SWEEP SPECIFICATION -- edit this block.
#
# `SWEEP_PRIORS` lists ONLY the parameters written to each member's TOML (every
# other model parameter keeps its model default). Then, in matching order:
#   * `SWEEP_GRIDS`  : per-parameter vectors of PHYSICAL values to sweep. The
#     ensemble is the FULL-FACTORIAL product of these, so the sweep size is
#     prod(length.(SWEEP_GRIDS)). Cost scales as k^P -- keep grids small.
#   * `SWEEP_DEFAULTS`: the default value of each parameter, used for the
#     near-default replicate cluster that measures the internal-variability
#     (chaotic) noise floor.
#
# The default below is a 2-D grid over a warm-rain parameter (moves lwp) and a
# mixing-length parameter (moves the ta/hur profiles), so a single run tells you
# which knob -- if any -- each observed variable actually responds to.
# ---------------------------------------------------------------------------

# IDENTIFIABILITY sweep: the two warm-rain autoconversion parameters, swept
# jointly on a 2-D grid, to answer whether precipitation (pr) breaks the
# q_liq/rain_tau degeneracy that LWP alone leaves. Both raise LWP by suppressing
# autoconversion (the reservoir); the hypothesis is that pr (the sink) responds
# to rain_tau in a way LWP does not, so the two observables together separate the
# parameters. analyze_sweep reports, per observed variable, the RMSE-to-obs
# response to EACH parameter (one column per param) plus reachability — read the
# lwp vs pr rows side by side to judge independence. σ just parameterizes the
# transform; explicit values round-trip exactly as long as they lie in (lower,
# upper).
const SWEEP_PRIORS = EKP.combine_distributions([
    PD.constrained_gaussian(
        "cloud_liquid_water_specific_humidity_autoconversion_threshold",
        5e-4, 4e-4, 1e-5, 3e-3,
    ),
    PD.constrained_gaussian("rain_autoconversion_timescale", 1500, 1000, 200, 4000),
])

# Grids (full-factorial → 3×3 = 9 sweep members). q_liq brackets the ~7e-4 the
# LWP-only runs converged to; rain_tau spans its physical prior range. Cost is
# k^P, so keep both at 3 points.
const SWEEP_GRIDS = [[2e-4, 7e-4, 1.2e-3], [500.0, 1500.0, 3000.0]]
const SWEEP_DEFAULTS = [5e-4, 1000.0]

# Full-factorial product of the per-parameter grids -> N_param x N_sweep matrix.
sweep_matrix =
    reduce(hcat, [collect(pt) for pt in vec(collect(Iterators.product(SWEEP_GRIDS...)))])

# Near-default replicate cluster (small ±1% jitter on both params) that anchors
# the baseline and gives a local-gradient noise band. NOTE: the model is
# deterministic over a single window, so this band is a proxy for "how flat is
# the response near default", not true internal variability — the real
# interannual noise floor is already encoded in the observation's SVDplusD
# covariance. The primary read is the RMSE-vs-parameter response curves and
# reachability, not this band.
jitter = [-0.01, 0.0, 0.01]
replicate_matrix = reduce(hcat, [SWEEP_DEFAULTS .* (1.0 .+ j) for j in jitter])

# Columns = members: sweep members first, then replicate members.
constrained_samples = hcat(sweep_matrix, replicate_matrix)
is_replicate =
    vcat(falses(size(sweep_matrix, 2)), trues(size(replicate_matrix, 2)))

const ITERATION = 1  # single-iteration "sweep"; mirrors calibration iteration 1
# Write to scratch, NOT the home checkout: the coupled model dumps GBs of
# output + JLD2/HDF5 checkpoints per member, which blows the 100 GiB home quota
# and makes checkpoint writes fail with EOFError. Scratch is large and is where
# the calibrations run.
const SWEEP_OUTPUT_DIR = "/glade/derecho/scratch/nefrathe/amip_parameter_sweep_lwp_pr"

# Rebuild the calibration config pointed at the sweep output directory so this
# never clobbers a real calibration.
(;
    config_file,
    short_names,
    minibatch_size,
    n_iterations,
    sample_date_ranges,
    extend,
    spinup,
    rng_seed,
) = CALIBRATE_CONFIG

sweep_config = CalibrationTools.CalibrateConfig(;
    config_file,
    short_names,
    minibatch_size,
    n_iterations,
    sample_date_ranges,
    extend,
    spinup,
    output_dir = SWEEP_OUTPUT_DIR,
    rng_seed,
)

if abspath(PROGRAM_FILE) == @__FILE__
    isdir(SWEEP_OUTPUT_DIR) || mkpath(SWEEP_OUTPUT_DIR)

    # generate_observations.jl writes the observation vector into the config's
    # output_dir (each setup keeps its own), so read it from there.
    observation_vector_filepath = joinpath(output_dir, "observation_vec.jld2")
    isfile(observation_vector_filepath) || error(
        "Observation vector not found at $observation_vector_filepath. " *
        "Generate it with generate_observations.jl (matching CALIBRATION_CONFIG) first.",
    )
    observation_vector = JLD2.load_object(observation_vector_filepath)

    obs_series = EKP.ObservationSeries(
        Dict(
            "observations" => observation_vector,
            "names" => [
                string(Dates.year(start_date)) for
                (start_date, _) in sample_date_ranges
            ],
            "minibatcher" => ClimaCalibrate.minibatcher_over_samples(
                length(observation_vector),
                minibatch_size,
            ),
        ),
    )

    n_ens = size(constrained_samples, 2)
    @info "Parameter sweep" n_ens n_sweep = count(!, is_replicate) n_replicate =
        count(is_replicate) output_dir = SWEEP_OUTPUT_DIR

    # Build the ensemble directly from our chosen values. EKP's Inversion process
    # accepts an arbitrary-size initial ensemble; we hand it exactly our grid
    # (transformed to the unconstrained space EKP stores internally). We never
    # call update_ensemble!, so these values are what actually get run.
    u = EKP.transform_constrained_to_unconstrained(SWEEP_PRIORS, constrained_samples)
    rng = Random.MersenneTwister(rng_seed)
    ekp = EKP.EnsembleKalmanProcess(u, obs_series, EKP.Inversion(); rng, verbose = true)

    interface = CouplerModelInterface(sweep_config)

    # Writes ekp + one parameter TOML per member under
    # SWEEP_OUTPUT_DIR/iteration_001/member_xxx/.
    ClimaCalibrate.initialize(ekp, SWEEP_PRIORS, SWEEP_OUTPUT_DIR)
    # HPCBackend job scripts load the interface from disk.
    JLD2.save_object(joinpath(SWEEP_OUTPUT_DIR, "interface.jld2"), interface)

    # --- Run the ensemble on the backend ------------------------------------
    env_vars = ["CLIMACOMMS_CONTEXT" => "SINGLETON", "CLIMACOMMS_DEVICE" => "CUDA"]
    backend = if ClimaCalibrate.get_backend() == ClimaCalibrate.DerechoBackend
        ClimaCalibrate.DerechoBackend(;
            directives = [
                :job_priority => "regular",
                :time => 720,
                :ntasks => 1,
                :cpus_per_task => 12,
                :gpus_per_task => 1,
            ],
            modules = ["climacommon/2025_02_25"],
            env_vars,
        )
    else
        error("Unsupported backend: $(ClimaCalibrate.get_backend())")
    end

    # run_iteration lives in the Calibration submodule (only `calibrate` is
    # re-exported to the ClimaCalibrate top level), so qualify it explicitly.
    ClimaCalibrate.Calibration.run_iteration(
        backend,
        ITERATION,
        n_ens,
        SWEEP_OUTPUT_DIR,
        abspath(ClimaCalibrate.model_interface_filepath(interface)),
        abspath(ClimaCalibrate.experiment_dir(interface)),
        ClimaCalibrate.exeflags(interface),
    )

    # --- Collect forward-model outputs G (N_obs x N_ens) and analyze --------
    g_ens = ClimaCalibrate.observation_map(interface, ITERATION)
    JLD2.save_object(joinpath(SWEEP_OUTPUT_DIR, "sweep_g_ensemble.jld2"), g_ens)

    analyze_sweep(
        obs_series,
        g_ens,
        constrained_samples,
        PD.get_name(SWEEP_PRIORS),
        is_replicate,
        ITERATION,
        SWEEP_OUTPUT_DIR,
    )
end
