# Zero-GPU pre-check of the G-ensemble date matching: build a SYNTHETIC
# member - an OutputVar with exactly the simulation-side conventions (sim
# grid, time in seconds, monthly slices stamped at period start, start_date
# attribute) for the run's full window - push it through the REAL
# preprocess_sim_vars (mask/coarsen/seasonal-mean hooks included) and the
# REAL GEnsembleBuilder fill against the REAL generated observations, and
# assert that dates match and every observation point is filled.
#
# This exercises the exact obs<->sim alignment that would otherwise first be
# tested hours into iteration 1. Values are dummies (a template slab
# replicated per month); only shapes, dates, units and conventions matter.
#
# Usage (after generate_observations.jl, since the builder needs the
# observation vector and preprocess_sim_vars reads the coverage masks):
#
#   CALIBRATION_CONFIG=<config.jl> julia --project=experiments/AMIP \
#       experiments/calibration/amip/check_g_dates.jl
#
# CHECK_TEMPLATE_SIMDIR may point to any member clima_atmos output directory
# to source the sim-side conventions; defaults to a completed run's member.

ENV["CLIMACOMMS_CONTEXT"] = get(ENV, "CLIMACOMMS_CONTEXT", "SINGLETON")

using Dates
import Random
import ClimaAnalysis
import ClimaCalibrate
import ClimaCalibrate: EnsembleBuilder
import ClimaCoupler
import ClimaCoupler: CalibrationTools
import EnsembleKalmanProcesses as EKP
import JLD2

include(joinpath(@__DIR__, "model_interface.jl"))
config_path = get(ENV, "CALIBRATION_CONFIG") do
    error("Set CALIBRATION_CONFIG")
end
include(config_path)

# --- the real observations and the real builder ------------------------------
(; output_dir, sample_date_ranges, minibatch_size, rng_seed) = CALIBRATE_CONFIG
observation_vector = JLD2.load_object(joinpath(output_dir, "observation_vec.jld2"))
obs_series = EKP.ObservationSeries(
    Dict(
        "observations" => observation_vector,
        "names" => [string(Dates.year(d1)) for (d1, _) in sample_date_ranges],
        "minibatcher" => ClimaCalibrate.minibatcher_over_samples(
            length(observation_vector),
            minibatch_size,
        ),
    ),
)
ekp = EKP.EnsembleKalmanProcess(
    obs_series,
    EKP.TransformUnscented(PRIORS, impose_prior = true);
    rng = Random.MersenneTwister(rng_seed),
    scheduler = EKP.DefaultScheduler(0.1),
)
g_ens_builder = EnsembleBuilder.GEnsembleBuilder(ekp)

# --- synthetic member with sim conventions -----------------------------------
template_dir = get(
    ENV,
    "CHECK_TEMPLATE_SIMDIR",
    "/glade/derecho/scratch/kphan/amip_calibration_rlut_pigroups_ocean_out/" *
    "iteration_001/member_001/amip_calibration_pigroups/output_0000/clima_atmos",
)
template = get(
    ClimaAnalysis.SimDir(template_dir);
    short_name = "rlut",
    reduction = "average",
    period = "1M",
)

(; spinup, extend) = CALIBRATE_CONFIG
run_start = first(sample_date_ranges[1]) - spinup
run_end = last(sample_date_ranges[1]) + extend
# Monthly slices stamped at period START, as ClimaDiagnostics writes them for
# a month-aligned start (verified on real member output).
months = collect(run_start:Dates.Month(1):(run_end - Dates.Day(1)))
times_s = [float(Dates.value(Dates.Second(m - run_start))) for m in months]
@info "Synthetic member" run_start run_end months = Dates.format.(months, "yyyy-mm")

tname = ClimaAnalysis.time_name(template)
ti = template.dim2index[tname]
slab = selectdim(template.data, ti, 1:1)
newdata = cat(fill(Array(slab), length(months))...; dims = ti)
newdims = copy(template.dims)
newdims[tname] = times_s
newattribs = Dict{String, Any}(string(k) => v for (k, v) in template.attributes)
newattribs["start_date"] = Dates.format(run_start, dateformat"yyyy-mm-ddTHH:MM:SS")
# One synthetic var per GRADED variable (a single-var synthetic can only fill
# its own block, so a multi-variable observation would always "fail"). All are
# clones of the same template slab: values are dummies; only shapes, dates,
# units and conventions matter. NOTE the template's units (W m^-2) must match
# every graded variable - true for rlut/rsut/rlutcs/rsutcs/lwcre/swcre; a
# future config grading a different-unit variable needs its own template.
synths = map(CALIBRATE_CONFIG.short_names) do name
    attribs = copy(newattribs)
    attribs["short_name"] = name
    ClimaAnalysis.remake(
        template;
        data = newdata,
        dims = newdims,
        attributes = attribs,
    )
end
@info "Synthetic slice dates as ClimaAnalysis sees them" ClimaAnalysis.dates(synths[1])

# --- the real pipeline --------------------------------------------------------
vars = preprocess_sim_vars(Any[synths...])
@info "After preprocess_sim_vars" ClimaAnalysis.dates(vars[1])

fills = map(vars) do v
    EnsembleBuilder.fill_g_ens_col!(
        g_ens_builder,
        1,
        v;
        checkers = (SequentialIndicesChecker(),),
        verbose = true,
    )
end
ok = all(fills)
missing_names = EnsembleBuilder.missing_short_names(g_ens_builder, 1)
g = EnsembleBuilder.get_g_ensemble(g_ens_builder)
n_filled = count(!isnan, view(g, :, 1))
n_obs = length(EKP.get_obs(observation_vector[1]))

println("=" ^ 72)
if ok && isempty(missing_names) && n_filled == n_obs
    println("G-DATE CHECK PASSED: synthetic member filled $n_filled/$n_obs ",
            "observation points; dates and flattening align")
else
    println("G-DATE CHECK FAILED: fill returned $ok, missing = $missing_names, ",
            "filled $n_filled of $n_obs")
    exit(1)
end
