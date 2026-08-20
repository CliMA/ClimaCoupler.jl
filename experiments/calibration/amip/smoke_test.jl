# Integration gate: run ONE ensemble member at the prior centre for a short
# simulated window and check it does not blow up.
#
# WHY THIS EXISTS. preflight.jl validates configuration, observations, run-
# directory hygiene and parameter wiring, but it never integrates the model,
# so it cannot see numerical instability. On 2026-08-14 the micro_edmf
# calibration passed preflight 25/25 and then lost ALL 17 members to NaN and
# CUDA.KernelException 18-35 minutes in - including the centre member, meaning
# the prior mean itself was an unstable model state. A single short member on
# the develop queue would have caught it in under an hour instead of burning a
# full 17-worker launch. Run this after preflight and before any ensemble.
#
# WHAT IT DOES. Rebuilds the config's CalibrateConfig with a private output
# directory and a short `extend`, writes the prior centre as member 1 of
# iteration 1, and calls the ordinary forward_model. Everything downstream of
# that - TOML merge, `_E#` splice, coupled-simulation construction, time
# stepping - is the real code path, not a mock.
#
# Usage (see smoke_test.sh for the PBS wrapper):
#
#   CALIBRATION_CONFIG=.../config/rlut_pigroups.jl \
#   SMOKE_DIR=/glade/derecho/scratch/$USER/pigroups_smoke \
#   SMOKE_DAYS=3 \
#   julia --project=experiments/AMIP experiments/calibration/amip/smoke_test.jl
#
# SMOKE_DAYS is the window graded AFTER spinup, so the total integration is
# spinup + SMOKE_DAYS. With the default 7-day spinup, SMOKE_DAYS=3 gives 10
# simulated days - long enough to get past initial-condition shock and into
# the regime where the micro_edmf blow-ups appeared, short enough to fit the
# develop queue.

ENV["CLIMACOMMS_CONTEXT"] = get(ENV, "CLIMACOMMS_CONTEXT", "SINGLETON")

import Dates
import TOML
import ClimaCalibrate
import ClimaCoupler
import ClimaCoupler: CalibrationTools
import EnsembleKalmanProcesses as EKP
import EnsembleKalmanProcesses.ParameterDistributions as PD
import Statistics

# run_calibration.jl defines the model interface and includes the config; we
# only need the interface and the config's globals, not the EKP loop.
include(joinpath(@__DIR__, "model_interface.jl"))

config_path = get(ENV, "CALIBRATION_CONFIG") do
    error("Set CALIBRATION_CONFIG to the calibration config to smoke test")
end
@info "Smoke testing $config_path"
include(config_path)

smoke_dir = get(ENV, "SMOKE_DIR") do
    error("Set SMOKE_DIR to a scratch directory for the smoke test output")
end
smoke_days = parse(Int, get(ENV, "SMOKE_DAYS", "3"))

cfg = CALIBRATE_CONFIG
smoke_config = CalibrationTools.CalibrateConfig(;
    config_file = cfg.config_file,
    short_names = cfg.short_names,
    minibatch_size = 1,
    n_iterations = 1,
    # Only the first entry is used (forward_model indexes [iter] with iter=1),
    # but keep two so any +1 indexing elsewhere stays in bounds.
    sample_date_ranges = cfg.sample_date_ranges[1:min(2, end)],
    extend = Dates.Day(smoke_days),
    spinup = cfg.spinup,
    output_dir = smoke_dir,
    rng_seed = cfg.rng_seed,
)

start_date = first(smoke_config.sample_date_ranges[1]) - smoke_config.spinup
end_date = last(smoke_config.sample_date_ranges[1]) + smoke_config.extend
@info "Smoke window" start_date end_date total_days = Dates.value(Dates.Day(end_date - start_date))

# Prior centre in constrained (physical) space - the same point member 1 of a
# TransformUnscented ensemble runs. Written under the SAMPLED names; the `_E#`
# splice happens inside forward_model, exactly as it will for a real member.
names = PD.get_name(PRIORS)
u0 = Statistics.mean(PRIORS)
u0 = u0 isa Number ? [u0] : collect(u0)
centre = PD.transform_unconstrained_to_constrained(PRIORS, u0)
centre = centre isa Number ? [centre] : collect(centre)

param_path = ClimaCalibrate.parameter_path(smoke_dir, 1, 1)
mkpath(dirname(param_path))
open(param_path, "w") do io
    TOML.print(
        io,
        Dict(n => Dict("value" => Float64(v), "type" => "float") for
             (n, v) in zip(names, centre)),
    )
end
@info "Prior centre written to $param_path" parameters = Dict(zip(names, centre))

ClimaCalibrate.forward_model(CouplerModelInterface(smoke_config), 1, 1)

# forward_model throws on failure, so reaching here means the member built a
# coupled simulation and integrated the whole window without a NaN abort.
@info "SMOKE TEST PASSED: $(Dates.value(Dates.Day(end_date - start_date))) " *
      "simulated days at the prior centre, no NaN"
