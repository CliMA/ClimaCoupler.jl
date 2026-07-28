# Example parameter sweep. Copy, edit, then from the repository root:
#
#   julia --project=experiments/AMIP \
#       experiments/calibration/param_sweep/run_sweep.jl my_sweep.jl
#
# Two inputs: config_file (which model) and parameters (what to vary).

config_file =
    joinpath(pkgdir(ClimaCoupler), "config", "amip_configs", "amip_calibration.yml")

# Absolute path on scratch. Model output is large, never use home.
output_dir = "/glade/derecho/scratch/CHANGE_ME/amip_my_sweep"

# Diagnostic variables for the response analysis (monthly means from the
# model output). Optional, this is the default.
short_names = ["lwp", "pr", "rsut", "rlut", "tas"]

# Per parameter: a vector sweeps exactly those values, a prior sweeps its
# quantiles (default 3, set prior_points to change), and a scalar pins the
# value in every member. The ensemble is the full factorial across the swept
# parameters, here 3 x 3 = 9 members.
parameters = Dict(
    "rain_autoconversion_timescale" => [900.0, 1800.0, 7200.0],
    "cloud_liquid_water_specific_humidity_autoconversion_threshold" =>
        PD.constrained_gaussian("q_liq_threshold", 5e-4, 3e-4, 0.0, 1.5e-3),
    "snow_autoconversion_timescale" => 1800.0,
)

# Near-copies of the first member that measure the internal-variability
# noise floor.
replicates = 3

# Options for ClimaCalibrate.add_workers. Optional, these are the defaults.
# `time` is the worker walltime in minutes. Anything add_workers does not
# consume goes to the scheduler, so `q = "preempt"` (20% charge) is the cheapest
# way to run: sweep members are independent and resumable, so preemption costs
# the least here.
worker_options = (; device = :gpu, time = 240)

# Number of workers. Optional, the default is one per member. Fewer workers run
# the members in waves.
n_workers = 12

