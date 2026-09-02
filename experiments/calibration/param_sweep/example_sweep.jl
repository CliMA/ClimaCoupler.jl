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

# Per parameter: a vector sweeps exactly those values, a prior sweeps evenly
# spaced quantiles of the prior (3 by default: the 25th, 50th, and 75th
# percentiles; set prior_points to change), and a scalar pins the value in every
# member. The members are every combination of the swept values, here
# 3 x 3 x 1 = 9 members.
parameters = Dict(
    "rain_autoconversion_timescale" => [900.0, 1800.0, 7200.0],
    "cloud_liquid_water_specific_humidity_autoconversion_threshold" =>
        PD.constrained_gaussian("q_liq_threshold", 5e-4, 3e-4, 0.0, 1.5e-3),
    "snow_autoconversion_timescale" => 1800.0,
)

# Options for ClimaCalibrate.add_workers. Optional, these are the defaults.
# `time` is the worker walltime in minutes. Anything add_workers does not
# consume goes to the scheduler, for example `q = "preempt"`.
worker_options = (; device = :gpu, time = 240)

# Number of workers, one member each. Optional, the default is one per member.
# Fewer workers run the members in waves, so pick a divisor of the member count:
# with 9 members, 9 or 3 workers keep every worker busy, 8 leaves one member to
# run alone in a second wave and roughly doubles the wall time.
n_workers = 9

# Raise the guard on the member count if a large sweep is what you want.
# Optional, the default is 64.
max_members = 64
