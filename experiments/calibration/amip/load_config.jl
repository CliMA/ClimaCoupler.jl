# Shared config selection logic. Include this file to get CALIBRATE_CONFIG,
# PRIORS, PRESSURE_LEVELS, NORMALIZATION_STATS_FP, etc.
#
# Set the CALIBRATE_CONFIG environment variable to the path of a config file to
# override the default (config/pressure_levels.jl).

using Dates
import EnsembleKalmanProcesses as EKP
import EnsembleKalmanProcesses.ParameterDistributions as PD
import ClimaCoupler
import ClimaCoupler: CalibrationTools

default_config = abspath(joinpath(@__DIR__, "config", "pressure_levels.jl"))
config_path = abspath(get(ENV, "CALIBRATE_CONFIG", default_config))
isfile(config_path) || error("Calibration config file not found: $config_path")
@info "Using calibration configuration in: $config_path"
include(config_path)

# Make output directory for generating observations and starting the calibration
(; output_dir) = CALIBRATE_CONFIG
isdir(output_dir) || mkdir(output_dir)

# TODO: Test calibration use be renamed to emulate diagnostics and as an
# environment variable
# TODO: Maybe use argparse for both cases
# TODO: Add ClimaGPUBackend as one of the option for testing calibration
