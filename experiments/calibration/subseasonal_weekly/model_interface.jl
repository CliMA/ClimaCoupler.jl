import ClimaCoupler

# Include our run_calibration.jl first to define CALIBRATE_CONFIG
include(joinpath(@__DIR__, "run_calibration.jl"))

# Reuse the forward_model from subseasonal pipeline
include(
    joinpath(
        pkgdir(ClimaCoupler),
        "experiments",
        "calibration",
        "subseasonal",
        "model_interface.jl",
    ),
)
