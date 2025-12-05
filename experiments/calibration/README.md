# ClimaCoupler Calibration Experiments

This folder contains coupled calibration experiment scripts.
The `coarse_amip` calibration uses the October monthly surface fluxes to calibrate the inverse entrainment timescale.

Components:
- generate_observations.jl: Script to generate observations for the calibration pipeline
- run_calibration.jl: Script for the overall calibration and postprocessing. Contains the expriment configuration, such as ensemble size and number of iterations.
- model_interface.jl: Contains `forward_model`, the function that gets run during calibration. This basically just uses the `setup_run` function.
- observation_map.jl: Contains the loss function and generalized plotting utilities.
