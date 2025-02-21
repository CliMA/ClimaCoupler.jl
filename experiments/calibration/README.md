# ClimaCoupler Calibration Experiments

This folder contains a trivial perfect-model calibration of the atmosphere coupled with the bucket model.
The calibration uses 30-day and lat/lon averages of top-of-atmosphere shortwave 
radiation to calibrate the `total_solar_irradiance` parameter in a perfect model setting. 
The current run script uses the `ClimaCalibrate.SlurmManager` to add Slurm workers
 which run each ensemble member in parallel.

To run this calibration on a Slurm cluster, ensure that `run_calibration.sh` is 
configured for your cluster and run `sbatch run_calibration.sh`. The output will 
be generated in `experiments/calibration/output`.

Components:
- run_calibration.sh: SBATCH script used to instantiate the project and run the calibration on a Slurm cluster.
- run_calibration.jl: Julia script for the overall calibration and postprocessing. Contains the expriment configuration, such as ensemble size and number of iterations.
- model_interface.jl: Contains `forward_model`, the function that gets run during calibration. This basically just uses the `setup_run` function.
- model_config.yml: Contains the configuration for the coupler
- 