## cloud droplet number concentration calibration

This folder contains materials needed to run a coupled calibration.

Before running, you will need to add two custom branches for packages to the `experiments/calibration` Project:
- `add EnsembleKalmanProcesses#orad/utki`
- `add ClimaAtmos#main`

To run from the coupler root dir: `sbatch experiments/calibration/cld_eff_rad/run_calibration.sh`

Files:
- `run_calibration.sh`: Top-level slurm script to be run from the ClimaCoupler.jl root directory, specifies resources and runs `run_calibration.jl`
- `run_calibration.jl`: Julia script which runs the calibration. Stores the experiment configuration (prior, number of iterations).
- `generate_observations.jl`: Generates observations for the calibration. Requires a diagnostic `OutputVar` with the correct resolution to resample to. Observations are saved to `observations.jld2`.
- `observation_map.jl`: Contains function `observation_map` to process a full ensemble of model outputs into a matrix for the EKP update step.
- `model_interface.jl`: Contains function `forward_model` to run the coupled model with `model_config.yml`.
- `model_config.yml`: YAML file storing coupler configuration.

Configuration:
- EKP: `EKP.TransformUnscented`. This enables us to have a large observation space and few ensemble members per parameter. This is an experimental `Process`. If there are issues, switch back to `EKP.TransformInversion` by uncommmenting the lines at the bottom of `run_calibration.jl`.
- Observations: one year of cloud-radiative effect data (`rsutcs - rsut`) from CERES (via the `radiation_obs` artifact). Covariance matrix is a diagonal with the temporal variance repeated 12 times.
- Prior: One prior for each of the parameters in [`ClimaAtmos.aerosol_ml_parameters`](https://github.com/CliMA/ClimaAtmos.jl/blob/main/src/parameters/create_parameters.jl#L151-L163). The standard deviation and mean likely need to be changed.
- Ensemble size: decided by UTKI, likely 15 (`2p+1` where p is the number of parameters).
- Iterations: 10. We likely don't need this many, but it is better to start too high

