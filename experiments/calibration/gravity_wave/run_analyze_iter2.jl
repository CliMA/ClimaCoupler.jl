# Script to manually run analyze_iteration for iteration 2
using Dates
import ClimaCalibrate
import ClimaAnalysis
import ClimaCoupler
import EnsembleKalmanProcesses as EKP
import EnsembleKalmanProcesses.ParameterDistributions as PD
import JLD2

include(joinpath(pkgdir(ClimaCoupler), "experiments", "calibration", "api.jl"))
include(joinpath(pkgdir(ClimaCoupler), "experiments/calibration/gravity_wave/observation_map.jl"))

years = 2010:2012
sample_date_ranges = [(DateTime(yr, 2, 1), DateTime(yr, 2, 1)) for yr in years]
const CALIBRATE_CONFIG = CalibrateConfig(;
    config_file = joinpath(pkgdir(ClimaCoupler), "config/amip_configs/amip_land.yml"),
    short_names = ["ta", "ua", "va"],
    minibatch_size = 1,
    n_iterations = 3,
    sample_date_ranges,
    extend = Dates.Month(1),
    spinup = Dates.Month(1),
    output_dir = "output/gravity_wave",
    rng_seed = 42,
)

# Load the saved EKP and prior
iteration = 2
output_dir = CALIBRATE_CONFIG.output_dir

ekp = JLD2.load_object(ClimaCalibrate.ekp_path(output_dir, iteration))
g_ensemble = JLD2.load_object(joinpath(ClimaCalibrate.path_to_iteration(output_dir, iteration), "G_ensemble.jld2"))

# Check how many iterations are in the EKP
@info "EKP state before update:" n_iterations=EKP.get_N_iterations(ekp)

# Manually update EKP with iteration 2's G_ensemble (this was never done due to crash)
EKP.update_ensemble!(ekp, g_ensemble)
@info "EKP state after update:" n_iterations=EKP.get_N_iterations(ekp)

# Save the updated EKP
JLD2.save_object(ClimaCalibrate.ekp_path(output_dir, iteration), ekp)
@info "Saved updated EKP"

# Reconstruct prior
priors = [PD.constrained_gaussian("nogw_Bt_0", 0.0043, 0.002, 0.001, 0.01)]
prior = EKP.combine_distributions(priors)

@info "Running analyze_iteration for iteration $iteration..."
ClimaCalibrate.analyze_iteration(ekp, g_ensemble, prior, output_dir, iteration)
@info "Done! Check $(ClimaCalibrate.path_to_iteration(output_dir, iteration)) for plots"
