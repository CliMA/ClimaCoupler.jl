using ClimaAnalysis, ClimaCoupler
using Statistics, Dates, LinearAlgebra
import JLD2
import EnsembleKalmanProcesses as EKP

FULL_CALIBRATION = false

include("generate_observations_helper.jl")

const start_date = DateTime(2010, 8, 1) # DateTime(2009, 12, 1)
const diagnostic_var =
    OutputVar(joinpath(pkgdir(ClimaCoupler), "rsut_1M_average.nc"))

include(joinpath(pkgdir(ClimaCoupler), "experiments/ClimaEarth/leaderboard/data_sources.jl"))
resample(ov) = resampled_as(shift_longitude(ov, -180.0, 180.0), diagnostic_var, dim_names = ["longitude", "latitude"])

rad_and_pr_obs_dict = get_obs_var_dict()
rsut = resample(rad_and_pr_obs_dict["rsut"](start_date))
rsutcs = resample(rad_and_pr_obs_dict["rsutcs"](start_date))
cre = rsutcs - rsut

# Create an EKP.Observation
const days_in_secs = 86_400
const month = 32 * days_in_secs

if FULL_CALIBRATION
    covariance_mat = compute_covariance(cre)
    cre_window = window(cre, "time"; left = 93days_in_secs, right = 426days_in_secs)
    # Observations info: https://clima.github.io/EnsembleKalmanProcesses.jl/dev/observations/
    cre_obs = EKP.Observation(vec(cre_window.data), covariance_mat, "cre")
else
    # Hardcoding to only do one month
    # TODO: Make it easier to change the number of months; not sure how to make it easier
    covariance_mat = compute_covariance(cre, full = false)
    cre_window = window(cre, "time"; left = 1month, right = 1month)
    cre_obs = EKP.Observation(vec(cre_window.data), covariance_mat, "cre")
end

observations = EKP.combine_observations([cre_obs])
observation_series = EKP.ObservationSeries(observations; metadata = cre_window)

JLD2.save_object("experiments/calibration/cld_eff_rad/obs_series_test.jld2", observation_series)
