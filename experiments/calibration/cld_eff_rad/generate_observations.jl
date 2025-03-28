using ClimaAnalysis, ClimaCoupler
using Statistics, Dates, LinearAlgebra
import JLD2
import EnsembleKalmanProcesses as EKP

const start_date = DateTime(2009, 10, 1)
const diagnostic_var =
    OutputVar(joinpath(pkgdir(ClimaCoupler), "experiments/calibration/output/iteration_000/member_005/model_config/output_active/clima_atmos/rsut_1M_average.nc"))

include(joinpath(pkgdir(ClimaCoupler), "experiments/ClimaEarth/leaderboard/data_sources.jl"))
resample(ov) = resampled_as(ov, diagnostic_var, dim_names = ["longitude", "latitude"])
function auto_covariance(output_var)
    lat, lon = size(output_var.data)[1:2]
    covariance = zeros(Float32, lat, lon)
    for i in 1:lat
        for j in 1:lon
            covariance[i, j] = var(output_var.data[i, j, :])
        end
    end
    return covariance
end

rad_and_pr_obs_dict = get_obs_var_dict()
rsut = resample(rad_and_pr_obs_dict["rsut"](start_date))
rsutcs = resample(rad_and_pr_obs_dict["rsutcs"](start_date))
cre = rsutcs - rsut

# Create an EKP.Observation
const days_in_secs = 86_400
cre_window = window(cre, "time"; left = 93days_in_secs, right = 426days_in_secs)
cre_obs = EKP.Observation(vec(cre_window.data), Diagonal(repeat(vec(auto_covariance(cre)), 12)), "cre")

observations = EKP.combine_observations([cre_obs])
JLD2.save_object("experiments/calibration/cld_eff_rad/observations.jld2", observations)
