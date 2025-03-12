using ClimaAnalysis, ClimaCoupler
using Statistics, Dates, LinearAlgebra
import JLD2
import EnsembleKalmanProcesses as EKP

const start_date = DateTime(2000, 3, 1)
const slice_time = 0.0
const diagnostic_var = OutputVar(joinpath(pkgdir(ClimaCoupler), "experiments/calibration/output/model_config/output_active/clima_atmos/rsut_1M_average.nc"))

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
cre = rsut - rsutcs

# Create an EKP.Observation
rsut_slice = slice(rsut; time = slice_time)
rsut_obs = EKP.Observation(
    vec(rsut_slice.data),
    Diagonal(vec(auto_covariance(rsut))),
    "rsut"
)
cre_slice = slice(cre; time = slice_time)
cre_obs = EKP.Observation(
    vec(cre_slice.data),
    Diagonal(vec(auto_covariance(cre))),
    "rsut"
)

observations = EKP.combine_observations([rsut_obs, cre_obs])
JLD2.save_object("experiments/calibration/cld_eff_rad/observations.jld2", observations)
