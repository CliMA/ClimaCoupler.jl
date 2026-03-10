# Estimate climate sensitivity from a pair of AMIP runs: a control
# run and a run with uniformly raised sea surface temperature (+2K SST). The difference
# in time-averaged TOA outgoing fluxes (rsut, rlut) between the two runs is used to
# compute the climate feedback parameter λ = ΔF_net / ΔTs, from which
# ECS ≈ ERF_2xCO2 / λ is derived.
#
# Usage: julia --project=experiments/ClimaEarth/ cess_climate_sensitivity.jl \
#          <ctrl_atmos_dir> <p2k_atmos_dir> <output_dir>

import ClimaAnalysis
import GeoMakie
import CairoMakie
import Dates

if length(ARGS) < 3
    error(
        "Usage: julia cess_climate_sensitivity.jl " *
        "<ctrl_atmos_dir> <p2k_atmos_dir> <output_dir>",
    )
end

ctrl_dir = ClimaAnalysis.SimDir(ARGS[1])
p2k_dir = ClimaAnalysis.SimDir(ARGS[2])
output_dir = ARGS[3]

mkpath(output_dir)

# Remove first 3 months of each run as spinup
const SPINUP_MONTHS = 3
# Radiative forcing for a doubling of CO2, derived from the simplified
# expression in Table 3 of Myhre et al. (1998, doi:10.1029/98GL01908): 
# ΔF = 5.35·ln(C/C₀) which gives 5.35·ln(2) ≈ 3.71 W m⁻²
const ERF_2xCO2 = 3.7  # W m⁻²
# SST perturbation applied in the p2k run
const DELTA_SST = 2.0  # K

# Load a monthly-averaged TOA flux variable from simdir, remove spinup months,
# and return the time-averaged field together with the dates of the averaging period.
function load_and_average(simdir, short_name)
    var = ClimaAnalysis.get(simdir; short_name, reduction = "average", period = "1M")
    spinup_cutoff = SPINUP_MONTHS * 31 * 86400.0
    if ClimaAnalysis.times(var)[end] >= spinup_cutoff
        var = ClimaAnalysis.window(var, "time", left = spinup_cutoff)
    end
    return ClimaAnalysis.average_time(var)
end

@info "Loading TOA fluxes"
rsut_ctrl = load_and_average(ctrl_dir, "rsut")
rlut_ctrl = load_and_average(ctrl_dir, "rlut")
rsut_p2k = load_and_average(p2k_dir, "rsut")
rlut_p2k = load_and_average(p2k_dir, "rlut")
var_dates = ClimaAnalysis.dates(ClimaAnalysis.get(ctrl_dir, "rlut"))

# Compute differences (p2k - ctrl)
delta_rsut = rsut_p2k - rsut_ctrl
delta_rlut = rlut_p2k - rlut_ctrl

# Set metadata on difference vars so plot titles are informative
delta_rsut.attributes["short_name"] = "delta_rsut"
delta_rsut.attributes["long_name"] = "TOA upwelling SW flux change (+2K SST − control)"
delta_rlut.attributes["short_name"] = "delta_rlut"
delta_rlut.attributes["long_name"] = "TOA upwelling LW flux change (+2K SST − control)"

# Compute area-weighted global means of each flux change
delta_rsut_global = ClimaAnalysis.weighted_average_lonlat(delta_rsut).data[]
delta_rlut_global = ClimaAnalysis.weighted_average_lonlat(delta_rlut).data[]

# Net TOA outgoing radiation change: positive means more energy lost to space
delta_toa = delta_rsut_global + delta_rlut_global

@info "Global mean TOA flux changes (p2K − control):"
@info "delta_rsut = $(round(delta_rsut_global, digits = 3)) W m⁻²"
@info "delta_rlut = $(round(delta_rlut_global, digits = 3)) W m⁻²"
@info "delta_toa = $(round(delta_toa, digits = 3)) W m⁻²"

# Climate feedback parameter λ = delta_toa / delta_sst
# ECS ≈ ERF_2xCO2 / λ 
λ = delta_toa / DELTA_SST
ecs_estimate = ERF_2xCO2 / λ

@info "Climate sensitivity estimate:"
@info "λ = delta_toa / delta_sst = $(round(λ, digits = 3)) W m⁻² K⁻¹"
@info "Approximate ECS = ERF_2xCO2 / λ ≈ $(round(ecs_estimate, digits = 2)) K"

# Build figure title from the actual averaging period and spinup
date_fmt = d -> Dates.format(d, "u yyyy")
fig_title =
    "TOA flux change (+2K SST − control) | " *
    "averaging period: $(date_fmt(var_dates[1])) – $(date_fmt(var_dates[end])) | " *
    "spinup: $SPINUP_MONTHS months"

fig = CairoMakie.Figure(; size = (700, 900))
CairoMakie.Label(fig[0, 1], fig_title; tellwidth = false)
ClimaAnalysis.Visualize.plot_bias_on_globe!(fig[1, 1], rsut_p2k, rsut_ctrl)
ClimaAnalysis.Visualize.plot_bias_on_globe!(fig[2, 1], rlut_p2k, rlut_ctrl)
output_path = joinpath(output_dir, "toa_flux_difference.png")
CairoMakie.save(output_path, fig)
@info "Saved plot to $output_path"
