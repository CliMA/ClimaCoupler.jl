#=
plot_calibration_2panel.jl

Publication-ready two-panel figure: parameter evolution + ensemble spread.

Usage:
    julia --project=experiments/ClimaEarth plot_calibration_2panel.jl [output_dir]
=#

import EnsembleKalmanProcesses as EKP
import EnsembleKalmanProcesses.ParameterDistributions as PD
using CairoMakie
using LaTeXStrings
using JLD2
using Statistics
using Printf

# ==========================================================================
# CONFIGURATION
# ==========================================================================
output_dir = length(ARGS) >= 1 ? ARGS[1] : "output/gw_calibration"
output_dir = abspath(output_dir)

plots_dir = joinpath(output_dir, "plots")
mkpath(plots_dir)

exp_suffix = basename(output_dir)

# Parameter display configuration
# Maps parameter names to (label, unit, scale_factor)
PARAM_CONFIG = Dict(
    "gw_Fs0" => (L"$F_{s0}$", "mPa", 1000.0),
    "ogw_mountain_height_width_exponent" => (L"$\gamma$", "", 1.0),
)
DEFAULT_CONFIG = ("Parameter", "", 1.0)

# ==========================================================================
# AUTO-DETECT ITERATIONS
# ==========================================================================
function detect_iterations(output_dir)
    iter_dirs = filter(readdir(output_dir)) do name
        startswith(name, "iteration_") && isdir(joinpath(output_dir, name))
    end
    iterations = sort([parse(Int, split(d, "_")[2]) for d in iter_dirs])

    if isempty(iterations)
        error("No iteration directories found in $output_dir")
    end

    last_complete = 0
    for iter in reverse(iterations)
        iter_path = joinpath(output_dir, @sprintf("iteration_%03d", iter))
        if isfile(joinpath(iter_path, "G_ensemble.jld2"))
            last_complete = iter
            break
        end
    end

    if last_complete == 0
        last_complete = maximum(iterations)
    end

    return iterations, last_complete
end

# ==========================================================================
# DATA LOADING
# ==========================================================================
function load_ekp_and_prior(output_dir, iteration)
    iter_str = @sprintf("iteration_%03d", iteration)
    ekp_path = joinpath(output_dir, iter_str, "eki_file.jld2")
    prior_path = joinpath(output_dir, "iteration_000", "prior.jld2")

    ekp = JLD2.load(ekp_path)["single_stored_object"]
    prior = JLD2.load(prior_path)["single_stored_object"]

    return ekp, prior
end

# ==========================================================================
# MAIN
# ==========================================================================
iterations, last_complete = detect_iterations(output_dir)
ekp, prior = load_ekp_and_prior(output_dir, last_complete)

n_ensemble = EKP.get_N_ens(ekp)
ϕ_all = EKP.get_ϕ(prior, ekp)
n_iters = length(ϕ_all)

# Get parameter name from prior and look up display config
param_name = PD.get_name(prior)[1]
param_label, param_unit, scale_factor = get(PARAM_CONFIG, param_name, DEFAULT_CONFIG)
ylabel_str = isempty(param_unit) ? param_label : LaTeXString(param_label.s * " ($param_unit)")

@info "Generating two-panel calibration figure..."
@info "Parameter: $param_name, Iterations: $n_iters, Ensemble members: $n_ensemble"

# Create figure with compact height for poster column
fig = Figure(size = (600, 160))

# Member colors
member_colors = [:dodgerblue, :darkorange, :seagreen]

# ==========================================================================
# Panel 1: Parameter Evolution
# ==========================================================================
ax1 = Axis(fig[1, 1],
    xlabel = "Iteration",
    ylabel = ylabel_str)

for m in 1:n_ensemble
    trajectory = [ϕ_all[iter][1, m] * scale_factor for iter in 1:n_iters]
    scatterlines!(ax1, 0:(n_iters-1), trajectory,
        color = member_colors[m],
        linewidth = 1,
        markersize = 8,
        label = "Member $m")
end

# Final mean horizontal line
final_mean = mean(ϕ_all[end][1, :]) * scale_factor
mean_label = isempty(param_unit) ? @sprintf("Mean: %.3f", final_mean) : @sprintf("Mean: %.1f %s", final_mean, param_unit)
hlines!(ax1, [final_mean],
    color = :black,
    linestyle = :dash,
    linewidth = 2,
    label = mean_label)

axislegend(ax1, position = :rt, framevisible = false,
    labelsize = 9, patchsize = (10, 5), rowgap = 0, padding = (2, 2, 2, 2))

# ==========================================================================
# Panel 2: Loss vs Iterations
# ==========================================================================
ax2 = Axis(fig[1, 2],
    xlabel = "Iteration",
    ylabel = "Loss")

errors = EKP.get_error(ekp)
iter_indices = 0:(length(errors) - 1)
scatterlines!(ax2, iter_indices, errors,
    color = :steelblue,
    linewidth = 2,
    markersize = 8)

# ==========================================================================
# Save with tight margins
# ==========================================================================
output_png = joinpath(plots_dir, "calibration_2panel_$(exp_suffix).png")
output_pdf = joinpath(plots_dir, "calibration_2panel_$(exp_suffix).pdf")
save(output_png, fig, px_per_unit = 2)
save(output_pdf, fig)

@info "Saved: $output_png"
@info "Saved: $output_pdf"

# Print prior vs posterior statistics
prior_samples = ϕ_all[1][1, :] * scale_factor
posterior_samples = ϕ_all[end][1, :] * scale_factor
unit_str = isempty(param_unit) ? "" : " $param_unit"

@info """
Prior vs Posterior ($param_name):
  Prior:     $(@sprintf("%.4g", mean(prior_samples))) ± $(@sprintf("%.4g", std(prior_samples)))$unit_str
  Posterior: $(@sprintf("%.4g", mean(posterior_samples))) ± $(@sprintf("%.4g", std(posterior_samples)))$unit_str
"""
