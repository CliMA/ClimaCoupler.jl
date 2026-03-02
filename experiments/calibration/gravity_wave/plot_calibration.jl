#=
plot_calibration.jl

Comprehensive plotting script for calibration results.
Standalone script - does not require the full calibration environment.

Usage:
    julia --project=experiments/ClimaEarth plot_calibration.jl [output_dir]

If output_dir is not provided, defaults to "output/gw_calibration"
=#

# ==========================================================================
# DEPENDENCIES
# ==========================================================================
import EnsembleKalmanProcesses as EKP
import EnsembleKalmanProcesses.ParameterDistributions as PD
using Makie  # Needed to load EKP's Makie extension
using CairoMakie
using JLD2
using Statistics
using Printf

# ==========================================================================
# CONFIGURATION
# ==========================================================================
# Parse command line or use default
output_dir = length(ARGS) >= 1 ? ARGS[1] : "output/gw_calibration"
output_dir = abspath(output_dir)

# Create plots subdirectory
plots_dir = joinpath(output_dir, "plots")
mkpath(plots_dir)

# Extract experiment name for file suffixes
exp_suffix = basename(output_dir)

@info "Calibration Plotting Script"
@info "Output directory: $output_dir"
@info "Plots will be saved to: $plots_dir"

# ==========================================================================
# AUTO-DETECT ITERATIONS
# ==========================================================================
"""
    detect_iterations(output_dir)

Auto-detect iteration directories and find the last complete iteration
(one that contains G_ensemble.jld2).

Returns (all_iterations, last_complete_iteration)
"""
function detect_iterations(output_dir)
    # Find all iteration_XXX directories
    iter_dirs = filter(readdir(output_dir)) do name
        startswith(name, "iteration_") && isdir(joinpath(output_dir, name))
    end

    # Extract iteration numbers
    iterations = sort([parse(Int, split(d, "_")[2]) for d in iter_dirs])

    if isempty(iterations)
        error("No iteration directories found in $output_dir")
    end

    # Find last complete iteration (has G_ensemble.jld2)
    last_complete = 0
    for iter in reverse(iterations)
        iter_path = joinpath(output_dir, @sprintf("iteration_%03d", iter))
        if isfile(joinpath(iter_path, "G_ensemble.jld2"))
            last_complete = iter
            break
        end
    end

    if last_complete == 0
        @warn "No complete iterations found (missing G_ensemble.jld2), using last available"
        last_complete = maximum(iterations)
    end

    @info "Detected $(length(iterations)) iterations, last complete: $last_complete"
    return iterations, last_complete
end

# ==========================================================================
# DATA LOADING UTILITIES
# ==========================================================================
"""
    load_ekp_and_prior(output_dir, iteration)

Load EKP object from specified iteration and prior from iteration 0.
"""
function load_ekp_and_prior(output_dir, iteration)
    iter_str = @sprintf("iteration_%03d", iteration)
    ekp_path = joinpath(output_dir, iter_str, "eki_file.jld2")
    prior_path = joinpath(output_dir, "iteration_000", "prior.jld2")

    ekp = JLD2.load(ekp_path)["single_stored_object"]
    prior = JLD2.load(prior_path)["single_stored_object"]

    return ekp, prior
end

"""
    load_g_ensemble(output_dir, iteration)

Load G_ensemble matrix from specified iteration.
Returns nothing if file doesn't exist.
"""
function load_g_ensemble(output_dir, iteration)
    iter_str = @sprintf("iteration_%03d", iteration)
    g_path = joinpath(output_dir, iter_str, "G_ensemble.jld2")

    if isfile(g_path)
        return JLD2.load(g_path)["single_stored_object"]
    end
    return nothing
end

# ==========================================================================
# MAIN EXECUTION
# ==========================================================================
# Detect iterations
iterations, last_complete = detect_iterations(output_dir)

# Load final EKP and prior
ekp, prior = load_ekp_and_prior(output_dir, last_complete)

# Get parameter info
n_params = PD.ndims(prior)
param_names = PD.get_name(prior)
n_ensemble = EKP.get_N_ens(ekp)

@info "Calibration summary:" n_params param_names n_ensemble last_complete

# Calculate grid layout for multi-parameter plots
n_cols = ceil(Int, sqrt(n_params))
n_rows = ceil(Int, n_params / n_cols)

# ==========================================================================
# PLOT 1: Parameter Evolution (all ensemble members)
# ==========================================================================
@info "Generating Plot 1: Parameter evolution over iterations..."

fig1 = Figure(size = (400 * n_cols, 350 * n_rows))
for i in 1:n_params
    row = div(i - 1, n_cols) + 1
    col = mod(i - 1, n_cols) + 1
    EKP.Visualize.plot_ϕ_over_iters(fig1[row, col], ekp, prior, i)
end
save(joinpath(plots_dir, "01_parameter_evolution_$(exp_suffix).png"), fig1)
@info "Saved: 01_parameter_evolution_$(exp_suffix).png"

# ==========================================================================
# PLOT 2: Parameter Mean Convergence with Std
# ==========================================================================
@info "Generating Plot 2: Parameter mean convergence with std..."

fig2 = Figure(size = (400 * n_cols, 350 * n_rows))
for i in 1:n_params
    row = div(i - 1, n_cols) + 1
    col = mod(i - 1, n_cols) + 1
    EKP.Visualize.plot_ϕ_mean_over_iters(fig2[row, col], ekp, prior, i; plot_std = true)
end
save(joinpath(plots_dir, "02_param_mean_convergence_$(exp_suffix).png"), fig2)
@info "Saved: 02_param_mean_convergence_$(exp_suffix).png"

# ==========================================================================
# PLOT 3: Error/Loss Over Iterations
# ==========================================================================
@info "Generating Plot 3: Error over iterations..."

fig3 = Figure(size = (800, 500))

# Left: scatter + line plot of errors
ax3a = Axis(fig3[1, 1],
    xlabel = "Iteration",
    ylabel = "Error (Loss)",
    title = "Loss Convergence")
errors = EKP.get_error(ekp)
iter_indices = 0:(length(errors) - 1)
scatter!(ax3a, iter_indices, errors, markersize = 12, color = :blue)
lines!(ax3a, iter_indices, errors, color = :blue, linewidth = 2)

# Right: log-scale error
ax3b = Axis(fig3[1, 2],
    xlabel = "Iteration",
    ylabel = "Log₁₀(Error)",
    title = "Loss Convergence (Log Scale)")
log_errors = log10.(max.(errors, 1e-10))
scatter!(ax3b, iter_indices, log_errors, markersize = 12, color = :red)
lines!(ax3b, iter_indices, log_errors, color = :red, linewidth = 2)

save(joinpath(plots_dir, "03_error_convergence_$(exp_suffix).png"), fig3)
@info "Saved: 03_error_convergence_$(exp_suffix).png"

# ==========================================================================
# PLOT 4: Ensemble Spread (Variance) Over Iterations
# ==========================================================================
@info "Generating Plot 4: Ensemble spread over iterations..."

fig4 = Figure(size = (500 * n_cols, 400 * n_rows))

# Get constrained parameters for all iterations
ϕ_all = EKP.get_ϕ(prior, ekp)  # List of matrices, one per iteration

for i in 1:n_params
    row = div(i - 1, n_cols) + 1
    col = mod(i - 1, n_cols) + 1

    ax = Axis(fig4[row, col],
        xlabel = "Iteration",
        ylabel = "Std Dev (constrained)",
        title = "$(param_names[i]) Spread")

    # Calculate std for each iteration
    stds = [std(ϕ[i, :]) for ϕ in ϕ_all]
    iter_idx = 0:(length(stds) - 1)

    barplot!(ax, collect(iter_idx), stds, color = :steelblue)

    # Add convergence percentage text
    if length(stds) > 1 && stds[1] > 0
        reduction = (1 - stds[end] / stds[1]) * 100
        text!(ax, length(stds) - 1, stds[end],
            text = @sprintf("%.0f%% reduction", reduction),
            align = (:center, :bottom),
            fontsize = 10)
    end
end

save(joinpath(plots_dir, "04_ensemble_spread_$(exp_suffix).png"), fig4)
@info "Saved: 04_ensemble_spread_$(exp_suffix).png"

# ==========================================================================
# PLOT 5: Prior vs Posterior Comparison
# ==========================================================================
@info "Generating Plot 5: Prior vs posterior comparison..."

fig5 = Figure(size = (500 * n_cols, 450 * n_rows))

ϕ_initial = ϕ_all[1]   # First iteration (prior samples)
ϕ_final = ϕ_all[end]   # Final iteration (posterior samples)

for i in 1:n_params
    row = div(i - 1, n_cols) + 1
    col = mod(i - 1, n_cols) + 1

    ax = Axis(fig5[row, col],
        xlabel = param_names[i],
        ylabel = "Density",
        title = "$(param_names[i]): Prior vs Posterior")

    # Sample values from prior and posterior
    prior_samples = ϕ_initial[i, :]
    posterior_samples = ϕ_final[i, :]

    # Determine range for histogram
    all_samples = vcat(prior_samples, posterior_samples)
    x_min, x_max = extrema(all_samples)
    x_range = x_max - x_min
    x_min -= 0.1 * x_range
    x_max += 0.1 * x_range

    # Plot histograms
    hist!(ax, prior_samples,
        bins = 15,
        normalization = :pdf,
        color = (:blue, 0.4),
        label = "Prior (iter 0)")
    hist!(ax, posterior_samples,
        bins = 15,
        normalization = :pdf,
        color = (:red, 0.6),
        label = "Posterior (final)")

    # Add mean lines
    vlines!(ax, [mean(prior_samples)], color = :blue, linestyle = :dash, linewidth = 2)
    vlines!(ax, [mean(posterior_samples)], color = :red, linestyle = :dash, linewidth = 2)

    xlims!(ax, x_min, x_max)
    axislegend(ax, position = :rt)
end

save(joinpath(plots_dir, "05_prior_vs_posterior_$(exp_suffix).png"), fig5)
@info "Saved: 05_prior_vs_posterior_$(exp_suffix).png"

# ==========================================================================
# PLOT 6: Parameter Correlation Matrix (if multiple parameters)
# ==========================================================================
if n_params > 1
    @info "Generating Plot 6: Parameter correlation matrix..."

    fig6 = Figure(size = (400 + 150 * n_params, 350 + 150 * n_params))

    ϕ_final_mat = ϕ_all[end]

    # Compute correlation matrix
    cor_matrix = cor(ϕ_final_mat')

    ax6 = Axis(fig6[1, 1],
        xlabel = "Parameter",
        ylabel = "Parameter",
        title = "Final Iteration Parameter Correlations",
        xticks = (1:n_params, param_names),
        yticks = (1:n_params, param_names),
        xticklabelrotation = π/4)

    hm = heatmap!(ax6, 1:n_params, 1:n_params, cor_matrix,
        colormap = :RdBu,
        colorrange = (-1, 1))
    Colorbar(fig6[1, 2], hm, label = "Correlation")

    # Add correlation values as text
    for i in 1:n_params, j in 1:n_params
        text!(ax6, i, j,
            text = @sprintf("%.2f", cor_matrix[i, j]),
            align = (:center, :center),
            fontsize = 12,
            color = abs(cor_matrix[i, j]) > 0.5 ? :white : :black)
    end

    save(joinpath(plots_dir, "06_parameter_correlation_$(exp_suffix).png"), fig6)
    @info "Saved: 06_parameter_correlation_$(exp_suffix).png"
else
    @info "Skipping Plot 6: Only one parameter, no correlation matrix needed"
end

# ==========================================================================
# PLOT 7: G Ensemble Spread Over Iterations
# ==========================================================================
@info "Generating Plot 7: G ensemble spread..."

fig7 = Figure(size = (800, 500))

# Collect G ensemble statistics across iterations
g_means = Float64[]
g_stds = Float64[]
g_spreads = Float64[]  # Mean spread across members

for iter in 0:last_complete
    g_ens = load_g_ensemble(output_dir, iter)
    if !isnothing(g_ens)
        # Mean of all G values
        push!(g_means, mean(g_ens))
        # Std across all values
        push!(g_stds, std(vec(g_ens)))
        # Mean spread per output dimension
        push!(g_spreads, mean([std(g_ens[i, :]) for i in 1:size(g_ens, 1)]))
    end
end

if !isempty(g_means)
    ax7a = Axis(fig7[1, 1],
        xlabel = "Iteration",
        ylabel = "G Ensemble Mean",
        title = "Forward Model Output Mean")
    iter_idx = 0:(length(g_means) - 1)
    scatter!(ax7a, collect(iter_idx), g_means, markersize = 12)
    lines!(ax7a, collect(iter_idx), g_means, linewidth = 2)

    ax7b = Axis(fig7[1, 2],
        xlabel = "Iteration",
        ylabel = "G Ensemble Spread (Std)",
        title = "Forward Model Output Spread")
    scatter!(ax7b, collect(iter_idx), g_spreads, markersize = 12, color = :orange)
    lines!(ax7b, collect(iter_idx), g_spreads, linewidth = 2, color = :orange)

    save(joinpath(plots_dir, "07_g_ensemble_convergence_$(exp_suffix).png"), fig7)
    @info "Saved: 07_g_ensemble_convergence_$(exp_suffix).png"
else
    @warn "No G ensemble data available for Plot 7"
end

# ==========================================================================
# PLOT 8: Ensemble Member Trajectories
# ==========================================================================
@info "Generating Plot 8: Ensemble member trajectories..."

fig8 = Figure(size = (500 * n_cols, 400 * n_rows))

n_iters = length(ϕ_all)

for i in 1:n_params
    row = div(i - 1, n_cols) + 1
    col = mod(i - 1, n_cols) + 1

    ax = Axis(fig8[row, col],
        xlabel = "Iteration",
        ylabel = "Parameter Value (constrained)",
        title = "$(param_names[i]) Trajectories")

    # Plot each ensemble member trajectory
    for m in 1:n_ensemble
        trajectory = [ϕ_all[iter][i, m] for iter in 1:n_iters]
        lines!(ax, 0:(n_iters-1), trajectory,
            color = (:gray, 0.5),
            linewidth = 1)
    end

    # Plot mean trajectory (thicker line)
    mean_trajectory = [mean(ϕ_all[iter][i, :]) for iter in 1:n_iters]
    lines!(ax, 0:(n_iters-1), mean_trajectory,
        color = :red,
        linewidth = 3,
        label = "Mean")

    # Add final mean value annotation
    text!(ax, n_iters - 1, mean_trajectory[end],
        text = @sprintf("%.4f", mean_trajectory[end]),
        align = (:left, :center),
        fontsize = 10,
        color = :red)
end

save(joinpath(plots_dir, "08_member_trajectories_$(exp_suffix).png"), fig8)
@info "Saved: 08_member_trajectories_$(exp_suffix).png"

# ==========================================================================
# SUMMARY
# ==========================================================================
@info """
=====================================
Calibration Plotting Complete!
=====================================
Output directory: $output_dir
Plots saved to: $plots_dir

Generated plots:
  01_parameter_evolution_$(exp_suffix).png     - All ensemble member parameter values
  02_param_mean_convergence_$(exp_suffix).png  - Mean with standard deviation bands
  03_error_convergence_$(exp_suffix).png       - Loss/error over iterations
  04_ensemble_spread_$(exp_suffix).png         - Variance collapse visualization
  05_prior_vs_posterior_$(exp_suffix).png      - Prior vs posterior distributions
  $(n_params > 1 ? "06_parameter_correlation_$(exp_suffix).png   - Parameter correlation matrix" : "(skipped correlation - single param)")
  07_g_ensemble_convergence_$(exp_suffix).png  - Forward model output convergence
  08_member_trajectories_$(exp_suffix).png     - Individual member trajectories

Final parameter estimates (mean +/- std):
"""

ϕ_final = ϕ_all[end]
for i in 1:n_params
    μ = mean(ϕ_final[i, :])
    σ = std(ϕ_final[i, :])
    @info "  $(param_names[i]): $(@sprintf("%.6f", μ)) +/- $(@sprintf("%.6f", σ))"
end
