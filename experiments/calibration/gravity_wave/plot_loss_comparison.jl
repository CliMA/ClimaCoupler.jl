#=
plot_loss_comparison.jl

Plot the EKP loss (the actual metric being optimized) and break it down by variable.

Usage:
    julia --project=experiments/ClimaEarth plot_loss_comparison.jl [output_dir]
=#

import EnsembleKalmanProcesses as EKP
import EnsembleKalmanProcesses.ParameterDistributions as PD
using CairoMakie
using JLD2
using Statistics
using Printf
using Dates
using LinearAlgebra

# Override JLD2's default_iotype to avoid Lustre issues
JLD2.default_iotype() = IOStream

# ==========================================================================
# CONFIGURATION
# ==========================================================================
output_dir = length(ARGS) >= 1 ? ARGS[1] : "output/gravity_wave"
output_dir = abspath(output_dir)

plots_dir = joinpath(output_dir, "plots")
mkpath(plots_dir)

exp_suffix = basename(output_dir)

@info "EKP Loss Analysis"
@info "Output directory: $output_dir"

# ==========================================================================
# LOAD DATA
# ==========================================================================

# Auto-detect iterations
function detect_iterations(output_dir)
    iter_dirs = filter(readdir(output_dir)) do name
        startswith(name, "iteration_") && isdir(joinpath(output_dir, name))
    end
    iterations = sort([parse(Int, split(d, "_")[2]) for d in iter_dirs])

    if isempty(iterations)
        error("No iteration directories found in $output_dir")
    end

    # Find last iteration with eki_file.jld2
    last_complete = 0
    for iter in reverse(iterations)
        iter_path = joinpath(output_dir, @sprintf("iteration_%03d", iter))
        if isfile(joinpath(iter_path, "eki_file.jld2"))
            last_complete = iter
            break
        end
    end

    return iterations, last_complete
end

iterations, last_complete = detect_iterations(output_dir)
@info "Found iterations: $iterations, last complete: $last_complete"

# Load EKP from last complete iteration
ekp_path = joinpath(output_dir, @sprintf("iteration_%03d", last_complete), "eki_file.jld2")
prior_path = joinpath(output_dir, "iteration_000", "prior.jld2")

@info "Loading EKP from $ekp_path"
ekp = JLD2.load(ekp_path)["single_stored_object"]
prior = JLD2.load(prior_path)["single_stored_object"]

# ==========================================================================
# EXTRACT LOSS DATA
# ==========================================================================

# Get number of completed iterations in the EKP object
n_iters = length(EKP.get_u(ekp))
@info "EKP has $n_iters iterations of data"

# Get EKP's internally computed error (this is what the constrained_params_and_error.png shows)
ekp_errors = EKP.get_error(ekp)
@info "EKP internal errors: $ekp_errors"

# Get the loss at each iteration (manually computed)
losses = Float64[]
losses_normalized = Float64[]  # Normalized by dim(y) like EKP does
raw_mses = Float64[]

# Get observation series
obs_series = EKP.get_observation_series(ekp)

# Get number of G ensembles actually stored
g_list = EKP.get_g(ekp)
n_g = length(g_list)
@info "Number of G ensembles stored: $n_g"

for i in 1:n_g
    g = g_list[i]  # G ensemble at iteration i
    g_mean = mean(g, dims=2)[:]  # Ensemble mean

    # Get observation for this iteration's minibatch
    y = EKP.get_obs(obs_series, i)  # Observation vector
    dim_y = length(y)

    # Get noise covariance - use build=false to avoid allocating full matrix
    # This returns a Vector{Diagonal} (one per observation in minibatch)
    Γ_vec = EKP.get_obs_noise_cov(obs_series, i; build=false)

    # Concatenate diagonals from each Diagonal matrix
    Γ_diag = vcat([diag(Γ) for Γ in Γ_vec]...)

    # Compute loss element-wise: sum((g - y)^2 / Γ_diag)
    diff = g_mean .- y
    loss = sum(diff .^ 2 ./ Γ_diag)
    push!(losses, loss)

    # Normalized loss (like EKP's avg_rmse): 1/dim(y) * sum(...)
    push!(losses_normalized, loss / dim_y)

    # Also compute raw MSE (unnormalized)
    mse = mean(diff .^ 2)
    push!(raw_mses, mse)

    @info @sprintf("  Iteration %d: Loss = %.4f, Normalized = %.4f, Raw MSE = %.4f", i-1, loss, loss/dim_y, mse)
end

@info """
Comparison with EKP internal errors:
  EKP internal:     $(ekp_errors)
  My normalized:    $(losses_normalized)
  Difference:       $(ekp_errors .- losses_normalized)
"""

# ==========================================================================
# PLOT 1: EKP Internal Error (what constrained_params_and_error.png shows)
# ==========================================================================

fig = Figure(size = (1600, 800))

# Panel 1: EKP's internal error
ax1 = Axis(fig[1, 1],
    xlabel = "Iteration",
    ylabel = "EKP Internal Error",
    title = "EKP get_error() - what the EKP plot shows")

if !isempty(ekp_errors)
    # EKP errors are indexed 1:N, corresponding to iterations 1:N (after updates)
    scatterlines!(ax1, 1:length(ekp_errors), ekp_errors,
        color = :green,
        linewidth = 2,
        markersize = 10,
        label = "EKP internal")

    for i in 2:length(ekp_errors)
        pct_change = 100 * (ekp_errors[i] - ekp_errors[i-1]) / ekp_errors[i-1]
        text!(ax1, i, ekp_errors[i],
            text = @sprintf("%+.1f%%", pct_change),
            fontsize = 10,
            align = (:center, :bottom),
            offset = (0, 5))
    end
end

# ==========================================================================
# PLOT 2: My normalized loss (should match EKP)
# ==========================================================================

ax2 = Axis(fig[1, 2],
    xlabel = "Iteration",
    ylabel = "Loss / dim(y)",
    title = "My Normalized Loss: (1/N) Σ (g-y)²/Γ")

if !isempty(losses_normalized)
    scatterlines!(ax2, 0:(length(losses_normalized)-1), losses_normalized,
        color = :blue,
        linewidth = 2,
        markersize = 10)

    for i in 2:length(losses_normalized)
        pct_change = 100 * (losses_normalized[i] - losses_normalized[i-1]) / losses_normalized[i-1]
        text!(ax2, i-1, losses_normalized[i],
            text = @sprintf("%+.1f%%", pct_change),
            fontsize = 10,
            align = (:center, :bottom),
            offset = (0, 5))
    end
end

# ==========================================================================
# PLOT 3: Raw squared error (without Γ normalization)
# ==========================================================================

ax3 = Axis(fig[1, 3],
    xlabel = "Iteration",
    ylabel = "MSE (unnormalized)",
    title = "Raw MSE: mean((G - y)²)")

# raw_mses already computed in the loss loop above

if !isempty(raw_mses)
    scatterlines!(ax3, 0:(length(raw_mses)-1), raw_mses,
        color = :red,
        linewidth = 2,
        markersize = 10)

    for i in 2:length(raw_mses)
        pct_change = 100 * (raw_mses[i] - raw_mses[i-1]) / raw_mses[i-1]
        text!(ax3, i-1, raw_mses[i],
            text = @sprintf("%+.1f%%", pct_change),
            fontsize = 10,
            align = (:center, :bottom),
            offset = (0, 5))
    end
end

# ==========================================================================
# PLOT 4: Noise covariance diagonal (shows relative weighting)
# ==========================================================================

ax4 = Axis(fig[2, 1:3],
    xlabel = "Observation index",
    ylabel = "Noise std dev (√Γᵢᵢ)",
    title = "Observation uncertainty (determines loss weighting)")

# Get Γ for the first iteration (build=false to avoid full matrix allocation)
Γ_vec = EKP.get_obs_noise_cov(obs_series, 1; build=false)
Γ_diag_full = vcat([diag(Γ) for Γ in Γ_vec]...)
noise_std = sqrt.(Γ_diag_full)

# Plot noise std dev
lines!(ax4, 1:length(noise_std), noise_std, color = :gray)

# Try to annotate regions (if we know the structure)
# Observation vector is likely: [var1_level1, var1_level2, ..., var2_level1, ...]
n_obs = length(noise_std)
@info "Total observation vector length: $n_obs"

# ==========================================================================
# SUMMARY STATISTICS
# ==========================================================================

if !isempty(losses)
    initial_loss = losses[1]
    final_loss = losses[end]
    total_change = 100 * (final_loss - initial_loss) / initial_loss

    @info """
    =====================================
    EKP Loss Summary
    =====================================
    Initial loss (iter 0): $(@sprintf("%.4f", initial_loss))
    Final loss (iter $(length(losses)-1)):   $(@sprintf("%.4f", final_loss))
    Total change: $(@sprintf("%+.2f%%", total_change))

    Raw MSE (unnormalized):
    Initial: $(@sprintf("%.4f", raw_mses[1]))
    Final:   $(@sprintf("%.4f", raw_mses[end]))
    Change:  $(@sprintf("%+.2f%%", 100 * (raw_mses[end] - raw_mses[1]) / raw_mses[1]))
    """
end

# Save figure
output_path = joinpath(plots_dir, "loss_comparison_$(exp_suffix).png")
save(output_path, fig)
@info "Saved: $output_path"
