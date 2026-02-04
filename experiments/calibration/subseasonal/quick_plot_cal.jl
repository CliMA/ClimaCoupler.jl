# quick_plot_cal.jl - Publication-quality calibration visualization
# Load EKP first so types are recognized when loading JLD2 files
import EnsembleKalmanProcesses as EKP
import EnsembleKalmanProcesses.ParameterDistributions as PD
using Makie
using CairoMakie
using JLD2
using Statistics
using KernelDensity

# ============================================================================
# Configuration
# ============================================================================
output_dir = "/glade/derecho/scratch/cchristo/calibration/exp25"
exp_suffix = basename(output_dir)
last_iter = 4  # adjust to your final iteration (displays as 0 to last_iter-1)

# Create output folder for plots
plots_dir = "calibration_plots"
mkpath(plots_dir)

# ============================================================================
# Publication-quality theme settings
# ============================================================================
publication_theme = Theme(
    fontsize = 14,
    Axis = (
        xlabelsize = 16,
        ylabelsize = 16,
        ylabelfont = :bold,
        titlesize = 14,
        xticklabelsize = 12,
        yticklabelsize = 12,
        spinewidth = 1.5,
        xgridvisible = true,
        ygridvisible = true,
        xgridcolor = (:gray, 0.2),
        ygridcolor = (:gray, 0.2),
        xgridwidth = 1,
        ygridwidth = 1,
    ),
    Lines = (linewidth = 3,),
    Scatter = (markersize = 12,),
    Legend = (framevisible = false, labelsize = 12),
)
set_theme!(publication_theme)

# Color palette (colorblind-friendly)
colors = [:royalblue, :coral, :seagreen, :mediumpurple, :goldenrod, :crimson]

# ============================================================================
# Load data
# ============================================================================
@info "Loading calibration data from $output_dir"
ekp = JLD2.load(joinpath(output_dir, "iteration_00$(last_iter)", "eki_file.jld2"))["single_stored_object"]
prior = JLD2.load(joinpath(output_dir, "iteration_000", "prior.jld2"))["single_stored_object"]

# Get parameter info
n_params = PD.ndims(prior)
param_names = PD.get_name(prior)
@info "Plotting $n_params parameters: $param_names"

# Calculate grid layout
n_cols = ceil(Int, sqrt(n_params))
n_rows = ceil(Int, n_params / n_cols)

# ============================================================================
# Helper function: Format parameter name for display
# ============================================================================
function format_param_name(name::String)
    # Convert snake_case to Title Case with spaces
    formatted = replace(name, "_" => " ")
    return titlecase(formatted)
end

# ============================================================================
# Figure 1: Parameter evolution over iterations (all ensemble members)
# ============================================================================
@info "Creating parameter evolution plot..."
fig1 = Figure(size = (420 * n_cols, 380 * n_rows), figure_padding = 25)

for i in 1:n_params
    row = div(i - 1, n_cols) + 1
    col = mod(i - 1, n_cols) + 1
    
    ax = Axis(
        fig1[row, col],
        xlabel = "Iteration",
        ylabel = "($(Char('a' + i - 1))) $(format_param_name(param_names[i]))",
    )
    
    # Get parameter values across iterations (in physical/constrained space)
    # Note: EKP uses 1-based indexing; iteration 1 = initial ensemble
    n_iters = EKP.get_N_iterations(ekp)
    for iter in 1:n_iters
        ϕ = EKP.get_ϕ(prior, ekp, iter)  # constrained parameters
        param_vals = ϕ[i, :]
        xs = fill(iter - 1, length(param_vals))  # Display as 0-indexed
        scatter!(ax, xs, param_vals, color = (colors[mod(i-1, length(colors)) + 1], 0.6), markersize = 14)
    end
    
    # Connect means with a line
    means = [mean(EKP.get_ϕ(prior, ekp, iter)[i, :]) for iter in 1:n_iters]
    lines!(ax, 0:(n_iters-1), means, color = :black, linewidth = 3.5)
end

save(joinpath(plots_dir, "evolution_$(exp_suffix).png"), fig1, px_per_unit = 3)
@info "Saved: $(plots_dir)/evolution_$(exp_suffix).png"

# ============================================================================
# Figure 2: Mean ± std convergence
# ============================================================================
@info "Creating mean convergence plot..."
fig2 = Figure(size = (420 * n_cols, 380 * n_rows), figure_padding = 25)

for i in 1:n_params
    row = div(i - 1, n_cols) + 1
    col = mod(i - 1, n_cols) + 1
    
    ax = Axis(
        fig2[row, col],
        xlabel = "Iteration",
        ylabel = "($(Char('a' + i - 1))) $(format_param_name(param_names[i]))",
    )
    
    n_iters = EKP.get_N_iterations(ekp)
    iters = collect(0:(n_iters-1))  # Display as 0-indexed
    means = Float64[]
    stds = Float64[]
    
    for iter in 1:n_iters
        ϕ = EKP.get_ϕ(prior, ekp, iter)
        push!(means, mean(ϕ[i, :]))
        push!(stds, std(ϕ[i, :]))
    end
    
    # Shaded uncertainty band
    band!(ax, iters, means .- stds, means .+ stds, 
          color = (colors[mod(i-1, length(colors)) + 1], 0.25))
    
    # Mean line with markers
    lines!(ax, iters, means, color = colors[mod(i-1, length(colors)) + 1], linewidth = 3)
    scatter!(ax, iters, means, color = colors[mod(i-1, length(colors)) + 1], markersize = 12)
end

save(joinpath(plots_dir, "convergence_$(exp_suffix).png"), fig2, px_per_unit = 3)
@info "Saved: $(plots_dir)/convergence_$(exp_suffix).png"

# ============================================================================
# Figure 3: Error convergence
# ============================================================================
@info "Creating error convergence plot..."
fig3 = Figure(size = (700, 500), figure_padding = 30)
ax3 = Axis(
    fig3[1, 1],
    xlabel = "Iteration",
    ylabel = "Normalized Error",
)

errors = EKP.get_error(ekp)
iterations = collect(0:(length(errors) - 1))

lines!(ax3, iterations, errors, color = :royalblue, linewidth = 3)
scatter!(ax3, iterations, errors, color = :royalblue, markersize = 14)

save(joinpath(plots_dir, "error_$(exp_suffix).png"), fig3, px_per_unit = 3)
@info "Saved: $(plots_dir)/error_$(exp_suffix).png"

# ============================================================================
# Figure 4: Prior vs Posterior distributions
# ============================================================================
@info "Creating prior vs posterior distribution plot..."
fig4 = Figure(size = (420 * n_cols, 380 * n_rows), figure_padding = 25)

# Get initial (prior) and final (posterior) ensemble
# Note: EKP uses 1-based indexing; iteration 1 = initial ensemble
n_iters_total = EKP.get_N_iterations(ekp)
ϕ_prior = EKP.get_ϕ(prior, ekp, 1)                # First iteration ensemble (initial)
ϕ_posterior = EKP.get_ϕ(prior, ekp, n_iters_total) # Final iteration ensemble

for i in 1:n_params
    row = div(i - 1, n_cols) + 1
    col = mod(i - 1, n_cols) + 1
    
    ax = Axis(
        fig4[row, col],
        xlabel = "($(Char('a' + i - 1))) $(format_param_name(param_names[i]))",
        ylabel = "Density",
    )
    
    prior_vals = vec(ϕ_prior[i, :])
    posterior_vals = vec(ϕ_posterior[i, :])
    
    # Compute kernel density estimates
    # Handle edge cases where all values might be the same
    prior_range = maximum(prior_vals) - minimum(prior_vals)
    posterior_range = maximum(posterior_vals) - minimum(posterior_vals)
    
    # Use histograms if range is too small for KDE
    if prior_range > 1e-10 && length(unique(prior_vals)) > 3
        kde_prior = kde(prior_vals)
        lines!(ax, kde_prior.x, kde_prior.density, 
               color = :gray, linewidth = 2.5, linestyle = :dash, label = "Initial (iter 1)")
        band!(ax, kde_prior.x, zeros(length(kde_prior.x)), kde_prior.density, 
              color = (:gray, 0.2))
    else
        # Fallback: plot as vertical lines at the value
        vlines!(ax, prior_vals, color = (:gray, 0.5), linewidth = 1, label = "Initial (iter 1)")
    end
    
    if posterior_range > 1e-10 && length(unique(posterior_vals)) > 3
        kde_posterior = kde(posterior_vals)
        lines!(ax, kde_posterior.x, kde_posterior.density, 
               color = :royalblue, linewidth = 3, label = "Final (iter $(last_iter))")
        band!(ax, kde_posterior.x, zeros(length(kde_posterior.x)), kde_posterior.density, 
              color = (:royalblue, 0.3))
    else
        vlines!(ax, posterior_vals, color = (:royalblue, 0.7), linewidth = 2, 
                label = "Final (iter $(last_iter))")
    end
    
    # Add legend only to first panel
    if i == 1
        axislegend(ax, position = :rt)
    end
end

save(joinpath(plots_dir, "distributions_$(exp_suffix).png"), fig4, px_per_unit = 3)
@info "Saved: $(plots_dir)/distributions_$(exp_suffix).png"

# ============================================================================
# Figure 5: Combined prior/posterior comparison (violin-style)
# ============================================================================
@info "Creating combined violin comparison plot..."
fig5 = Figure(size = (max(800, 150 * n_params), 550), figure_padding = 30)
ax5 = Axis(
    fig5[1, 1],
    xlabel = "Parameter",
    ylabel = "Normalized Value",
    xticks = (1:n_params, [format_param_name(p) for p in param_names]),
    xticklabelrotation = π/4,
)

# Normalize each parameter to [0,1] for comparison
for i in 1:n_params
    prior_vals = vec(ϕ_prior[i, :])
    posterior_vals = vec(ϕ_posterior[i, :])
    
    # Normalize based on combined range
    all_vals = vcat(prior_vals, posterior_vals)
    val_min, val_max = extrema(all_vals)
    range_val = val_max - val_min
    
    if range_val > 1e-10
        prior_norm = (prior_vals .- val_min) ./ range_val
        posterior_norm = (posterior_vals .- val_min) ./ range_val
    else
        prior_norm = fill(0.5, length(prior_vals))
        posterior_norm = fill(0.5, length(posterior_vals))
    end
    
    # Offset for prior (left) and posterior (right)
    offset = 0.15
    
    # Plot initial ensemble points
    scatter!(ax5, fill(i - offset, length(prior_norm)), prior_norm,
             color = (:gray, 0.6), markersize = 10, label = i == 1 ? "Initial" : nothing)
    
    # Plot final ensemble points  
    scatter!(ax5, fill(i + offset, length(posterior_norm)), posterior_norm,
             color = (:royalblue, 0.8), markersize = 10, label = i == 1 ? "Final" : nothing)
    
    # Add mean indicators
    hlines!(ax5, [mean(prior_norm)], xmin = i - 0.25, xmax = i - 0.05, 
            color = :gray, linewidth = 2)
    hlines!(ax5, [mean(posterior_norm)], xmin = i + 0.05, xmax = i + 0.25, 
            color = :royalblue, linewidth = 2)
end

axislegend(ax5, position = :rt)

save(joinpath(plots_dir, "comparison_$(exp_suffix).png"), fig5, px_per_unit = 3)
@info "Saved: $(plots_dir)/comparison_$(exp_suffix).png"

# ============================================================================
# Summary
# ============================================================================
@info """
Plots saved to $(plots_dir)/:
  1. evolution_$(exp_suffix).png     - Parameter evolution (all members)
  2. convergence_$(exp_suffix).png   - Mean ± std convergence  
  3. error_$(exp_suffix).png         - Calibration error
  4. distributions_$(exp_suffix).png - Initial vs final distributions
  5. comparison_$(exp_suffix).png    - Normalized comparison
"""
