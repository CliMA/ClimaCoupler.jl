# plot_calibration.jl
# Load EKP first so types are recognized when loading JLD2 files
import EnsembleKalmanProcesses as EKP
import EnsembleKalmanProcesses.ParameterDistributions as PD
using Makie  # Needed to load EKP's Makie extension
using CairoMakie
using JLD2

output_dir = "/glade/derecho/scratch/cchristo/calibration/exp17"

# Extract directory name for suffix
exp_suffix = basename(output_dir)

# Load the final EKP object and prior
last_iter = 4  # adjust to your final iteration
ekp = JLD2.load(joinpath(output_dir, "iteration_00$(last_iter)", "eki_file.jld2"))["single_stored_object"]
prior = JLD2.load(joinpath(output_dir, "iteration_000", "prior.jld2"))["single_stored_object"]

# Get number of parameters
n_params = PD.ndims(prior)
param_names = PD.get_name(prior)
@info "Plotting $n_params parameters: $param_names"

# Calculate grid layout (aim for roughly square grid)
n_cols = ceil(Int, sqrt(n_params))
n_rows = ceil(Int, n_params / n_cols)

# Plot parameter evolution over iterations (all params)
fig = Figure(size = (400 * n_cols, 350 * n_rows))
for i in 1:n_params
    row = div(i - 1, n_cols) + 1
    col = mod(i - 1, n_cols) + 1
    EKP.Visualize.plot_ϕ_over_iters(fig[row, col], ekp, prior, i)
end
save("calibration_results_$(exp_suffix).png", fig)

# Plot mean with std (all params)
fig2 = Figure(size = (400 * n_cols, 350 * n_rows))
for i in 1:n_params
    row = div(i - 1, n_cols) + 1
    col = mod(i - 1, n_cols) + 1
    EKP.Visualize.plot_ϕ_mean_over_iters(fig2[row, col], ekp, prior, i; plot_std = true)
end
save("param_mean_convergence_$(exp_suffix).png", fig2)

# Plot error over iterations including iteration 0 (manual plot)
fig3 = Figure(size = (600, 400))
ax3 = Axis(fig3[1, 1], xlabel = "Iteration", ylabel = "Error", title = "Error over iterations (including iter 0)")
errors = EKP.get_error(ekp)
iterations = 0:(length(errors) - 1)
scatter!(ax3, iterations, errors)
lines!(ax3, iterations, errors)
save("errors_new_$(exp_suffix).png", fig3)

@info "Plots saved to current directory with suffix: $exp_suffix"