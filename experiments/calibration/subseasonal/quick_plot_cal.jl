# plot_calibration.jl
# Load EKP first so types are recognized when loading JLD2 files
import EnsembleKalmanProcesses as EKP
import EnsembleKalmanProcesses.ParameterDistributions as PD
using Makie  # Needed to load EKP's Makie extension
using CairoMakie
using JLD2

output_dir = "/glade/derecho/scratch/cchristo/calibration/exp9"

# Extract directory name for suffix
exp_suffix = basename(output_dir)

# Load the final EKP object and prior
last_iter = 3  # adjust to your final iteration
ekp = JLD2.load(joinpath(output_dir, "iteration_00$(last_iter)", "eki_file.jld2"))["single_stored_object"]
prior = JLD2.load(joinpath(output_dir, "iteration_000", "prior.jld2"))["single_stored_object"]

# Plot parameter evolution over iterations
fig = Figure(size = (1000, 400))
EKP.Visualize.plot_ϕ_over_iters(fig[1, 1], ekp, prior, 1)  # dim 1 = entr_inv_tau
EKP.Visualize.plot_error_over_iters(fig[1, 2], ekp)
# save(joinpath(output_dir, "calibration_results.png"), fig)
save("calibration_results_$(exp_suffix).png", fig)

# Plot mean with std
fig2 = Figure(size = (600, 400))
EKP.Visualize.plot_ϕ_mean_over_iters(fig2[1, 1], ekp, prior, 1; plot_std = true)
# save(joinpath(output_dir, "param_mean_convergence.png"), fig2)
save("param_mean_convergence_$(exp_suffix).png", fig2)

# Plot error over iterations including iteration 0 (manual plot)
fig3 = Figure(size = (600, 400))
ax3 = Axis(fig3[1, 1], xlabel = "Iteration", ylabel = "Error", title = "Error over iterations (including iter 0)")
errors = EKP.get_error(ekp)
iterations = 0:(length(errors) - 1)
scatter!(ax3, iterations, errors)
lines!(ax3, iterations, errors)
save("errors_new_$(exp_suffix).png", fig3)

@info "Plots saved to $output_dir"