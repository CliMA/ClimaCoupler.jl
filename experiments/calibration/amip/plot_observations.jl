# Round-trip check of the generated observations: load observation_vec.jld2,
# reconstruct OutputVars from the EKP Observation objects with
# ClimaCalibrate.ObservationRecipe.reconstruct_vars, and plot them on the
# globe with ClimaAnalysis. This verifies what EKP will actually be shown -
# masking, coarsening, units and orientation - after all preprocessing, not
# what the preprocessing intended. Run it after generate_observations.jl and
# eyeball observation_check.png before submitting anything.
#
# Usage, from the repository root (output_dir = the run's output directory):
#
#   julia --project=experiments/AMIP \
#       experiments/calibration/amip/plot_observations.jl <output_dir>

import ClimaAnalysis
import ClimaCalibrate
import EnsembleKalmanProcesses as EKP
import JLD2
import Statistics
import CairoMakie, GeoMakie
import Printf: @sprintf

output_dir =
    length(ARGS) >= 1 ? ARGS[1] : error("Usage: plot_observations.jl <output_dir>")

obs_vec = JLD2.load_object(joinpath(output_dir, "observation_vec.jld2"))
@info "Loaded $(length(obs_vec)) observations from $output_dir"

# All samples share preprocessing; the first is the calibration target.
vars = ClimaCalibrate.ObservationRecipe.reconstruct_vars(obs_vec[1])

fig = CairoMakie.Figure(size = (1200, 600 * length(vars)))
for (i, var) in enumerate(vars)
    # Reconstructed vars keep a singleton time dimension (the sample date);
    # the on-globe plotters require strictly 2D, so slice it off.
    if ClimaAnalysis.has_time(var)
        var = ClimaAnalysis.slice(var, time = ClimaAnalysis.times(var)[end])
    end
    sn = ClimaAnalysis.short_name(var)
    data = filter(isfinite, vec(var.data))
    n_nan = length(var.data) - length(data)
    stats = @sprintf(
        "%s: %d points (%d masked, %.1f%%)  min %.1f  mean %.1f  max %.1f %s",
        sn,
        length(data),
        n_nan,
        100 * n_nan / length(var.data),
        minimum(data),
        Statistics.mean(data),
        maximum(data),
        ClimaAnalysis.units(var),
    )
    @info stats
    ClimaAnalysis.Visualize.heatmap2D_on_globe!(
        fig[i, 1],
        var;
        more_kwargs = Dict(
            :axis => Dict(:title => stats),
            :coast => Dict(:color => :black),
        ),
    )
end
figpath = joinpath(output_dir, "observation_check.png")
CairoMakie.save(figpath, fig)
@info "Saved $figpath"
