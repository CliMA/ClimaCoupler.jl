import ClimaAnalysis
import Dates
import Test: @test, @testset

include("data_sources.jl")

"""
    test_rmse_thresholds(diagnostics_folder_path, spinup)

Test that the annual RMSE values for specific variables have not increased
beyond acceptable thresholds. The variables tested are:
- pr (precipitation)
- rsut (top of atmosphere outgoing shortwave radiation)
- rsutcs (clear-sky top of atmosphere outgoing shortwave radiation)

The spinup is the number of months to remove from the beginning of the
simulation.

More variables can be added by adding the short name and RMSE pair to the
dictionary returned by `get_rmse_thresholds`.

If this test fails, it indicates a regression in the model's physics, resulting
in a higher RMSE. If this increased RMSE is considered acceptable, then the
thresholds should be updated accordingly.
"""
function test_rmse_thresholds(diagnostics_folder_path, spinup)
    sim_var_dict = get_sim_var_dict(diagnostics_folder_path)
    obs_var_dict = get_obs_var_dict()
    rmse_thresholds = get_rmse_thresholds()

    sim_vars = (sim_var_dict[short_name]() for short_name in keys(rmse_thresholds))
    obs_vars = (
        obs_var_dict[ClimaAnalysis.short_name(sim_var)](sim_var.attributes["start_date"]) for sim_var in sim_vars
    )
    short_names = (ClimaAnalysis.short_name(var) for var in sim_vars)

    rmses = map(sim_vars, obs_vars) do sim_var, obs_var
        # Remove first spin_up_months from simulation
        spinup_cutoff = spinup * 31 * 86400.0
        ClimaAnalysis.times(sim_var)[end] >= spinup_cutoff &&
            (sim_var = ClimaAnalysis.window(sim_var, "time", left = spinup_cutoff))

        obs_var = ClimaAnalysis.resampled_as(obs_var, sim_var)
        obs_var = ClimaAnalysis.average_time(obs_var)
        sim_var = ClimaAnalysis.average_time(sim_var)

        ClimaAnalysis.global_rmse(sim_var, obs_var)
    end

    @testset "RMSE thresholds" begin
        for (short_name, rmse) in zip(short_names, rmses)
            @info "RMSE for $short_name: $rmse"
            @test rmse < rmse_thresholds[short_name]
        end
    end
end

"""
    get_rmse_thresholds()

Return a dictionary mapping short names to maximum acceptable RMSE values.
"""
function get_rmse_thresholds()
    rmse_thresholds = Dict(
        "pr" => 3.0,      # mm/day
        "rsut" => 29.0,   # W/m²
        "rsutcs" => 10.8,  # W/m²
    )
    return rmse_thresholds
end

if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) != 1
        error("Usage: julia test_rmses.jl <Filepath to simulation data>")
    end
    leaderboard_base_path = ARGS[begin]
    spinup = 3
    test_rmse_thresholds(leaderboard_base_path, spinup)
end
