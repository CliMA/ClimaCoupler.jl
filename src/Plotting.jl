"""
    Plotting

Module providing plotting functionality for ClimaCoupler simulations.

The actual implementations are provided by the `ClimaCouplerMakieExt` extension
when Makie.jl and related packages are loaded.
"""
module Plotting

import ClimaComms
import ClimaAnalysis as CAN

import ..Interfacer
import ..SimOutput

export postprocess

function make_diagnostics_plots end

function make_ocean_diagnostics_plots end

function debug end

function debug_plot_fields end

function debug_plot! end

function print_extrema end

function plot_global_conservation end

function compute_leaderboard end

function compute_pfull_leaderboard end

# Maps required packages (as a tuple) to the functions provided by that extension
extension_fns = [
    (:Makie, :CairoMakie, :ClimaCoreMakie, :GeoMakie, :Poppler_jll, :Printf) => [
        :make_diagnostics_plots,
        :make_ocean_diagnostics_plots,
        :debug,
        :debug_plot_fields,
        :debug_plot!,
        :print_extrema,
        :plot_global_conservation,
        :compute_leaderboard,
        :compute_pfull_leaderboard,
    ],
    (
        :Makie,
        :CairoMakie,
        :ClimaCoreMakie,
        :GeoMakie,
        :Poppler_jll,
        :Printf,
        :Oceananigans,
    ) => [:debug_plot!, :print_extrema],
]

"""
    is_pkg_loaded(pkg::Symbol)

Check if `pkg` is loaded or not.
"""
function is_pkg_loaded(pkg::Symbol)
    return any(k -> Symbol(k.name) == pkg, keys(Base.loaded_modules))
end

function __init__()
    # Register error hint if a package is not loaded
    if isdefined(Base.Experimental, :register_error_hint)
        Base.Experimental.register_error_hint(MethodError) do io, exc, _argtypes, _kwargs
            for (required_pkgs, fns) in extension_fns
                if Symbol(exc.f) in fns
                    missing_pkgs = [pkg for pkg in required_pkgs if !is_pkg_loaded(pkg)]
                    if !isempty(missing_pkgs)
                        pkg_list = join(string.(missing_pkgs), ", ")
                        print(io, "\nImport $pkg_list to enable `$(exc.f)`.";)
                    end
                end
            end
        end
    end
end

"""
    postprocess(cs; conservation_softfail = false, rmse_check = false)

Process the results after a simulation has completed, including generating
plots, checking conservation, and other diagnostics.
All postprocessing is performed using the root process only, if applicable.

When `conservation_softfail` is true, throw an error if conservation of water
and/or energy is not respected.

When `rmse_check` is true, compute the RMSE against observations and test
that it is below a certain threshold.

The postprocessing includes:
- Energy and water conservation checks (if running SlabPlanet with checks enabled)
- Animations (if not running in MPI)
- AMIP plots of the final state of the model
- Error against observations
- Optional additional atmosphere diagnostics plots
- Plots of useful coupler and component model fields for debugging
"""
function postprocess(
    cs::Interfacer.CoupledSimulation;
    conservation_softfail = false,
    rmse_check = false,
)
    # Only perform postprocessing on root process and if diagnostics handler exists
    if !ClimaComms.iamroot(ClimaComms.context(cs)) || isnothing(cs.diags_handler)
        return nothing
    end

    (;
        coupler_output_dir,
        atmos_output_dir,
        land_output_dir,
        ocean_output_dir,
        artifacts_dir,
    ) = cs.dir_paths

    # Plot generic diagnostics
    @info "Plotting diagnostics for coupler, atmos, land, and ocean"
    make_diagnostics_plots(coupler_output_dir, artifacts_dir, output_prefix = "coupler_")
    make_diagnostics_plots(atmos_output_dir, artifacts_dir, output_prefix = "atmos_")
    make_diagnostics_plots(land_output_dir, artifacts_dir, output_prefix = "land_")

    # Note: slab ocean doesn't have diagnostics, so we only handle Oceananigans here
    make_ocean_diagnostics_plots(ocean_output_dir, artifacts_dir, output_prefix = "ocean_")

    # Plot all model states and coupler fields (useful for debugging)
    ClimaComms.context(cs) isa ClimaComms.SingletonCommsContext && debug(cs, artifacts_dir)

    # If we have enough data (in time, but also enough variables), plot the leaderboard.
    # We need pressure to compute the leaderboard.
    simdir = CAN.SimDir(atmos_output_dir)
    if !isempty(simdir)
        pressure_in_output = "pfull" in CAN.available_vars(simdir)
        t_end = float(cs.t[])
        month = 1
        if t_end > 84600 * 31 * month # 3 months for spin up
            leaderboard_base_path = artifacts_dir
            compute_leaderboard(leaderboard_base_path, atmos_output_dir, month)
            rmse_check && SimOutput.test_rmse_thresholds(atmos_output_dir, month)
            pressure_in_output &&
                compute_pfull_leaderboard(leaderboard_base_path, atmos_output_dir, month)
        end
    end

    # Perform conservation checks if they exist
    if !isnothing(cs.conservation_checks)
        @info "Conservation Check Plots"
        plot_global_conservation(
            cs.conservation_checks.energy,
            cs,
            conservation_softfail,
            figname1 = joinpath(artifacts_dir, "total_energy_bucket.png"),
            figname2 = joinpath(artifacts_dir, "total_energy_log_bucket.png"),
        )
        plot_global_conservation(
            cs.conservation_checks.water,
            cs,
            conservation_softfail,
            figname1 = joinpath(artifacts_dir, "total_water_bucket.png"),
            figname2 = joinpath(artifacts_dir, "total_water_log_bucket.png"),
        )
    end

    # Close all diagnostics file writers
    !isnothing(cs.diags_handler) &&
        map(diag -> close(diag.output_writer), cs.diags_handler.scheduled_diagnostics)

    return nothing
end

end
