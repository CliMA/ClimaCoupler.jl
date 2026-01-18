"""
    Plotting

Module providing plotting functionality for ClimaCoupler simulations.

The actual implementations are provided by the `ClimaCouplerMakieExt` extension
when Makie.jl and related packages are loaded.
"""
module Plotting

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

end
