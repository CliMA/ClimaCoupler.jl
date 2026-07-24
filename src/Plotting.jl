"""
    Plotting

Module providing plotting functionality for ClimaCoupler simulations.

The actual implementations are provided by the `ClimaCouplerMakieExt` extension
when Makie.jl and related packages are loaded.
"""
module Plotting

import ClimaComms
import ClimaAnalysis as CAN
import Dates

import ..Interfacer
import ..SimOutput

export postprocess

function make_diagnostics_plots end

"""
    make_ocean_diagnostics_plots(output_path, plot_path; output_prefix = "", plot_diagnostics = :all)

Plot the Oceananigans ocean surface diagnostics found under `output_path` to
`plot_path`, one global-map figure per time step (surface speed, SST, SSS, SSH).
`plot_diagnostics` selects `:all` time steps (default) or only the `:last`. The
implementation is provided by the `ClimaCouplerCMIPMakieExt` extension (it needs
Oceananigans loaded).
"""
function make_ocean_diagnostics_plots end

"""
    make_seaice_diagnostics_plots(output_path, plot_path; output_prefix = "", plot_diagnostics = :all)

Plot the Oceananigans sea ice surface diagnostics found under `output_path` to
`plot_path`, one global-map figure per time step (drift speed, concentration,
thickness, top temperature). `plot_diagnostics` selects `:all` time steps
(default) or only the `:last`. The implementation is provided by the
`ClimaCouplerCMIPMakieExt` extension (it needs Oceananigans loaded).
"""
function make_seaice_diagnostics_plots end

"""
    plot_snapshots(cs::Interfacer.CoupledSimulation, date)

Create instantaneous ("snapshot") plots of each component model's key surface
fields at the given `date` and save them, one figure per component, to
`artifacts/<component>/snapshot_<date>.png`.

Snapshots read directly from each component's live state (so they are produced
even for runs that later crash) and are meant to be generated periodically
during a run by a callback (see the `plot_interval` configuration flag). The
implementation is provided by the `ClimaCouplerMakieExt` extension.
"""
function plot_snapshots end

"""
    snapshot_plot(place, var; p_loc = (1, 1))

Plot a single 2-D lat/lon `var` on a global map at grid position `p_loc` using
the shared plot styling (per-variable colormap, Robinson projection, and
coastlines). This is the common rendering path shared by the snapshot callback
and the diagnostics postprocessing. The implementation is provided by the
`ClimaCouplerMakieExt` extension.
"""
function snapshot_plot end

"""
    snapshot_plot_fields(sim)

Return the (small, cheap, curated) set of field names to include in
instantaneous snapshot plots for the component model `sim`. This can be
specialized per component model to control which fields are plotted; the default
is provided by the `ClimaCouplerMakieExt` extension.
"""
function snapshot_plot_fields end

# ============================================================================
# Shared plot styling
#
# The functions below define plot styling (colormaps, map projection,
# coastlines) that is shared between the instantaneous snapshot callback and
# the end-of-run `postprocess` diagnostics. They live in the base package (with
# no Makie dependency) so that the styling policy is owned in one place and both
# plotting paths render with a consistent look. The Makie extension consumes
# these to build the keyword arguments it passes to the actual plotting calls.
# ============================================================================

"""
    PROJECTION

PROJ string for the map projection used by all global (lat/lon) plots. We use
the Robinson projection (`"+proj=robin"`), a compromise projection that shows
the whole globe with reasonable shape and area distortion.
"""
const PROJECTION = "+proj=robin"

# Map from a lowercased *exact* short name to a colormap. Checked before the
# substring rules so that short, ambiguous names (like the diagnostic "pr" for
# precipitation) are matched exactly rather than as substrings.
const _COLORMAP_EXACT = [
    "pr" => :rain,   # precipitation (diagnostic short name)
    "sss" => :haline,
    "sst" => :thermal,
]

# Map from a lowercased substring of a variable's short name to a colormap. The
# first entry whose key is contained in the short name wins, so more specific
# keys should come before more general ones. All colormaps are cmocean schemes
# (available through ColorSchemes.jl) chosen to suit each field type. Keys are
# kept long enough to avoid spurious substring matches.
const _COLORMAP_SUBSTRINGS = [
    # Precipitation and other water fluxes
    "precip" => :rain,
    "rain" => :rain,
    "snow" => :ice,
    # Salinity
    "salinity" => :haline,
    # Temperature
    "temperature" => :thermal,
    # Speeds (non-negative magnitudes)
    "speed" => :speed,
    "wind" => :speed,
    # Depth-like / positive-definite water storage
    "depth" => :deep,
    "thickness" => :deep,
    # Ice concentration / area fractions in [0, 1]
    "concentration" => :ice,
    "fraction" => :dense,
]

# Exact (lowercased) short names for inherently signed velocity-like components.
# These are matched exactly, not as substrings, so that e.g. the "w" of "snow"
# or the "v" of "soil_water" does not trigger a divergent colormap.
const _DIVERGENT_EXACT = ["u", "v", "w"]

# Substrings (lowercased) whose presence in a short name marks the variable as
# inherently signed (can be positive or negative). Unlike `_DIVERGENT_EXACT`,
# these are long enough to match unambiguously as substrings.
const _DIVERGENT_SUBSTRINGS =
    ["velocity", "flux", "anomaly", "tendency", "stress", "curl", "vort"]

"""
    is_divergent(short_name, data)

Return `true` if a variable should be plotted with a divergent colormap.

A variable is treated as divergent if its (lowercased) `short_name` is one of the
signed velocity components (`u`, `v`, `w`), contains one of the signed-quantity
substrings (fluxes, anomalies, ...), or if its `data` spans zero (has both
positive and negative finite values).
"""
function is_divergent(short_name, data)
    name = lowercase(string(short_name))
    name in _DIVERGENT_EXACT && return true
    any(sub -> occursin(sub, name), _DIVERGENT_SUBSTRINGS) && return true
    finite = filter(isfinite, data)
    isempty(finite) && return false
    return minimum(finite) < 0 < maximum(finite)
end

"""
    colormap_for(short_name, data)

Return a sensible colormap (a `Symbol`) for a variable given its `short_name`
and `data`.

The colormap is chosen as follows:
1. If the short name exactly matches a known field type (see `_COLORMAP_EXACT`),
   or contains one of the known field-type substrings (precipitation,
   temperature, salinity, ...; see `_COLORMAP_SUBSTRINGS`), use the colormap for
   that type.
2. Otherwise, if the variable is divergent (see [`is_divergent`](@ref)), use
   the divergent colormap `:balance`.
3. Otherwise, fall back to `:viridis`.
"""
function colormap_for(short_name, data)
    name = lowercase(string(short_name))
    for (key, cmap) in _COLORMAP_EXACT
        name == key && return cmap
    end
    for (key, cmap) in _COLORMAP_SUBSTRINGS
        occursin(key, name) && return cmap
    end
    is_divergent(short_name, data) && return :balance
    return :viridis
end

"""
    symmetric_colorrange(data)

Return a `(-m, m)` colorrange centered on zero, where `m` is the largest finite
absolute value in `data`, or `nothing` if `data` has no usable finite values.
Used to center divergent colormaps on zero so that the zero level maps to the
neutral middle of the colormap.
"""
function symmetric_colorrange(data)
    finite = filter(isfinite, data)
    isempty(finite) && return nothing
    m = maximum(abs, finite)
    m == 0 && return nothing
    return (-m, m)
end

"""
    COLORBAR_WIDTH

Width (in pixels) of the colorbars on global plots. Kept narrow so the colorbar
does not dominate each panel; the variable name and units live in the panel
title instead of on the colorbar.
"""
const COLORBAR_WIDTH = 12

"""
    geo_plot_kwargs(short_name, data; title = nothing)

Build the `more_kwargs` dictionary passed to `ClimaAnalysis`'s global plotting
functions for a variable with the given `short_name` and `data`.

The returned dictionary selects the [`colormap_for`](@ref) the variable, sets
the [`PROJECTION`](@ref) on the axis, hides the lat/lon gridlines, enables black
coastlines, uses a narrow ([`COLORBAR_WIDTH`](@ref)) unlabeled colorbar (the
variable name and units belong in the panel `title`), and, for divergent fields,
centers the colorrange on zero with [`symmetric_colorrange`](@ref).

If `title` is given, it is used as the panel (axis) title.
"""
function geo_plot_kwargs(short_name, data; title = nothing)
    plot_kwargs = Dict{Symbol, Any}(:colormap => colormap_for(short_name, data))
    if is_divergent(short_name, data)
        crange = symmetric_colorrange(data)
        isnothing(crange) || (plot_kwargs[:colorrange] = crange)
    end
    axis_kwargs = Dict{Symbol, Any}(
        :dest => PROJECTION,
        # Hide the lat/lon graticule and its tick labels entirely.
        :xgridvisible => false,
        :ygridvisible => false,
        :xticklabelsvisible => false,
        :yticklabelsvisible => false,
        :xticksvisible => false,
        :yticksvisible => false,
    )
    isnothing(title) || (axis_kwargs[:title] = title)
    return Dict(
        :plot => plot_kwargs,
        :axis => axis_kwargs,
        :coast => Dict(:color => :black),
        # Narrow, label-free colorbar (label = "" overrides the default that
        # `heatmap2D_on_globe!` would otherwise build from the short name/units).
        :cb => Dict(:label => "", :width => COLORBAR_WIDTH),
    )
end

"""
    component_artifacts_dir(artifacts_dir, sim_name)

Return (creating it if necessary) the per-component subdirectory of
`artifacts_dir` where plots for the component named `sim_name` are stored, e.g.
`artifacts/atmos_sim`. Grouping plots by component keeps the top-level
`artifacts` directory manageable once monthly snapshots are added.
"""
function component_artifacts_dir(artifacts_dir, sim_name)
    dir = joinpath(artifacts_dir, string(sim_name))
    mkpath(dir)
    return dir
end

function plot_global_conservation end

function compute_leaderboard end

function compute_pfull_leaderboard end

# Maps required packages (as a tuple) to the functions provided by that extension
extension_fns = [
    (:Makie, :CairoMakie, :ClimaCoreMakie, :GeoMakie, :Poppler_jll, :Printf) => [
        :make_diagnostics_plots,
        :plot_snapshots,
        :snapshot_plot,
        :snapshot_plot_fields,
        :plot_global_conservation,
        :compute_leaderboard,
        :compute_pfull_leaderboard,
    ],
    # Ocean/sea-ice diagnostics read Oceananigans JLD2 output, so they also need
    # Oceananigans loaded (provided by ClimaCouplerCMIPMakieExt).
    (:Makie, :CairoMakie, :GeoMakie, :Oceananigans, :Poppler_jll, :Printf) =>
        [:make_ocean_diagnostics_plots, :make_seaice_diagnostics_plots],
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
    postprocess(cs; conservation_softfail = false, rmse_check = false, plot_diagnostics = :all)

Process the results after a simulation has completed, including generating
plots, checking conservation, and other diagnostics.
All postprocessing is performed using the root process only, if applicable.

When `conservation_softfail` is true, throw an error if conservation of water
and/or energy is not respected.

When `rmse_check` is true, compute the RMSE against observations and test
that it is below a certain threshold.

`plot_diagnostics` controls how much of the diagnostics time series is plotted:
- `:all` (default): plot every saved time step, producing one summary file per
  time step per component (named with the date).
- `:last`: plot only the last saved time step (i.e. the final averaging window),
  producing a single summary file per component.

The postprocessing includes:
- Diagnostics plots for the coupler and each component model, written to
  per-component subdirectories of `artifacts` (e.g. `artifacts/atmos_sim`).
  Two-dimensional lat/lon fields use the shared plot styling (per-variable
  colormap, Robinson projection, and coastlines).
- The atmosphere leaderboard comparing against observations (global runs only)
- Optionally, an RMSE check against observations (`rmse_check`)
- Energy and water conservation checks (if running SlabPlanet with checks enabled)

Note that instantaneous "snapshot" plots of the component and coupler fields are
produced *during* the run by a callback (see [`plot_snapshots`](@ref) and the
`plot_interval` configuration flag), so that they are available even for runs
that crash before completion. `postprocess` covers the end-of-run diagnostics
that require the full time series.
"""
function postprocess(
    cs::Interfacer.CoupledSimulation;
    conservation_softfail = false,
    rmse_check = false,
    plot_diagnostics = :all,
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
        ice_output_dir,
        artifacts_dir,
    ) = cs.dir_paths

    # Plot generic diagnostics into per-component subdirectories of `artifacts`
    # (e.g. `artifacts/atmos_sim`) so that the top-level artifacts directory
    # stays manageable alongside the periodic snapshot plots.
    @info "Plotting diagnostics for coupler, atmos, land, and ocean"
    make_diagnostics_plots(
        coupler_output_dir,
        component_artifacts_dir(artifacts_dir, "coupler");
        plot_diagnostics,
    )
    make_diagnostics_plots(
        atmos_output_dir,
        component_artifacts_dir(artifacts_dir, "atmos_sim");
        plot_diagnostics,
    )
    make_diagnostics_plots(
        land_output_dir,
        component_artifacts_dir(artifacts_dir, "land_sim");
        plot_diagnostics,
    )

    # Note: slab ocean/prescribed sea ice don't have diagnostics, so these only
    # produce plots for the Oceananigans ocean and ClimaSeaIce sea ice models.
    make_ocean_diagnostics_plots(
        ocean_output_dir,
        component_artifacts_dir(artifacts_dir, "ocean_sim");
        plot_diagnostics,
    )
    make_seaice_diagnostics_plots(
        ice_output_dir,
        component_artifacts_dir(artifacts_dir, "ice_sim");
        plot_diagnostics,
    )

    # Helper function to find a tuple of (short_name, reduction, period, coord_type)
    # whose period is "1M".
    function find_first_1M_tuple(simdir)
        for short_name in CAN.available_vars(simdir)
            for reduction in CAN.available_reductions(simdir; short_name)
                for period in CAN.available_periods(simdir; short_name, reduction)
                    period == "1M" || continue
                    for coord_type in
                        CAN.available_coord_types(simdir; short_name, reduction, period)
                        return (short_name, reduction, period, coord_type)
                    end
                end
            end
        end
        return nothing
    end

    # Leaderboard requires global spatial data; skip in column / SCM mode.
    if !Interfacer.is_column_mode(cs)
        simdir = CAN.SimDir(atmos_output_dir)
        first_1M_tuple = find_first_1M_tuple(simdir)
        if !isempty(simdir) && !isnothing(first_1M_tuple)
            # We need pressure to compute the leaderboard.
            pressure_in_output = "pfull" in CAN.available_vars(simdir)

            # Make sure we have enough data to compute the leaderboard.
            short_name, reduction, period, coord_type = first_1M_tuple
            var_1M = get(simdir; short_name, reduction, period, coord_type)
            start_date = first(CAN.dates(var_1M))
            end_date = last(CAN.dates(var_1M))
            spinup = 3
            if end_date >= start_date + Dates.Month(spinup)
                leaderboard_base_path = artifacts_dir
                compute_leaderboard(leaderboard_base_path, atmos_output_dir, spinup)
                rmse_check && SimOutput.test_rmse_thresholds(atmos_output_dir, spinup)
                pressure_in_output && compute_pfull_leaderboard(
                    leaderboard_base_path,
                    atmos_output_dir,
                    spinup,
                )
            end
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
