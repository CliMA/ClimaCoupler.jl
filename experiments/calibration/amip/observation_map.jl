using Dates: Week, Month, Year, Day, Millisecond
import JLD2
import ClimaAnalysis
import ClimaCoupler
import ClimaCalibrate
import ClimaCalibrate: EnsembleBuilder
import ClimaCalibrate.Checker: SequentialIndicesChecker
import CairoMakie
import GeoMakie

# Override JLD2's default_iotype to use IOStream instead of MmapIO
# This avoids Bus errors from memory-mapped files on Lustre filesystem
JLD2.default_iotype() = IOStream

include(
    joinpath(
        pkgdir(ClimaCoupler),
        "experiments",
        "calibration",
        "amip",
        "preprocessing.jl",
    ),
)

"""
    preprocess_sim_vars(vars)

Preprocess sim variables before passing `vars` to `GEnsembleBuilder`.

This should be nearly identical to preprocessing the observational data.

Note that it is a little inefficient to keep these all in memory, but it
simplifies reasoning about the code.
"""
function preprocess_sim_vars(vars)
    # Data loader uses unitless as the units so we do the same here
    vars = set_unitless_units!.(vars)

    vars = select_pressure_levels.(vars, Ref(PRESSURE_LEVELS))
    vars = select_altitude_levels.(vars, Ref(ALTITUDE_LEVELS))
    # Model cloud fraction `cl` is in %, but the CALIPSO observation is a fraction
    # in [0, 1]. Convert sim `cl` to a fraction (and mark it unitless) to match.
    # Model precipitation `pr` is in kg m^-2 s^-1, but the GPCP observation is in
    # mm/day (see GPCPDataLoader). 1 kg m^-2 s^-1 of water = 86400 mm/day, so
    # multiply the sim by 86400 to match. Both obs and sim carry CliMA's
    # downward-negative sign convention (GPCP is flip_sign'd on load), so no sign
    # change is needed here.
    vars = map(vars) do var
        if ClimaAnalysis.short_name(var) == "cl"
            var = ClimaAnalysis.remake(var; data = var.data ./ 100)
            var = ClimaAnalysis.set_units(var, "unitless")
        elseif ClimaAnalysis.short_name(var) == "pr"
            var = ClimaAnalysis.remake(var; data = var.data .* 86400)
            var = ClimaAnalysis.set_units(var, "mm/day")
        end
        return var
    end
    # We do not resample since the simulation variables are already on the
    # simulation grid
    lat_left = -90
    lat_right = 90
    vars = apply_lat_window.(vars, lat_left, lat_right)

    # Restrict the simulation to the points the observation actually covers, BEFORE
    # the zonal mean. Satellite retrievals are not global (MAC `lwp` is ocean-only,
    # NaN over ~54% of grid points), and since zonal_average ignores NaNs the
    # observed zonal mean is an ocean-only average. Without this mask the simulated
    # zonal mean would average all longitudes instead — a different spatial sample.
    # For lwp that is not a small effect: it flips the sign of the area-weighted
    # bias (11.5% low vs 9.1% high), because modelled LWP over land is far lower
    # than over ocean. See apply_coverage_mask.
    coverage_fp = coverage_mask_path(CALIBRATE_CONFIG.output_dir)
    if isfile(coverage_fp)
        coverage_masks = JLD2.load_object(coverage_fp)
        vars = map(vars) do var
            entry = get(coverage_masks, ClimaAnalysis.short_name(var), nothing)
            isnothing(entry) && return var
            return apply_coverage_mask(var, entry...)
        end
    else
        @warn "No coverage masks found; simulation zonal means will average all \
               longitudes even where the observation has no data. Regenerate the \
               observations to create them." coverage_fp
    end

    # Zonal (longitude) mean — must match the obs-side preprocessing in
    # generate_observations.jl so the flattened G aligns with the observation
    # vector (see zonal_average).
    vars = zonal_average.(vars)

    if isfile(NORMALIZATION_STATS_FP)
        # Note: This should not be used with SVDplusDCovariance matrix
        normalization_stats = JLD2.load_object(NORMALIZATION_STATS_FP)
        apply_normalization!.(Ref(normalization_stats), vars)
    end

    # Note: We also do not process the time dimension either since we can rely
    # on GEnsembleBuilder to pick out the right times for us
    return vars
end

"""
    load_and_preprocess_vars(simdir, short_names)

Load and preprocess variables from `simdir`.

The short names `swcre` and `lwcre` are also available.
"""
function load_and_preprocess_vars(simdir::ClimaAnalysis.SimDir, short_names)
    vars = []
    for short_name in short_names
        if short_name == "swcre"
            rsut = get(simdir; short_name = "rsut", reduction = "average", period = "1M")
            rsutcs =
                get(simdir; short_name = "rsutcs", reduction = "average", period = "1M")
            var = rsutcs - rsut
            ClimaAnalysis.set_short_name!(var, "swcre")
            var = ClimaAnalysis.set_units(var, "W m^-2")
            push!(vars, var)
            continue
        elseif short_name == "lwcre"
            rlut = get(simdir; short_name = "rlut", reduction = "average", period = "1M")
            rlutcs =
                get(simdir; short_name = "rlutcs", reduction = "average", period = "1M")
            var = rlutcs - rlut
            ClimaAnalysis.set_short_name!(var, "lwcre")
            var = ClimaAnalysis.set_units(var, "W m^-2")
            push!(vars, var)
            continue
        end
        coord_types = ClimaAnalysis.available_coord_types(
            simdir;
            short_name = short_name,
            reduction = "average",
            period = "1M",
        )

        for coord_type in coord_types
            # Instead of searching though the observation series to determine if
            # the vertical coordinate is pressure or z, we process both of them
            # and let GEnsembleBuilder handle it
            var = get(simdir; short_name, reduction = "average", period = "1M", coord_type)
            push!(vars, var)
        end
    end

    vars = preprocess_sim_vars(vars)
    return vars
end

"""
    process_member_data!(g_ens_builder, diagnostics_folder_path, col_idx)

Process the `col_idx`th member of the G ensemble matrix using `g_ens_builder`
by loading the diagnostics at `diagnostics_folder_path` and processing the
variables.
"""
function process_member_data!(g_ens_builder, diagnostics_folder_path, col_idx)
    short_names = EnsembleBuilder.missing_short_names(g_ens_builder, col_idx)
    @info "Short names: $short_names"

    simdir = ClimaAnalysis.SimDir(diagnostics_folder_path)
    vars = load_and_preprocess_vars(simdir, short_names)

    for variable in vars
        EnsembleBuilder.fill_g_ens_col!(
            g_ens_builder,
            col_idx,
            variable;
            checkers = (SequentialIndicesChecker(),),
            verbose = true,
        )
    end

    return vars
end

"""
    get_job_id(config::ClimaCoupler.CalibrationTools.CalibrateConfig)

Get the job ID from the calibration config.
"""
function get_job_id(config::ClimaCoupler.CalibrationTools.CalibrateConfig)
    (; config_file) = config
    return replace(basename(config_file), ".yml" => "")
end

"""
    ClimaCalibrate.observation_map(interface::CouplerModelInterface, iteration)

Compute the G ensemble matrix from the forward model outputs.
"""
function ClimaCalibrate.observation_map(interface::CouplerModelInterface, iteration)
    (; config) = interface
    (; output_dir) = config
    ekp = JLD2.load_object(ClimaCalibrate.ekp_path(output_dir, iteration))

    g_ens_builder = EnsembleBuilder.GEnsembleBuilder(ekp)
    job_id = get_job_id(config)

    for m in 1:EKP.get_N_ens(ekp)
        member_path = ClimaCalibrate.path_to_ensemble_member(output_dir, iteration, m)
        simdir_path = joinpath(member_path, job_id, "output_active")
        @info "Processing member $m: $simdir_path"
        try
            process_member_data!(g_ens_builder, simdir_path, m)
        catch e
            @error "Ensemble member $m failed" exception = (e, catch_backtrace())
            # Fill failed member column with NaN so EKP can handle the failure
            EnsembleBuilder.fill_g_ens_col!(g_ens_builder, m, NaN)
        end
    end

    g_ens = EnsembleBuilder.get_g_ensemble(g_ens_builder)
    # Too many NaNs - abort (90% threshold like subseasonal)
    if count(isnan, g_ens) > 0.9 * length(g_ens)
        error("Too many NaNs")
    end
    return EnsembleBuilder.is_complete(g_ens_builder) ? g_ens :
           error("G ensemble matrix is not completed")
end
