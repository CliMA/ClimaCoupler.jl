using Dates: Week, Month, Year, Day, Millisecond
import JLD2
import ClimaAnalysis
import ClimaCoupler
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
include(
    joinpath(
        pkgdir(ClimaCoupler),
        "experiments",
        "calibration",
        "amip",
        "post_analyze_iteration.jl",
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
    # We do not resample since the simulation variables are already on the
    # simulation grid
    lat_left = -60
    lat_right = 60
    vars = apply_lat_window.(vars, lat_left, lat_right)

    # Note: We also do not process the time dimension either since we can rely
    # on GEnsembleBuilder to pick out the right times for us
    return var
end

function process_member_data!(g_ens_builder, diagnostics_folder_path, col_idx, iteration)
    short_names = EnsembleBuilder.missing_short_names(g_ens_builder, col_idx)
    sample_date_ranges = CALIBRATE_CONFIG.sample_date_ranges[iteration + 1]
    @info "Short names: $short_names"

    simdir = ClimaAnalysis.SimDir(diagnostics_folder_path)
    for short_name in short_names

        coord_types = available_coord_types(
            simdir;
            short_name = short_name,
            reduction = "average",
            period = "1m",
        )

        vars = []
        for coord_type in coord_types
            # Instead of searching though the observation series to determine if
            # the vertical coordinate is pressure or z, we process both of them
            # and let GEnsembleBuilder handle it
            var = get(simdir; short_name, reduction = "average", period = "1m", coord_type)
            push!(vars, var)
        end

        vars = preprocess_sim_vars(vars)

        for var in vars
            EnsembleBuilder.fill_g_ens_col!(
                g_ens_builder,
                col_idx,
                var;
                checkers = (SequentialIndicesChecker(),),
                verbose = true,
            )
        end
    end

    return nothing
end

# Override observation_map to use correct job_id path
function ClimaCalibrate.observation_map(iteration)
    output_dir = CALIBRATE_CONFIG.output_dir
    ekp = JLD2.load_object(ClimaCalibrate.ekp_path(output_dir, iteration))

    g_ens_builder = EnsembleBuilder.GEnsembleBuilder(ekp)
    job_id = get_job_id()

    for m in 1:EKP.get_N_ens(ekp)
        member_path = ClimaCalibrate.path_to_ensemble_member(output_dir, iteration, m)
        simdir_path = joinpath(member_path, job_id, "output_active")
        @info "Processing member $m: $simdir_path"
        try
            process_member_data!(g_ens_builder, simdir_path, m, iteration)
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
