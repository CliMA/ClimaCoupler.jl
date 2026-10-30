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
    preprocess_sim_vars(vars, config)

Preprocess sim variables before passing `vars` to `GEnsembleBuilder`.

This should be nearly identical to preprocessing the observational data.

Note that it is a little inefficient to keep these all in memory, but it
simplifies reasoning about the code.
"""
function preprocess_sim_vars(vars, config)
    # Data loader uses unitless as the units so we do the same here
    vars = set_unitless_units!.(vars)

    # We do not resample since the simulation variables are already on the
    # simulation grid
    vars = select_pressure_levels.(vars, Ref(config.pressure_levels))

    # Where each observation has data, saved by `generate_observations.jl`.
    coverage_fp = coverage_mask_path(config.output_dir)
    coverage_masks =
        isfile(coverage_fp) ? JLD2.load_object(coverage_fp) : Dict{String, Any}()
    isfile(coverage_fp) ||
        @warn "No coverage masks found, so the simulation will be averaged \
over points the observation does not cover. Regenerate the observations to create \
them." coverage_fp

    vars = loss_values(vars, coverage_masks, config)

    # Note: We also do not process the time dimension either since we can rely
    # on GEnsembleBuilder to pick out the right times for us
    return vars
end

"""
    load_and_preprocess_vars(simdir, short_names, config)

Load the model output in `simdir` for the loss fields `short_names`, from
`CalibrationTools.model_var`, and preprocess it.
"""
function load_and_preprocess_vars(simdir::ClimaAnalysis.SimDir, short_names, config)
    vars = []
    for short_name in short_names
        # A field on both pressure and altitude gives both; GEnsembleBuilder picks the one
        # the observation uses.
        append!(vars, CalibrationTools.model_var(simdir, short_name))
    end

    vars = preprocess_sim_vars(vars, config)
    return vars
end

"""
    process_member_data!(g_ens_builder, diagnostics_folder_path, col_idx, config)

Process the `col_idx`th member of the G ensemble matrix using `g_ens_builder`
by loading the diagnostics at `diagnostics_folder_path` and processing the
variables.
"""
function process_member_data!(g_ens_builder, diagnostics_folder_path, col_idx, config)
    short_names = EnsembleBuilder.missing_short_names(g_ens_builder, col_idx)
    @info "Short names: $short_names"

    simdir = ClimaAnalysis.SimDir(diagnostics_folder_path)
    vars = load_and_preprocess_vars(simdir, short_names, config)

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
            process_member_data!(g_ens_builder, simdir_path, m, config)
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
