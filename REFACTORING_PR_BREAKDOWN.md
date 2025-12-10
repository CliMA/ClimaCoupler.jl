# Refactoring PR Breakdown

This document outlines how to break down the large refactoring into smaller, manageable PRs.

## Suggested PR Breakdown:

### PR 1: Move `restore.jl` to Checkpointer
- Move `restore!` functions to `Checkpointer.jl`
- Update includes in `climaland_bucket.jl` and `climaatmos.jl`
- Update calls to `Checkpointer.restore!`
- **Why first**: Self-contained, minimal dependencies

### PR 2: Create Input module and move `cli_options.jl`
- Create `Input.jl` module
- Move `argparse_settings()` and `parse_commandline()` from `cli_options.jl`
- Update references to use `Input.argparse_settings` and `Input.parse_commandline`
- **Note**: Can keep `cli_options.jl` temporarily for backward compatibility (or delete if nothing else uses it)

### PR 3: Move `arg_parsing.jl` functions to Input
- Move `get_coupler_config_dict()` and `get_coupler_args()` to `Input.jl`
- Update `setup_run.jl` to use `Input.get_coupler_config_dict` and `Input.get_coupler_args`
- Update other files that use these functions
- Delete `arg_parsing.jl` after migration

### PR 4: Create Postprocessor module structure
- Create `Postprocessor.jl` with basic structure
- Move `postprocess()` and `postprocess_sim()` from `setup_run.jl`
- Move `simulated_years_per_day`, `walltime_per_coupling_step`, `save_sypd_walltime_to_disk`
- Update `SimCoordinator` to call `Postprocessor.simulated_years_per_day`, etc.
- **Why now**: Establishes the module without the large files

### PR 5: Move diagnostics functions to Postprocessor
- Move `coupler_diagnostics.jl` functions to `Postprocessor.jl`
- Update `SimCoordinator` to use `Postprocessor.CD.orchestrate_diagnostics`
- **Why now**: Self-contained, relatively small

### PR 6: Move plotting functions to Postprocessor
- Create `Postprocessor/` directory
- Move `diagnostics_plots.jl` and `debug_plots.jl` to `Postprocessor/`
- Update includes in `Postprocessor.jl`
- Update `postprocess_sim` to call these functions

### PR 7: Move leaderboard functions to Postprocessor
- Move `leaderboard/` directory to `Postprocessor/leaderboard/`
- Update includes in `Postprocessor.jl`
- Update `postprocess_sim` to call leaderboard functions

### PR 8: Move benchmarks to Postprocessor
- Move `benchmarks.jl` to `Postprocessor/`
- Update includes in `Postprocessor.jl`

### PR 9: Move Postprocessor.jl into Postprocessor/ folder
- Move `Postprocessor.jl` to `Postprocessor/Postprocessor.jl`
- Update include path in `ClimaCoupler.jl`
- Cleanup

### PR 10: Cleanup and final touches
- Delete old files (`cli_options.jl`, etc.)
- Update any remaining references
- Move `get_field` methods to appropriate component files
- Final cleanup

## Tips for Managing Interdependencies:

1. **Use feature flags or temporary compatibility shims** if needed
2. **Keep old files temporarily** and add deprecation warnings
3. **Test each PR independently** - each should be functional on its own
4. **Update imports gradually** - you can keep both old and new imports working temporarily

The key is to move in dependency order: start with leaf modules (Checkpointer), then modules that depend on them (Input), then modules that depend on those (Postprocessor).
