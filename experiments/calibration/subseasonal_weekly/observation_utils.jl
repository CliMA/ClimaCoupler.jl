# Import shared definitions from subseasonal pipeline:
#   var_units, remove_global_mean, get_var, preprocess_var
include(
    joinpath(
        @__DIR__,
        "..",
        "subseasonal",
        "observation_utils.jl",
    ),
)

# =============================================================================
# Weekly-specific overrides and additions
# =============================================================================

# Weekly pipeline uses mm/day for precipitation (ERA5 data is converted)
var_units["pr"] = "mm d^-1"

# ERA5 precipitation conversion:
# - ERA5 'tp' (total precipitation) is in meters per hour (for daily-mean files)
# - Convert to mm/day: multiply by 24000
# - 1000 mm/m * 24 hours/day = 24000
const ERA5_PR_CONVERSION = 24000.0
