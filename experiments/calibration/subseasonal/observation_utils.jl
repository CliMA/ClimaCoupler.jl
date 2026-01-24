import ClimaCoupler
using Statistics
import Dates
using ClimaAnalysis
import ClimaCalibrate

# Units for each short name used in calibration (target units for comparison)
const var_units = Dict(
    "pr" => "mm d^-1",  # Precipitation in mm/day for comparison
    "mslp" => "Pa",
    "tas" => "K",
    "tas - ta" => "K",
    "ts" => "K",
    "sp" => "Pa",
    "hfls" => "W m^-2",
    "hfss" => "W m^-2",
    "rsus" => "W m^-2",
    "rlus" => "W m^-2",
)

# =============================================================================
# UNIT CONVERSION FACTORS FOR PRECIPITATION
# Both model and ERA5 data need to be converted to the same units (mm/day)
# =============================================================================

# Model precipitation conversion:
# - Model stores pr as NEGATIVE kg/m²/s (downward flux convention)
# - Convert to positive mm/day: multiply by -86400
# - 1 kg/m² = 1 mm of water, 86400 s/day, negative to flip sign
const MODEL_PR_CONVERSION = -86400.0

# ERA5 precipitation conversion:
# - ERA5 'tp' (total precipitation) is in meters per hour (for daily-mean files)
# - Convert to mm/day: multiply by 24000
# - 1000 mm/m * 24 hours/day = 24000
const ERA5_PR_CONVERSION = 24000.0

"""
    apply_pr_conversion(data, source::Symbol)

Apply precipitation unit conversion based on the data source.
- source = :model → multiply by MODEL_PR_CONVERSION
- source = :era5 → multiply by ERA5_PR_CONVERSION
"""
function apply_pr_conversion(data, source::Symbol)
    if source == :model
        return data .* MODEL_PR_CONVERSION
    elseif source == :era5
        return data .* ERA5_PR_CONVERSION
    else
        error("Unknown source: $source. Use :model or :era5")
    end
end

function remove_global_mean(var)
    mean_var = ClimaAnalysis.Var.average_lonlat(var; weighted = true)
    mean_data = mean_var.data
    return ClimaAnalysis.replace(val -> val - mean_data[1], var)
end

"""
    compute_nh_average(var, land_mask, use_land_mask)

Compute area-weighted average over the Northern Hemisphere (lat > 0).

Arguments:
- var: OutputVar with (lon, lat) data
- land_mask: 2D boolean array (lon, lat), true = land
- use_land_mask: if true, only average over land points in NH

Returns: scalar (the weighted average)
"""
function compute_nh_average(var, land_mask, use_land_mask::Bool)
    lats = ClimaAnalysis.latitudes(var)
    lons = ClimaAnalysis.longitudes(var)
    data = var.data  # (lon, lat)
    
    # Find NH indices (lat > 0)
    nh_lat_indices = findall(lat -> lat > 0, lats)
    
    # Cosine latitude weights for NH
    nh_weights = cos.(deg2rad.(lats[nh_lat_indices]))
    
    # Collect weighted data
    weighted_sum = 0.0
    weight_total = 0.0
    
    for (j_local, j_global) in enumerate(nh_lat_indices)
        w = nh_weights[j_local]
        for i in eachindex(lons)
            val = data[i, j_global]
            # Skip NaN values
            if isnan(val)
                continue
            end
            # Apply land mask if requested
            if use_land_mask && !land_mask[i, j_global]
                continue
            end
            weighted_sum += val * w
            weight_total += w
        end
    end
    
    if weight_total == 0.0
        return NaN
    end
    
    return weighted_sum / weight_total
end
