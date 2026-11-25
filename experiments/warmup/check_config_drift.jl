# Check that the warmup configs still match the production configs.
#
#     julia --project=experiments/warmup experiments/warmup/check_config_drift.jl
#
# Model selections are type parameters, so a drifted config yields images production
# cannot use. Exits non-zero on any difference not in `ATMOS_EXPECTED` or `CMIP_EXPECTED`.

using YAML

const REPO = normpath(joinpath(@__DIR__, "..", ".."))

# Deviations we intend. Each key maps to the reason it differs; anything else is drift.
const ATMOS_EXPECTED = Dict(
    "initial_condition" => "WeatherModel needs an ERA5 artifact that is cluster-only; see README",
    "t_end" => "the workload runs briefly",
    # The coupler applies these itself in `get_atmos_config_dict`, so a standalone
    # atmosphere workload has to replicate them or it caches the wrong types.
    "surface_setup" => "coupler override",
    "toml" => "coupler override",
    "albedo_model" => "coupler override",
    "log_progress" => "coupler override",
    "output_dir_style" => "coupler override",
    "FLOAT_TYPE" => "coupler override",
    "topography" => "coupler override",
    "topo_smoothing" => "coupler override",
)

const CMIP_EXPECTED = Dict(
    "atmos_config_file" => "points at the warmup atmosphere config",
    "t_end" => "the workload runs briefly",
    "rmse_check" => "only gates `postprocess`, which the workload never reaches",
    # Shortened so every schedule fires at least once during the short run; `time_to_period`
    # maps all of these to `Dates.Millisecond`, so the schedule type is preserved.
    "checkpoint_dt" => "shortened so the schedule fires",
    "coupler_diagnostics_period" => "shortened so the schedule fires",
    "dt_cloud_fraction" => "shortened so the schedule fires",
    "dt_rad" => "shortened so the schedule fires",
    "land_diagnostics_period" => "shortened so the schedule fires",
    "ocean_diagnostic_interval" => "shortened so the schedule fires",
    "seaice_diagnostic_interval" => "shortened so the schedule fires",
    "extra_atmos_diagnostics" => "periods shortened so the schedule fires",
)

"""
    drift(production, warmup, expected)

Return sorted `(key, reason)` pairs where the config dicts differ, skipping `expected` keys.
"""
function drift(production, warmup, expected)
    issues = Tuple{String, String}[]
    for key in union(keys(production), keys(warmup))
        haskey(expected, key) && continue
        inprod, inwarm = haskey(production, key), haskey(warmup, key)
        if inprod && !inwarm
            push!(
                issues,
                (
                    key,
                    "missing from the warmup config (production has $(repr(production[key])))",
                ),
            )
        elseif inwarm && !inprod
            push!(issues, (key, "only in the warmup config ($(repr(warmup[key])))"))
        elseif production[key] != warmup[key]
            push!(
                issues,
                (key, "production $(repr(production[key])), warmup $(repr(warmup[key]))"),
            )
        end
    end
    sort!(issues, by = first)
end

const PAIRS = [
    (
        "atmosphere",
        joinpath(REPO, "config", "atmos_configs", "climaatmos_progedmf.yml"),
        joinpath(REPO, "experiments", "warmup", "configs", "atmos_warmup.yml"),
        ATMOS_EXPECTED,
    ),
    (
        "coupled",
        joinpath(REPO, "config", "longrun_configs", "cmip_progedmf_land.yml"),
        joinpath(REPO, "experiments", "warmup", "configs", "cmip_warmup.yml"),
        CMIP_EXPECTED,
    ),
]

"""
    coupler_override_drift(atmos, coupled)

Return `(key, why)` pairs where a value the atmosphere config replicates from the coupler
disagrees with the coupled config.

`drift` deliberately skips these keys: the coupler sets them itself, so the atmosphere
config must differ from the *production atmosphere* config. That makes them invisible to
the comparison above, which is how `FLOAT_TYPE` and `topography` once defaulted to
`Float32` and `Val{:NoWarp}` in a standalone atmosphere workload. They are not free to be
anything, though -- they have to agree with the coupled config, which is itself checked
against production. So cross-check them here.
"""
function coupler_override_drift(atmos, coupled)
    issues = Tuple{String, String}[]
    for (key, reason) in ATMOS_EXPECTED
        reason == "coupler override" || continue
        haskey(coupled, key) || continue   # not every override is a coupled-config key
        haskey(atmos, key) || begin
            push!(
                issues,
                (
                    key,
                    "coupled config sets $(repr(coupled[key])), atmosphere config omits it",
                ),
            )
            continue
        end
        atmos[key] == coupled[key] || push!(
            issues,
            (key, "atmosphere $(repr(atmos[key])), coupled $(repr(coupled[key]))"),
        )
    end
    sort!(issues, by = first)
end

total = sum(PAIRS) do (label, prodpath, warmpath, expected)
    issues = drift(YAML.load_file(prodpath), YAML.load_file(warmpath), expected)
    println(
        label,
        ": ",
        isempty(issues) ? "matches production" : "$(length(issues)) drifted",
    )
    for (key, why) in issues
        println("  ", rpad(key, 32), why)
    end
    length(issues)
end

# The cross-check above runs on the warmup configs themselves, not against production.
let atmos = YAML.load_file(PAIRS[1][3]), coupled = YAML.load_file(PAIRS[2][3])
    issues = coupler_override_drift(atmos, coupled)
    println(
        "coupler overrides: ",
        isempty(issues) ? "atmosphere config agrees with the coupled config" :
        "$(length(issues)) disagree",
    )
    for (key, why) in issues
        println("  ", rpad(key, 32), why)
    end
    global total += length(issues)
end

println()
if total == 0
    println("Warmup configs match production apart from the recorded deviations.")
else
    println(
        "$total drifted setting(s). Either update the warmup config to match production,",
    )
    println(
        "or, if the difference is intentional, record it in this script with its reason.",
    )
    println(
        "Settings that are type parameters (gravity-wave model, microphysics, turbconv,",
    )
    println("radiation mode) silently invalidate the whole cache when they drift.")
    exit(1)
end
