# Verify the warmup configs still describe the production configuration.
#
#     julia --project=experiments/warmup experiments/warmup/check_config_drift.jl
#
# A warmup image is only useful if the workload builds the same types a production run
# builds. Several model selections are type parameters -- the gravity-wave model is a type
# parameter of `AtmosModel`, for instance -- so a config that drifts from production yields
# an image the run cannot use. That failure is silent: the run simply recompiles, which
# looks exactly like "precompilation didn't help".
#
# This happened. Production switched `orographic_gravity_wave` from `raw_topo` to
# `raw_topo_online` on 2026-10-02 and the warmup config kept the old value for a month.
#
# `CouplerWarmupCommon.check_matches_coupled` guards the other half of the problem -- that
# the component workloads build what the *coupler* builds. Nothing guarded drift against
# the production config itself, which is what this script does.
#
# Exits non-zero on any difference that is not an explicitly recorded deviation, including
# keys added upstream that the warmup config has not picked up.

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

"Differences between a production config and its warmup counterpart, excluding `expected`."
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
