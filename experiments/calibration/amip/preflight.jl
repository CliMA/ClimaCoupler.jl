# Pre-flight validator for a calibration launch. See PREFLIGHT_PLAN.md.
#
# Runs four check groups and prints one line per check (PASS, WARN, or FAIL).
# Exits nonzero if any check fails. The launch procedure is: generate
# observations, run preflight, submit only on success.
#
#   A. Config consistency: priors, ensemble size, walltime estimate, dates,
#      grid divisibility.
#   B. Observation vector: sample count, block names and lengths, noise
#      levels, coverage masks.
#   C. Run directory hygiene: partial members, CALIBRATION_CONFIG env.
#   D. Parameter wiring: each prior must change the ClimaAtmos parameter
#      struct when perturbed. This is the check that catches parameters the
#      model silently ignores.
#
# Usage, from the repository root:
#
#   julia --project=experiments/AMIP experiments/calibration/amip/preflight.jl \
#       experiments/calibration/amip/config/pipeline_test.jl [--skip-wiring] [--sypd 2.0]
#
# A config may define PREFLIGHT_NON_ATMOS_PARAMS, a vector of parameter names
# that do not enter through the ClimaAtmos parameter struct (for example land
# parameters). Those parameters get WARN (unchecked) instead of FAIL in check D.

ENV["CLIMACOMMS_CONTEXT"] = "SINGLETON"

import Dates
import JLD2
import Statistics
import TOML
import ClimaAnalysis
import ClimaCalibrate
import ClimaCoupler
import ClimaCoupler: CalibrationTools
import ClimaParams as CP
import EnsembleKalmanProcesses as EKP
import EnsembleKalmanProcesses.ParameterDistributions as PD

# ---------------------------------------------------------------------------
# Arguments and config
# ---------------------------------------------------------------------------

function parse_preflight_args(args)
    config_path = ""
    skip_wiring = false
    sypd = 2.0
    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--skip-wiring"
            skip_wiring = true
        elseif a == "--sypd"
            i += 1
            sypd = parse(Float64, args[i])
        elseif startswith(a, "--")
            error("Unknown flag $a")
        else
            config_path = abspath(a)
        end
        i += 1
    end
    isempty(config_path) &&
        error("Usage: preflight.jl <config.jl> [--skip-wiring] [--sypd X]")
    isfile(config_path) || error("Config file $config_path does not exist")
    return (; config_path, skip_wiring, sypd)
end

# ---------------------------------------------------------------------------
# Result reporting
# ---------------------------------------------------------------------------

const PREFLIGHT_RESULTS = NamedTuple[]

function record!(status, name, msg = "")
    push!(PREFLIGHT_RESULTS, (; status, name, msg))
    line = rpad(status, 4) * " | " * name
    isempty(msg) || (line *= " | " * msg)
    println(line)
    return nothing
end

pass!(name, msg = "") = record!("PASS", name, msg)
warn!(name, msg = "") = record!("WARN", name, msg)
fail!(name, msg = "") = record!("FAIL", name, msg)

# ---------------------------------------------------------------------------
# A. Config consistency
# ---------------------------------------------------------------------------

"Worker walltime in run_calibration.jl and the PBS job walltime, in hours."
const JOB_WALLTIME_HOURS = 12.0

function check_config(sypd)
    names = PD.get_name(Main.PRIORS)
    if isempty(names)
        fail!("priors", "no parameters in PRIORS")
    elseif length(unique(names)) != length(names)
        fail!("priors", "duplicate parameter names")
    else
        pass!("priors", "$(length(names)) parameters, unique names")
    end

    n_par = ndims(Main.PRIORS)
    n_ens = 2 * n_par + 1
    pass!("ensemble", "TransformUnscented needs $n_ens members ($n_ens workers)")

    cfg = Main.CALIBRATE_CONFIG
    sim_days = maximum(cfg.sample_date_ranges) do (s, e)
        Dates.value((e + cfg.extend) - (s - cfg.spinup)) / 86_400_000
    end
    sim_days = round(sim_days; digits = 1)
    hours_per_iter = sim_days / 365.25 / sypd * 24
    total_hours = cfg.n_iterations * hours_per_iter + 1.0
    msg =
        "$(sim_days) sim days per member, " *
        "$(round(hours_per_iter; digits = 1)) h per iteration, " *
        "$(round(total_hours; digits = 1)) h total at $(sypd) SYPD"
    if total_hours > JOB_WALLTIME_HOURS
        warn!("walltime", msg * " exceeds the $(JOB_WALLTIME_HOURS) h job")
    else
        pass!("walltime", msg)
    end

    if isdefined(Main, :COVARIANCE_DATE_RANGES)
        covranges = Main.COVARIANCE_DATE_RANGES
        if length(covranges) < 2
            fail!("dates", "need at least 2 covariance date ranges")
        else
            outside = filter(r -> !(r in covranges), cfg.sample_date_ranges)
            if isempty(outside)
                pass!(
                    "dates",
                    "$(length(covranges)) covariance samples contain all targets",
                )
            else
                fail!(
                    "dates",
                    "sample ranges $(outside) are not covariance dates " *
                    "(SVDplusD requires this)",
                )
            end
        end
    else
        warn!("dates", "config defines no COVARIANCE_DATE_RANGES")
    end

    nlon, nlat = comparison_grid_size(cfg.config_file)
    if isdefined(Main, :COARSEN_FACTOR)
        f = Main.COARSEN_FACTOR
        if nlon % f == 0 && nlat % f == 0
            pass!(
                "grid",
                "$(nlon)x$(nlat) divides by $f -> $(div(nlon, f))x$(div(nlat, f))",
            )
        else
            fail!("grid", "$(nlon)x$(nlat) does not divide by COARSEN_FACTOR=$f")
        end
    else
        pass!("grid", "zonal mean over $(nlon)x$(nlat)")
    end
    return nothing
end

"Comparison grid size, matching get_lonlat_regridder in preprocessing.jl."
function comparison_grid_size(config_file)
    config_dict = ClimaCoupler.Input.get_coupler_config_dict(config_file)
    pts = get(config_dict, "netcdf_interpolation_num_points", nothing)
    isnothing(pts) || return (pts[1], pts[2])
    h_elem = get(config_dict, "h_elem", 12)
    nlon = h_elem * 4 * 3
    return (nlon, nlon ÷ 2)
end

# ---------------------------------------------------------------------------
# B. Observation vector
# ---------------------------------------------------------------------------

function check_observations()
    cfg = Main.CALIBRATE_CONFIG
    obs_path = joinpath(cfg.output_dir, "observation_vec.jld2")
    if !isfile(obs_path)
        warn!(
            "obs vector",
            "not generated yet. Run generate_observations.jl, then preflight again",
        )
        return nothing
    end
    obs_vec = JLD2.load_object(obs_path)

    n_samples = length(cfg.sample_date_ranges)
    if length(obs_vec) == n_samples
        pass!("obs samples", "$(length(obs_vec)) samples match sample_date_ranges")
    else
        fail!(
            "obs samples",
            "$(length(obs_vec)) samples but $(n_samples) sample_date_ranges. " *
            "Regenerate the observations",
        )
    end

    obs = first(obs_vec)
    y = EKP.get_obs(obs)
    if all(isfinite, y)
        pass!("obs values", "$(length(y)) constraints, all finite")
    else
        fail!("obs values", "$(count(!isfinite, y)) non-finite entries")
    end

    obs_vars = ClimaCalibrate.ObservationRecipe.reconstruct_vars(obs)
    names = ClimaAnalysis.short_name.(obs_vars)
    lens = [length(ClimaAnalysis.flatten(v).data) for v in obs_vars]
    if sort(names) == sort(cfg.short_names)
        pass!("obs blocks", join(("$n=$l" for (n, l) in zip(names, lens)), ", "))
    else
        fail!("obs blocks", "blocks $(names) do not match short_names $(cfg.short_names)")
    end
    if sum(lens) != length(y)
        fail!("obs layout", "block lengths sum to $(sum(lens)) but y has $(length(y))")
    end

    cov = EKP.get_obs_noise_cov(obs)
    diag = [cov[i, i] for i in 1:size(cov, 1)]
    if all(d -> isfinite(d) && d > 0, diag)
        pass!("obs noise", "covariance diagonal finite and positive")
    else
        fail!("obs noise", "covariance diagonal has zero or non-finite entries")
    end

    off = 0
    for (n, l) in zip(names, lens)
        rng = (off + 1):(off + l)
        off += l
        last(rng) <= length(y) || continue
        rms = sqrt(Statistics.mean(abs2, y[rng]))
        rel = Statistics.median(sqrt.(diag[rng])) / rms
        msg = "relative noise $(round(100 * rel; digits = 1))% of signal RMS"
        (0.05 <= rel <= 2.0) ? pass!("noise $n", msg) : warn!("noise $n", msg)
    end

    check_coverage_masks(cfg)
    return nothing
end

function check_coverage_masks(cfg)
    mask_path = Main.coverage_mask_path(cfg.output_dir)
    if !isfile(mask_path)
        warn!("coverage masks", "no coverage_masks.jld2 next to the observations")
        return nothing
    end
    masks = JLD2.load_object(mask_path)
    spatial = ("longitude", "latitude", "altitude")
    for name in cfg.short_names
        if !haskey(masks, name)
            fail!("mask $name", "no coverage mask saved for this variable")
            continue
        end
        dims, mask = masks[name]
        canon = Main.canonical_dim_name.(dims)
        bad = filter(c -> !(c in spatial), canon)
        if !isempty(bad)
            fail!(
                "mask $name",
                "dimensions $(dims) canonicalize to $(canon). " *
                "$(bad) are not spatial, the mask cannot apply to the simulation",
            )
        elseif all(mask)
            fail!(
                "mask $name",
                "mask is saturated (every point missing). " *
                "Check the covariance date ranges against the data record",
            )
        else
            frac = round(count(mask) / length(mask); digits = 3)
            pass!("mask $name", "masks fraction $frac of points")
        end
    end
    return nothing
end

# ---------------------------------------------------------------------------
# C. Run directory hygiene
# ---------------------------------------------------------------------------

function check_run_dir(config_path)
    cfg = Main.CALIBRATE_CONFIG
    partial = CalibrationTools.find_partial_members(cfg.output_dir)
    if isempty(partial)
        pass!("run dir", "no partial members in $(cfg.output_dir)")
    else
        fail!(
            "run dir",
            "$(length(partial)) partial members (first: $(first(partial))). " *
            "Resuming over them writes wrong-dated diagnostics. " *
            "Remove everything in each member directory except parameters.toml",
        )
    end

    env_cfg = get(ENV, "CALIBRATION_CONFIG", nothing)
    if isnothing(env_cfg)
        pass!("env", "CALIBRATION_CONFIG not set (set it in the submit script)")
    elseif abspath(env_cfg) == config_path
        pass!("env", "CALIBRATION_CONFIG matches this config")
    else
        fail!("env", "CALIBRATION_CONFIG points at $(env_cfg), not $(config_path)")
    end
    return nothing
end

# ---------------------------------------------------------------------------
# D. Parameter wiring
# ---------------------------------------------------------------------------

"""
Write a member-style parameter TOML with the given names and values.
"""
function write_param_toml(path, names, values)
    open(path, "w") do io
        for (n, v) in zip(names, values)
            println(io, "[", n, "]")
            println(io, "value = ", Float64(v))
            println(io, "type = \"float\"")
            println(io)
        end
    end
    return path
end

"""
Build the ClimaAtmos parameter struct through the same TOML merge a member
uses: the config's coupler_toml files, then the member parameter file on top.
No spaces are built and no GPU is needed.
"""
function atmos_params_from_toml(CA, base_tomls, member_toml, FT, gw_flags)
    override_file = CP.merge_toml_files(vcat(base_tomls, [member_toml]); override = true)
    toml_dict = CP.create_toml_dict(FT; override_file)
    return CA.ClimaAtmosParameters(
        toml_dict;
        has_non_orographic_gw = gw_flags.non_orographic,
        has_orographic_gw = gw_flags.orographic,
    )
end

function check_wiring(CA, scratch)
    cfg = Main.CALIBRATE_CONFIG
    config_dict = ClimaCoupler.Input.get_coupler_config_dict(cfg.config_file)
    base_tomls = map(config_dict["coupler_toml"]) do f
        isfile(f) ? abspath(f) : joinpath(pkgdir(ClimaCoupler), f)
    end
    for f in base_tomls
        isfile(f) || (fail!("wiring", "coupler_toml $f does not exist"); return nothing)
    end
    FT = get(config_dict, "FLOAT_TYPE", "Float32") == "Float64" ? Float64 : Float32
    gw_flags = (;
        non_orographic = get(config_dict, "non_orographic_gravity_wave", false),
        orographic = !isnothing(get(config_dict, "orographic_gravity_wave", nothing)),
    )

    names = PD.get_name(Main.PRIORS)
    # For a single scalar parameter these statistics come back as scalars.
    _asvec(x) = x isa Number ? [x] : collect(x)
    u0 = _asvec(Statistics.mean(Main.PRIORS))
    su = sqrt.(_asvec(Statistics.var(Main.PRIORS)))
    center = _asvec(PD.transform_unconstrained_to_constrained(Main.PRIORS, u0))

    unchecked =
        isdefined(Main, :PREFLIGHT_NON_ATMOS_PARAMS) ? Main.PREFLIGHT_NON_ATMOS_PARAMS :
        String[]

    # `<base>_E<index>` priors calibrate single elements of a vector
    # parameter. write_param_toml emits them as scalars under their sampled
    # names, which ClimaAtmos has never heard of, so the wiring check would
    # report "changes nothing" and FAIL every element prior. Splice them into
    # the base vector exactly as model_interface.jl does before a member runs.
    # This is a no-op for configs without `_E#` names: the splice returns the
    # input path unchanged when there is nothing to assemble.
    base_params = CalibrationTools.parameter_dict(config_dict)
    splice(path) = CalibrationTools.write_spliced_parameter_file(
        path,
        base_params,
        path * ".spliced.toml",
    )

    center_toml = splice(write_param_toml(joinpath(scratch, "center.toml"), names, center))
    base = atmos_params_from_toml(CA, base_tomls, center_toml, FT, gw_flags)

    for (i, name) in enumerate(names)
        if name in unchecked
            warn!("wiring $name", "listed in PREFLIGHT_NON_ATMOS_PARAMS, not checked")
            continue
        end
        u = copy(u0)
        u[i] += su[i]
        perturbed = _asvec(PD.transform_unconstrained_to_constrained(Main.PRIORS, u))
        pert_toml = splice(
            write_param_toml(joinpath(scratch, "perturbed_$i.toml"), names, perturbed),
        )
        params = atmos_params_from_toml(CA, base_tomls, pert_toml, FT, gw_flags)
        diff = CalibrationTools.numeric_leaf_diff(base, params)
        if isempty(diff)
            fail!(
                "wiring $name",
                "perturbing $(round(center[i]; sigdigits = 4)) -> " *
                "$(round(perturbed[i]; sigdigits = 4)) changes nothing in the " *
                "ClimaAtmos parameters. The model ignores this parameter",
            )
        else
            shown = join(first(diff, 3), ", ")
            extra = length(diff) > 3 ? " and $(length(diff) - 3) more" : ""
            pass!("wiring $name", "reaches " * shown * extra)
        end
    end
    return nothing
end

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if abspath(PROGRAM_FILE) == @__FILE__
    args = parse_preflight_args(ARGS)
    # ClimaCoupler.Input parses the global command line through ArgParse and
    # exits on the unknown positional argument, so drop the arguments now that
    # they are parsed.
    empty!(ARGS)

    include(
        joinpath(
            pkgdir(ClimaCoupler),
            "experiments",
            "calibration",
            "amip",
            "preprocessing.jl",
        ),
    )
    @info "Preflight for $(args.config_path)"
    include(args.config_path)

    println("\n--- A. Config consistency ---")
    check_config(args.sypd)

    println("\n--- B. Observation vector ---")
    check_observations()

    println("\n--- C. Run directory ---")
    check_run_dir(args.config_path)

    println("\n--- D. Parameter wiring ---")
    if args.skip_wiring
        warn!("wiring", "skipped (--skip-wiring)")
    else
        @info "Loading ClimaAtmos for the wiring check"
        @eval import ClimaAtmos
        scratch = mktempdir()
        Base.invokelatest(check_wiring, ClimaAtmos, scratch)
    end

    n_fail = count(r -> r.status == "FAIL", PREFLIGHT_RESULTS)
    n_warn = count(r -> r.status == "WARN", PREFLIGHT_RESULTS)
    println(
        "\n=== Preflight: $(n_fail) failed, $(n_warn) warnings, " *
        "$(length(PREFLIGHT_RESULTS) - n_fail - n_warn) passed ===",
    )
    n_fail > 0 && exit(1)
end
