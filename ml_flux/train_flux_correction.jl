#!/usr/bin/env julia
#=
Offline training of a neural network flux correction for ClimaAtmos EDMF.

The NN predicts vertical flux profiles F(z) at cell faces for temperature and
moisture. Training minimises the MSE between -dF/dz and the residual tendency
(ERA5 - model), so the learned correction is conservative by construction.

Usage (CPU, single-threaded):
    julia --project=. train_flux_correction.jl
Usage (CPU, multi-threaded — recommended):
    julia -t 24 --project=. train_flux_correction.jl
Usage (GPU via SLURM):
    sbatch run_training.sbatch
Usage (interactive GPU session):
    see gpu_session.sh

Dependency note:
    WeightInitializers must be pinned to v1.0.5. Newer versions have a broken
    CUDA extension (cuRAND.RNG removed in CUDA.jl v6). After any Pkg.add or
    Pkg.update, re-pin with:
        Pkg.add(name="WeightInitializers", version="1.0.5")
=#

using NCDatasets
using Dates
using Statistics
using LinearAlgebra
using Printf
using Random
using Interpolations
using Lux
using Optimisers
using Zygote
using BSON
using SHA
using CairoMakie

# GPU support — set USE_GPU=true via environment variable to enable.
# On SLURM GPU jobs, the sbatch script sets this before launching.
# Note: we use LuxCUDA (which bundles cuDNN) and set local_toolkit=false in
# LocalPreferences.toml so cuDNN is loaded from Julia artifacts, not /usr/local/cuda.
const USE_GPU = get(ENV, "USE_GPU", "false") == "true" && let
    local gpu_avail = false
    try
        @eval using LuxCUDA
        gpu_avail = CUDA.functional()
    catch e
        @warn "CUDA loading failed, falling back to CPU" exception=e
    end
    gpu_avail
end

if USE_GPU
    println("GPU: ", CUDA.device(), " — ", CUDA.name(CUDA.device()))
    const DEV = gpu_device()
else
    println("Running on CPU")
    const DEV = cpu_device()
end

# ─────────────────────────────────────────────────────────────────────────────
# Configuration
# ─────────────────────────────────────────────────────────────────────────────
# const MODEL_DIR = "/home/cchristo/clima/ClimaCoupler.jl/output/wxquest_progedmf/output_0004/clima_atmos"
const MODEL_DIR = "/net/sampo/data1/cchristo/coupler_runs/copies/wxquest_progedmf/output_0008/clima_atmos"
const ERA5_DIR  = "/net/sampo/data1/wxquest_data/hourly_stats/ml_correct_v1"

const DT = 3600.0f0  # hourly data

const Z_MAX = 20_000.0f0  # only train on levels below this height (m)

const N_EPOCHS       = 100
const BATCH_SIZE     = 1024
const LEARNING_RATE  = 1f-3
const LR_MIN         = 1f-5    # cosine annealing floor
const TRAIN_FRACTION = 0.8
const HIDDEN_DIM     = 128
const PATIENCE       = 15     # early stopping: stop after this many epochs without improvement
const DROPOUT_RATE   = 0.1f0  # spatial dropout between conv layers

# const WEIGHT_DECAY       = 1f-3  # AdamW L2 regularisation
const WEIGHT_DECAY       = 3f-4  # AdamW L2 regularisation
const FLUX_SMOOTH_WEIGHT = 1f-4  # light penalty on flux curvature

const ARCHITECTURE       = :cnn    # :unet, :mlp, or :cnn
const TEMPORAL_AVERAGING = :daily  # :hourly or :daily
const PREDICT_MODE       = :direct # :flux (predict fluxes, loss on -dF/dz) or :direct (predict tendencies)
const TARGET_VARS        = :T      # :both, :T (temperature only), or :q (moisture only)

# ── Input features ───────────────────────────────────────────────────────────
# Each entry is (netcdf_varname, transform) where transform is applied per-element
# after loading. Add/remove entries here to change what the NN sees.
const INPUT_FEATURES = [
    ("ta",    identity),       # temperature (K) — the variable we're correcting
    ("hus",   identity),       # specific humidity (kg/kg) — the other target variable
    ("pfull", x -> log(max(x, 1f0))),  # log-pressure — encodes altitude / stratification
    # ("hur",   identity),       # relative humidity — key for cloud/convection regimes
    ("tke",   identity),       # turbulent kinetic energy — characterises SGS turbulence
    # ("arup",  identity),       # updraft area fraction — convective activity indicator
    # ("waup",  identity),       # updraft vertical velocity — convective intensity
]
const N_FEATURES = length(INPUT_FEATURES)

n_target_vars() = TARGET_VARS == :both ? 2 : 1

# ── Data subsetting ──────────────────────────────────────────────────────────
# Time window: only use model timesteps within this day range (1-indexed from
# the start of the run). Set to (1, Inf) to use all available timesteps.
const TIME_DAY_START = 3       # skip first 2 days (spinup)
const TIME_DAY_END   = 8      # through end of day 8

# Spatial subsampling: randomly keep this fraction of columns per timestep.
# Set to 1.0 to use all columns (no subsampling).
const COLUMN_SUBSAMPLE_FRACTION = 0.10   # 20% of 192×96 ≈ 3686 columns per step

const DATASET_CACHE_DIR = joinpath(@__DIR__, "cached_datasets")
const RUNS_DIR = joinpath(@__DIR__, "runs")

# ─────────────────────────────────────────────────────────────────────────────
# Run directory management
# ─────────────────────────────────────────────────────────────────────────────

function next_run_index()
    isdir(RUNS_DIR) || return 1
    existing = filter(d -> occursin(r"^\d{3}_", d), readdir(RUNS_DIR))
    isempty(existing) && return 1
    max_idx = maximum(parse(Int, m.match) for d in existing for m in eachmatch(r"^(\d{3})", d))
    return max_idx + 1
end

function make_run_dir()
    idx = next_run_index()
    tavg = TEMPORAL_AVERAGING == :daily ? "day" : "1h"
    name = @sprintf("%03d_%s_h%d_%s", idx, ARCHITECTURE, HIDDEN_DIM, tavg)
    path = joinpath(RUNS_DIR, name)
    mkpath(path)
    println("Run directory: $path")
    return path
end

# ─────────────────────────────────────────────────────────────────────────────
# Dataset caching — skip the expensive ERA5 regridding pipeline on repeat runs
# ─────────────────────────────────────────────────────────────────────────────

function dataset_cache_key()
    parts = string(
        MODEL_DIR, "|", ERA5_DIR, "|",
        Z_MAX, "|", TIME_DAY_START, "-", TIME_DAY_END, "|",
        COLUMN_SUBSAMPLE_FRACTION, "|",
        join([f[1] for f in INPUT_FEATURES], ","), "|",
        TEMPORAL_AVERAGING, "|", TARGET_VARS,
    )
    return bytes2hex(sha256(parts))
end

function save_dataset_cache(X, Y, dz, nz)
    mkpath(DATASET_CACHE_DIR)
    key = dataset_cache_key()
    path = joinpath(DATASET_CACHE_DIR, "dataset_$(key[1:12]).bson")
    BSON.@save path X Y dz nz
    println("  Cached dataset to $path")
end

function load_dataset_cache()
    key = dataset_cache_key()
    path = joinpath(DATASET_CACHE_DIR, "dataset_$(key[1:12]).bson")
    if isfile(path)
        println("  Loading cached dataset from $path")
        BSON.@load path X Y dz nz
        return X, Y, dz, nz
    end
    return nothing
end

# ─────────────────────────────────────────────────────────────────────────────
# Data loading helpers
# ─────────────────────────────────────────────────────────────────────────────

"""Load a 4D model variable (time, lon, lat, z) from the hourly-instantaneous file."""
function load_model_var(varname::String; dir = MODEL_DIR)
    path = joinpath(dir, "$(varname)_1h_inst.nc")
    NCDataset(path) do ds
        Float32.(Array(ds[varname]))  # (time, lon, lat, z)
    end
end

"""Load model grid info: lon, lat, z_reference, z_physical."""
function load_model_grid(; dir = MODEL_DIR)
    NCDataset(joinpath(dir, "ta_1h_inst.nc")) do ds
        lon = Float64.(Array(ds["lon"]))
        lat = Float64.(Array(ds["lat"]))
        z_ref = Float64.(Array(ds["z_reference"]))
        z_phys = Float32.(Array(ds["z_physical"]))  # (lon, lat, z)
        dates = Array(ds["date"])
        return (; lon, lat, z_ref, z_phys, dates)
    end
end

"""
Load one ERA5 pressure-level file matching a DateTime.
Returns (t, q, z_geopot, pressure_levels, era5_lon, era5_lat).
"""
function load_era5_pressure(dt::DateTime; dir = ERA5_DIR)
    fname = @sprintf("era5_pressure_levels_%04d%02d%02d_%02d00.nc",
        year(dt), month(dt), day(dt), hour(dt))
    path = joinpath(dir, fname)
    if !isfile(path)
        return nothing
    end
    NCDataset(path) do ds
        t   = Float32.(nomissing(Array(ds["t"][:, :, :, 1]), NaN32))   # (lon, lat, plev)
        q   = Float32.(nomissing(Array(ds["q"][:, :, :, 1]), NaN32))
        z   = Float32.(nomissing(Array(ds["z"][:, :, :, 1]), NaN32))   # geopotential m²/s²
        plev = Float64.(nomissing(Array(ds["pressure_level"])))        # hPa
        elon = Float64.(nomissing(Array(ds["longitude"])))
        elat = Float64.(nomissing(Array(ds["latitude"])))
        return (; t, q, z, plev, lon = elon, lat = elat)
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Regridding: ERA5 → model horizontal grid (bilinear interpolation)
# ─────────────────────────────────────────────────────────────────────────────

"""
Bilinear interpolation of a 3D field (lon, lat, plev) from ERA5 to model grid.
Handles the periodic longitude wrapping.
"""
function regrid_horizontal(field_era5::Array{Float32, 3},
                           era5_lon::Vector{Float64}, era5_lat::Vector{Float64},
                           model_lon::Vector{Float64}, model_lat::Vector{Float64})
    nlev = size(field_era5, 3)
    nlon_m = length(model_lon)
    nlat_m = length(model_lat)
    out = zeros(Float32, nlon_m, nlat_m, nlev)

    era5_lat_sorted = sort(era5_lat)
    flip_lat = era5_lat[1] > era5_lat[end]

    Threads.@threads for k in 1:nlev
        slab = @view field_era5[:, :, k]
        if flip_lat
            slab_sorted = slab[:, end:-1:1]
        else
            slab_sorted = slab
        end
        itp = interpolate((era5_lon, era5_lat_sorted), slab_sorted,
                          (Gridded(Linear()), Gridded(Linear())))
        etp = extrapolate(itp, (Periodic(), Flat()))
        for j in 1:nlat_m, i in 1:nlon_m
            out[i, j, k] = etp(model_lon[i], model_lat[j])
        end
    end
    return out
end

# ─────────────────────────────────────────────────────────────────────────────
# Vertical interpolation: ERA5 pressure levels → model z levels
# ─────────────────────────────────────────────────────────────────────────────

"""
For a single column, interpolate a profile from ERA5 pressure levels to model
levels using geopotential height as the common coordinate.

era5_z_col: geopotential heights on ERA5 pressure levels (m), length nplev
era5_var_col: variable values on ERA5 pressure levels, length nplev
model_z_col: physical heights on model levels (m), length nz_model

ERA5 geopotential → height: z_m = z_geopot / 9.80665
"""
function interp_vertical_column(era5_z_col::AbstractVector{Float32},
                                era5_var_col::AbstractVector{Float32},
                                model_z_col::AbstractVector{Float32})
    nz = length(model_z_col)

    mask = .!isnan.(era5_z_col) .& .!isnan.(era5_var_col)
    z_good = era5_z_col[mask]
    v_good = era5_var_col[mask]

    if length(z_good) < 2
        return fill(NaN32, nz)
    end

    perm = sortperm(z_good)
    z_sorted = z_good[perm]
    v_sorted = v_good[perm]

    itp = interpolate((z_sorted,), v_sorted, Gridded(Linear()))
    etp = extrapolate(itp, Flat())

    return Float32[etp(model_z_col[k]) for k in 1:nz]
end

"""
Interpolate a full 3D ERA5 field (already regridded horizontally) from ERA5
pressure levels to model z levels.

field_era5_regridded: (nlon, nlat, nplev) on model horiz grid
era5_z_regridded:     (nlon, nlat, nplev) geopotential/g heights on model horiz grid
model_z_phys:         (nlon, nlat, nz_model) physical heights
"""
function interp_vertical_3d(field_era5::Array{Float32, 3},
                            era5_z::Array{Float32, 3},
                            model_z_phys::Array{Float32, 3})
    nlon, nlat, _ = size(field_era5)
    nz_model = size(model_z_phys, 3)
    out = zeros(Float32, nlon, nlat, nz_model)
    col_indices = [(i, j) for j in 1:nlat for i in 1:nlon]
    Threads.@threads for idx in col_indices
        i, j = idx
        era5_z_col = era5_z[i, j, :] ./ 9.80665f0
        era5_v_col = field_era5[i, j, :]
        model_z_col = model_z_phys[i, j, :]
        out[i, j, :] = interp_vertical_column(era5_z_col, era5_v_col, model_z_col)
    end
    return out
end

# ─────────────────────────────────────────────────────────────────────────────
# Build training dataset
# ─────────────────────────────────────────────────────────────────────────────

function build_dataset()
    if TEMPORAL_AVERAGING == :daily
        return build_dataset_daily_mean()
    else
        return build_dataset_hourly()
    end
end

"""
Build the full training dataset (hourly instantaneous mode):
  - features X: (n_features × n_levels, n_columns)
  - targets  Y: (2 × n_levels, n_columns)   [T tendency, q tendency]
  - dz:         (n_levels, n_columns)
"""
function build_dataset_hourly()
    rng_data = Random.MersenneTwister(123)

    println("Loading model grid...")
    grid = load_model_grid()
    nlon = length(grid.lon)
    nlat = length(grid.lat)
    nz_full = length(grid.z_ref)
    model_dates = grid.dates
    nt = length(model_dates)

    z_mask = grid.z_ref .<= Z_MAX
    z_idx  = findall(z_mask)
    nz     = length(z_idx)

    # Time subsetting: convert day range to hourly timestep indices
    t0 = model_dates[1]
    ti_start = max(1, findfirst(d -> (d - t0).value / (3600_000) >= 24 * (TIME_DAY_START - 1), model_dates))
    ti_end   = something(findlast(d -> (d - t0).value / (3600_000) < 24 * TIME_DAY_END, model_dates), nt)
    # Need ti+1, so last usable pair index is ti_end-1
    ti_end = min(ti_end, nt) - 1

    println("  Grid: $(nlon) lon × $(nlat) lat × $(nz_full) z ($(nz) below $(Z_MAX)m), $(nt) timesteps")
    println("  Date range: $(model_dates[1]) to $(model_dates[end])")
    println("  Time subset: steps $(ti_start)-$(ti_end) ($(model_dates[ti_start]) to $(model_dates[ti_end]))")
    println("  Column subsample: $(COLUMN_SUBSAMPLE_FRACTION * 100)% of $(nlon * nlat) = ~$(round(Int, COLUMN_SUBSAMPLE_FRACTION * nlon * nlat)) columns/step")

    feat_names = [f[1] for f in INPUT_FEATURES]

    println("Loading model variables...")
    println("  Input features: ", join(feat_names, ", "), " ($N_FEATURES total)")

    # Load all unique variables needed (features + ta/hus for targets)
    vars_needed = unique(vcat(feat_names, ["ta", "hus"]))
    feat_data = Dict{String, Array{Float32, 4}}()
    for name in vars_needed
        println("    loading $name ...")
        feat_data[name] = load_model_var(name)
    end
    ta_model  = feat_data["ta"]
    hus_model = feat_data["hus"]

    n_valid_pairs = 0
    all_X = Vector{Matrix{Float32}}()
    all_Y = Vector{Matrix{Float32}}()
    all_dz = Vector{Matrix{Float32}}()

    # Build (i,j) column index list for subsampling
    all_col_ij = [(i, j) for j in 1:nlat for i in 1:nlon]
    ncols_total = length(all_col_ij)
    ncols_sample = max(1, round(Int, COLUMN_SUBSAMPLE_FRACTION * ncols_total))

    for ti in ti_start:ti_end
        dt_now  = model_dates[ti]
        dt_next = model_dates[ti + 1]

        era5_now  = load_era5_pressure(dt_now)
        era5_next = load_era5_pressure(dt_next)
        if era5_now === nothing || era5_next === nothing
            continue
        end

        if ti == ti_start || ti % 24 == 0
            println("  Processing timestep $ti / $ti_end: $dt_now")
            flush(stdout)
        end

        t_now  = regrid_horizontal(era5_now.t,  era5_now.lon, era5_now.lat,
                                   grid.lon, grid.lat)
        t_next = regrid_horizontal(era5_next.t, era5_next.lon, era5_next.lat,
                                   grid.lon, grid.lat)
        q_now  = regrid_horizontal(era5_now.q,  era5_now.lon, era5_now.lat,
                                   grid.lon, grid.lat)
        q_next = regrid_horizontal(era5_next.q, era5_next.lon, era5_next.lat,
                                   grid.lon, grid.lat)
        z_now  = regrid_horizontal(era5_now.z,  era5_now.lon, era5_now.lat,
                                   grid.lon, grid.lat)

        t_era5_now  = interp_vertical_3d(t_now,  z_now, grid.z_phys)
        t_era5_next = interp_vertical_3d(t_next, z_now, grid.z_phys)
        q_era5_now  = interp_vertical_3d(q_now,  z_now, grid.z_phys)
        q_era5_next = interp_vertical_3d(q_next, z_now, grid.z_phys)

        # Slice to levels below Z_MAX
        dT_era5 = (t_era5_next[:, :, z_idx] .- t_era5_now[:, :, z_idx]) ./ DT
        dq_era5 = (q_era5_next[:, :, z_idx] .- q_era5_now[:, :, z_idx]) ./ DT

        ta_now_z  = ta_model[ti, :, :, z_idx]
        ta_next_z = ta_model[ti + 1, :, :, z_idx]
        hus_now_z  = hus_model[ti, :, :, z_idx]
        hus_next_z = hus_model[ti + 1, :, :, z_idx]

        dT_model = (ta_next_z .- ta_now_z) ./ DT
        dq_model = (hus_next_z .- hus_now_z) ./ DT

        dT_corr = dT_era5 .- dT_model
        dq_corr = dq_era5 .- dq_model

        # Randomly sample columns for this timestep
        col_sample = if COLUMN_SUBSAMPLE_FRACTION >= 1.0
            all_col_ij
        else
            all_col_ij[randperm(rng_data, ncols_total)[1:ncols_sample]]
        end
        nc = length(col_sample)

        nv = n_target_vars()
        X_step = zeros(Float32, N_FEATURES * nz, nc)
        Y_step = zeros(Float32, nv * nz, nc)
        dz_step = zeros(Float32, nz, nc)

        Threads.@threads for col in 1:nc
            i, j = col_sample[col]
            z_col = grid.z_phys[i, j, z_idx]

            for k in 1:nz
                if k < nz
                    dz_step[k, col] = z_col[k + 1] - z_col[k]
                else
                    dz_step[k, col] = dz_step[k - 1, col]
                end
            end

            for (fi, (name, transform)) in enumerate(INPUT_FEATURES)
                offset = (fi - 1) * nz
                raw = feat_data[name][ti, i, j, z_idx]
                X_step[offset+1:offset+nz, col] = transform.(raw)
            end

            if TARGET_VARS == :T || TARGET_VARS == :both
                Y_step[1:nz, col] = dT_corr[i, j, :]
            end
            if TARGET_VARS == :q
                Y_step[1:nz, col] = dq_corr[i, j, :]
            elseif TARGET_VARS == :both
                Y_step[nz+1:2*nz, col] = dq_corr[i, j, :]
            end
        end

        push!(all_X, X_step)
        push!(all_Y, Y_step)
        push!(all_dz, dz_step)
        n_valid_pairs += 1
    end

    println("  Built $n_valid_pairs valid timestep pairs, total columns: $(sum(size.(all_X, 2)))")

    X = hcat(all_X...)
    Y = hcat(all_Y...)
    dz = hcat(all_dz...)

    nan_cols = vec(any(isnan.(X), dims = 1) .| any(isnan.(Y), dims = 1))
    good = .!nan_cols
    X = X[:, good]
    Y = Y[:, good]
    dz = dz[:, good]
    println("  After NaN filtering: $(size(X, 2)) columns remain")

    return X, Y, dz, nz
end

"""
Build training dataset using daily-mean averaging.

For each day in the time window, averages both model state inputs and
ERA5-minus-model tendency corrections over all valid hourly pairs within
that day. This removes gravity waves, tides, and diurnal cycle noise,
leaving the mean bias signal that a column NN can learn.

Returns the same format as build_dataset_hourly():
  - features X: (n_features × n_levels, n_columns)
  - targets  Y: (2 × n_levels, n_columns)
  - dz:         (n_levels, n_columns)
"""
function build_dataset_daily_mean()
    rng_data = Random.MersenneTwister(123)

    println("Loading model grid...")
    grid = load_model_grid()
    nlon = length(grid.lon)
    nlat = length(grid.lat)
    nz_full = length(grid.z_ref)
    model_dates = grid.dates
    nt = length(model_dates)

    z_mask = grid.z_ref .<= Z_MAX
    z_idx  = findall(z_mask)
    nz     = length(z_idx)

    # Time subsetting
    t0 = model_dates[1]
    ti_start = max(1, findfirst(d -> (d - t0).value / (3600_000) >= 24 * (TIME_DAY_START - 1), model_dates))
    ti_end   = something(findlast(d -> (d - t0).value / (3600_000) < 24 * TIME_DAY_END, model_dates), nt)
    ti_end = min(ti_end, nt) - 1

    # Group timestep indices into days
    day_groups = Dict{Date, Vector{Int}}()
    for ti in ti_start:ti_end
        d = Date(model_dates[ti])
        if !haskey(day_groups, d)
            day_groups[d] = Int[]
        end
        push!(day_groups[d], ti)
    end
    days_sorted = sort(collect(keys(day_groups)))

    println("  Grid: $(nlon) lon × $(nlat) lat × $(nz_full) z ($(nz) below $(Z_MAX)m), $(nt) timesteps")
    println("  Date range: $(model_dates[1]) to $(model_dates[end])")
    println("  Daily-mean mode: $(length(days_sorted)) days ($(days_sorted[1]) to $(days_sorted[end]))")

    col_subsample = TEMPORAL_AVERAGING == :daily ? max(COLUMN_SUBSAMPLE_FRACTION, 0.25) : COLUMN_SUBSAMPLE_FRACTION
    println("  Column subsample: $(col_subsample * 100)% of $(nlon * nlat) = ~$(round(Int, col_subsample * nlon * nlat)) columns/day")

    feat_names = [f[1] for f in INPUT_FEATURES]

    println("Loading model variables...")
    println("  Input features: ", join(feat_names, ", "), " ($N_FEATURES total)")

    vars_needed = unique(vcat(feat_names, ["ta", "hus"]))
    feat_data = Dict{String, Array{Float32, 4}}()
    for name in vars_needed
        println("    loading $name ...")
        feat_data[name] = load_model_var(name)
    end
    ta_model  = feat_data["ta"]
    hus_model = feat_data["hus"]

    all_col_ij = [(i, j) for j in 1:nlat for i in 1:nlon]
    ncols_total = length(all_col_ij)
    ncols_sample = max(1, round(Int, col_subsample * ncols_total))

    n_valid_days = 0
    all_X = Vector{Matrix{Float32}}()
    all_Y = Vector{Matrix{Float32}}()
    all_dz = Vector{Matrix{Float32}}()

    for day in days_sorted
        hours = day_groups[day]
        println("  Processing day $day ($(length(hours)) hourly pairs)")
        flush(stdout)

        # Accumulators over the full grid
        X_accum  = zeros(Float32, N_FEATURES, nlon, nlat, nz)
        dT_accum = zeros(Float32, nlon, nlat, nz)
        dq_accum = zeros(Float32, nlon, nlat, nz)
        n_hours_valid = 0

        for ti in hours
            dt_now  = model_dates[ti]
            dt_next = model_dates[ti + 1]

            era5_now  = load_era5_pressure(dt_now)
            era5_next = load_era5_pressure(dt_next)
            if era5_now === nothing || era5_next === nothing
                continue
            end

            t_now  = regrid_horizontal(era5_now.t,  era5_now.lon, era5_now.lat, grid.lon, grid.lat)
            t_next = regrid_horizontal(era5_next.t, era5_next.lon, era5_next.lat, grid.lon, grid.lat)
            q_now  = regrid_horizontal(era5_now.q,  era5_now.lon, era5_now.lat, grid.lon, grid.lat)
            q_next = regrid_horizontal(era5_next.q, era5_next.lon, era5_next.lat, grid.lon, grid.lat)
            z_now  = regrid_horizontal(era5_now.z,  era5_now.lon, era5_now.lat, grid.lon, grid.lat)

            t_era5_now  = interp_vertical_3d(t_now,  z_now, grid.z_phys)
            t_era5_next = interp_vertical_3d(t_next, z_now, grid.z_phys)
            q_era5_now  = interp_vertical_3d(q_now,  z_now, grid.z_phys)
            q_era5_next = interp_vertical_3d(q_next, z_now, grid.z_phys)

            dT_era5 = (t_era5_next[:, :, z_idx] .- t_era5_now[:, :, z_idx]) ./ DT
            dq_era5 = (q_era5_next[:, :, z_idx] .- q_era5_now[:, :, z_idx]) ./ DT

            ta_now_z   = ta_model[ti, :, :, z_idx]
            ta_next_z  = ta_model[ti + 1, :, :, z_idx]
            hus_now_z  = hus_model[ti, :, :, z_idx]
            hus_next_z = hus_model[ti + 1, :, :, z_idx]

            dT_model = (ta_next_z .- ta_now_z) ./ DT
            dq_model = (hus_next_z .- hus_now_z) ./ DT

            dT_accum .+= (dT_era5 .- dT_model)
            dq_accum .+= (dq_era5 .- dq_model)

            # Accumulate input features
            for (fi, (name, transform)) in enumerate(INPUT_FEATURES)
                raw = feat_data[name][ti, :, :, z_idx]
                X_accum[fi, :, :, :] .+= transform.(raw)
            end

            n_hours_valid += 1
        end

        if n_hours_valid == 0
            continue
        end

        # Average over valid hours
        inv_n = 1.0f0 / n_hours_valid
        X_accum  .*= inv_n
        dT_accum .*= inv_n
        dq_accum .*= inv_n

        # Subsample columns
        col_sample = if col_subsample >= 1.0
            all_col_ij
        else
            all_col_ij[randperm(rng_data, ncols_total)[1:ncols_sample]]
        end
        nc = length(col_sample)

        nv = n_target_vars()
        X_day  = zeros(Float32, N_FEATURES * nz, nc)
        Y_day  = zeros(Float32, nv * nz, nc)
        dz_day = zeros(Float32, nz, nc)

        Threads.@threads for col in 1:nc
            i, j = col_sample[col]
            z_col = grid.z_phys[i, j, z_idx]

            for k in 1:nz
                if k < nz
                    dz_day[k, col] = z_col[k + 1] - z_col[k]
                else
                    dz_day[k, col] = dz_day[k - 1, col]
                end
            end

            for fi in 1:N_FEATURES
                offset = (fi - 1) * nz
                X_day[offset+1:offset+nz, col] = X_accum[fi, i, j, :]
            end

            if TARGET_VARS == :T || TARGET_VARS == :both
                Y_day[1:nz, col] = dT_accum[i, j, :]
            end
            if TARGET_VARS == :q
                Y_day[1:nz, col] = dq_accum[i, j, :]
            elseif TARGET_VARS == :both
                Y_day[nz+1:2*nz, col] = dq_accum[i, j, :]
            end
        end

        push!(all_X, X_day)
        push!(all_Y, Y_day)
        push!(all_dz, dz_day)
        n_valid_days += 1
    end

    println("  Built $n_valid_days valid days, total columns: $(sum(size.(all_X, 2)))")

    X = hcat(all_X...)
    Y = hcat(all_Y...)
    dz = hcat(all_dz...)

    nan_cols = vec(any(isnan.(X), dims = 1) .| any(isnan.(Y), dims = 1))
    good = .!nan_cols
    X = X[:, good]
    Y = Y[:, good]
    dz = dz[:, good]
    println("  After NaN filtering: $(size(X, 2)) columns remain")

    return X, Y, dz, nz
end

# ─────────────────────────────────────────────────────────────────────────────
# Normalisation
# ─────────────────────────────────────────────────────────────────────────────

struct NormStats
    mean::Vector{Float32}
    std::Vector{Float32}
end

"""Per-dimension normalisation for input features."""
function compute_norm(X::Matrix{Float32})
    m = vec(mean(X, dims = 2))
    s = vec(std(X, dims = 2))
    s[s .< 1f-8] .= 1f0
    NormStats(m, s)
end

"""
Per-variable normalisation for targets: one (mean, std) per target variable
across all levels. This preserves vertical structure while bringing each to O(1).
"""
function compute_target_norm(Y::Matrix{Float32}, nz::Int)
    nv = n_target_vars()
    m = Float32[]
    s = Float32[]
    for v in 1:nv
        rows = (v-1)*nz+1 : v*nz
        push!(m, fill(Float32(mean(Y[rows, :])), nz)...)
        push!(s, fill(Float32(std(Y[rows, :])), nz)...)
    end
    NormStats(m, s)
end

normalise(X, ns::NormStats)   = (X .- ns.mean) ./ ns.std
denormalise(X, ns::NormStats) = X .* ns.std .+ ns.mean

# ─────────────────────────────────────────────────────────────────────────────
# Neural network architectures
# ─────────────────────────────────────────────────────────────────────────────

"""
Build the correction NN. Dispatches on ARCHITECTURE and PREDICT_MODE.

Input:  (n_features × nz, batch)
Output: depends on PREDICT_MODE:
  :flux   → (2*(nz-1), batch) — interior face fluxes for T and q
  :direct → (2*nz, batch)     — tendency corrections for T and q at cell centres
"""
function build_model(input_dim::Int, nz::Int; predict_mode = PREDICT_MODE)
    if ARCHITECTURE == :unet
        return build_unet_model(input_dim, nz; predict_mode)
    elseif ARCHITECTURE == :mlp
        return build_mlp_model(input_dim, nz; predict_mode)
    elseif ARCHITECTURE == :cnn
        return build_cnn_model(input_dim, nz; predict_mode)
    else
        error("Unknown ARCHITECTURE: $ARCHITECTURE. Use :unet, :mlp, or :cnn")
    end
end

function build_mlp_model(input_dim::Int, nz::Int; predict_mode = PREDICT_MODE)
    nv = n_target_vars()
    output_dim = predict_mode == :direct ? nv * nz : nv * (nz - 1)

    model = Chain(
        Dense(input_dim => HIDDEN_DIM, gelu),
        Dense(HIDDEN_DIM => HIDDEN_DIM, gelu),
        Dense(HIDDEN_DIM => HIDDEN_DIM, gelu),
        Dense(HIDDEN_DIM => output_dim),
    )
    return model
end

"""
1D CNN over the vertical column — no pooling, no skip connections.

Stacked Conv1d layers with increasing dilation to grow the receptive field
without losing resolution. The final layer maps to 2 channels (T, q) at
cell centres, then centres are averaged to interior faces.

With 6 layers of kernel=3 and dilations [1,1,2,4,1,1], the receptive field
is 19 levels — enough that every level sees ~40% of the 49-level column,
and the network remains simple and robust.
"""
function build_cnn_model(input_dim::Int, nz::Int; predict_mode = PREDICT_MODE)
    n_feat = N_FEATURES
    nv = n_target_vars()
    @assert input_dim == n_feat * nz "input_dim ($input_dim) != N_FEATURES ($n_feat) * nz ($nz)"

    nf = nz - 1
    ch = max(HIDDEN_DIM ÷ 4, 16)
    direct = predict_mode == :direct

    dilations = [1, 1, 2, 4, 1, 1]

    model = @compact(
        conv1 = Conv((3,), n_feat => ch, gelu; pad = 1 * dilations[1], dilation = dilations[1]),
        conv2 = Conv((3,), ch => ch, gelu; pad = 1 * dilations[2], dilation = dilations[2]),
        drop1 = Dropout(DROPOUT_RATE),
        conv3 = Conv((3,), ch => 2ch, gelu; pad = 1 * dilations[3], dilation = dilations[3]),
        conv4 = Conv((3,), 2ch => 2ch, gelu; pad = 1 * dilations[4], dilation = dilations[4]),
        drop2 = Dropout(DROPOUT_RATE),
        conv5 = Conv((3,), 2ch => ch, gelu; pad = 1 * dilations[5], dilation = dilations[5]),
        conv6 = Conv((3,), ch => ch, gelu; pad = 1 * dilations[6], dilation = dilations[6]),
        out_conv = Conv((1,), ch => nv),
    ) do x
        B = size(x, 2)

        h = reshape(x, nz, n_feat, B)

        h = conv1(h)
        h = conv2(h)
        h = drop1(h)
        h = conv3(h)
        h = conv4(h)
        h = drop2(h)
        h = conv5(h)
        h = conv6(h)

        out = out_conv(h)  # (nz, nv, B)

        if direct
            @return vcat([reshape(out[:, v:v, :], nz, B) for v in 1:nv]...)
        else
            faces = (out[1:nf, :, :] .+ out[2:nz, :, :]) ./ 2f0
            @return vcat([reshape(faces[:, v:v, :], nf, B) for v in 1:nv]...)
        end
    end

    return model
end

"""
1D UNet over the vertical column.

Encoder:  two downsampling stages via Conv + MaxPool
Decoder:  two upsampling stages via ConvTranspose + skip connections
Output:   Conv to 2 channels (T, q) at cell centres, then interpolated to faces

The vertical dimension is padded to the next multiple of 4 for clean pooling,
then cropped back after the decoder.
"""
function build_unet_model(input_dim::Int, nz::Int; predict_mode = PREDICT_MODE)
    n_feat = N_FEATURES
    nv = n_target_vars()
    @assert input_dim == n_feat * nz "input_dim ($input_dim) != N_FEATURES ($n_feat) * nz ($nz)"

    nf = nz - 1
    nz_pad = nz % 4 == 0 ? nz : nz + (4 - nz % 4)
    ch = max(HIDDEN_DIM ÷ 4, 16)
    direct = predict_mode == :direct

    model = @compact(
        enc1 = Chain(Conv((3,), n_feat => ch, gelu; pad = 1),
                     Conv((3,), ch => ch, gelu; pad = 1)),
        enc2 = Chain(Conv((3,), ch => 2ch, gelu; pad = 1),
                     Conv((3,), 2ch => 2ch, gelu; pad = 1)),
        bneck = Chain(Conv((3,), 2ch => 4ch, gelu; pad = 1),
                      Conv((3,), 4ch => 4ch, gelu; pad = 1)),
        up2     = ConvTranspose((2,), 4ch => 2ch; stride = 2),
        dec2    = Chain(Conv((3,), 4ch => 2ch, gelu; pad = 1),
                        Conv((3,), 2ch => 2ch, gelu; pad = 1)),
        up1     = ConvTranspose((2,), 2ch => ch; stride = 2),
        dec1    = Chain(Conv((3,), 2ch => ch, gelu; pad = 1),
                        Conv((3,), ch => ch, gelu; pad = 1)),
        out_conv = Conv((1,), ch => nv),
        pool     = MaxPool((2,)),
    ) do x
        B = size(x, 2)

        x3d = reshape(x, nz, n_feat, B)

        if nz_pad > nz
            x3d = cat(x3d, x3d[1:nz_pad - nz, :, :] .* 0f0; dims = 1)
        end

        e1 = enc1(x3d)
        e2 = enc2(pool(e1))
        b  = bneck(pool(e2))

        d2 = dec2(cat(up2(b), e2; dims = 2))
        d1 = dec1(cat(up1(d2), e1; dims = 2))

        out = out_conv(d1)  # (nz_pad, nv, B)
        c = out[1:nz, :, :]  # (nz, nv, B)

        if direct
            @return vcat([reshape(c[:, v:v, :], nz, B) for v in 1:nv]...)
        else
            faces = (c[1:nf, :, :] .+ c[2:nz, :, :]) ./ 2f0
            @return vcat([reshape(faces[:, v:v, :], nf, B) for v in 1:nv]...)
        end
    end

    return model
end

# ─────────────────────────────────────────────────────────────────────────────
# Loss functions
# ─────────────────────────────────────────────────────────────────────────────

"""Direct tendency MSE loss: prediction and target are both (2*nz, batch)."""
function direct_tendency_loss(pred, target_norm, dz, nz)
    return mean((pred .- target_norm) .^ 2)
end

"""
Flux divergence loss: predict face fluxes, compare -dF/dz to target tendency.

flux_interior: (nv*(nz-1), batch)
target_norm:   (nv*nz, batch)
dz:            (nz, batch)
"""
function flux_divergence_loss(flux_interior, target_norm, dz, nz)
    nf = nz - 1
    nv = n_target_vars()
    z_row = flux_interior[1:1, :] .* 0f0
    total_loss = 0f0
    smooth = 0f0

    for v in 1:nv
        F_int = flux_interior[(v-1)*nf+1 : v*nf, :]
        F_full = vcat(z_row, F_int, z_row)
        dFdz = -(F_full[2:nz+1, :] .- F_full[1:nz, :]) ./ dz
        target_v = target_norm[(v-1)*nz+1 : v*nz, :]
        total_loss += mean((dFdz .- target_v) .^ 2)

        if nf >= 3
            d2F = F_int[1:nf-2, :] .- 2f0 .* F_int[2:nf-1, :] .+ F_int[3:nf, :]
            smooth += mean(d2F .^ 2)
        end
    end

    return total_loss + FLUX_SMOOTH_WEIGHT * smooth
end

function compute_loss(pred, target_norm, dz, nz)
    if PREDICT_MODE == :direct
        return direct_tendency_loss(pred, target_norm, dz, nz)
    else
        return flux_divergence_loss(pred, target_norm, dz, nz)
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Live training progress plots
# ─────────────────────────────────────────────────────────────────────────────

function save_loss_plot(train_losses, val_losses, path)
    fig = Figure(size = (700, 400))
    ax = Axis(fig[1, 1]; xlabel = "Epoch", ylabel = "Loss",
              title = "Training progress", yscale = log10)
    lines!(ax, 1:length(train_losses), train_losses; label = "Train", linewidth = 2)
    lines!(ax, 1:length(val_losses), val_losses; label = "Val", linewidth = 2)
    axislegend(ax; position = :rt)
    save(path, fig; px_per_unit = 2)
end

"""
Save snapshot of predicted vs target profiles for a few validation columns.
Adapts to TARGET_VARS — shows whichever variables are being trained.
"""
function save_profile_snapshot(model, ps, st, X_val, Y_val, dz_val,
                               nz, Y_norm_stats, z_km, epoch, path)
    rng_plot = Random.MersenneTwister(77)
    n_val = size(X_val, 2)
    n_cols = min(4, n_val)
    col_idx = randperm(rng_plot, n_val)[1:n_cols]

    Xb = X_val[:, col_idx]
    Yb = Y_val[:, col_idx]
    dzb = dz_val[:, col_idx]

    Xb_cpu = Xb |> cpu_device()
    Yb_cpu = Yb |> cpu_device()
    dzb_cpu = dzb |> cpu_device()

    st_test = Lux.testmode(st)
    raw_pred, _ = Lux.apply(model, Xb, ps, st_test)
    raw_pred_cpu = raw_pred |> cpu_device()

    nv = n_target_vars()

    # Convert model output to tendency predictions in normalised space per variable
    pred_norms = Vector{Matrix{Float32}}(undef, nv)
    if PREDICT_MODE == :direct
        for v in 1:nv
            pred_norms[v] = raw_pred_cpu[(v-1)*nz+1 : v*nz, :]
        end
    else
        nf = nz - 1
        for v in 1:nv
            F_int = raw_pred_cpu[(v-1)*nf+1 : v*nf, :]
            F_full = vcat(zeros(Float32, 1, n_cols), F_int, zeros(Float32, 1, n_cols))
            pred_norms[v] = -(F_full[2:nz+1, :] .- F_full[1:nz, :]) ./ dzb_cpu
        end
    end

    var_labels = TARGET_VARS == :both ? ["T", "q"] :
                 TARGET_VARS == :T   ? ["T"] : ["q"]
    var_units  = TARGET_VARS == :both ? ["K/hr", "g/kg/hr"] :
                 TARGET_VARS == :T   ? ["K/hr"] : ["g/kg/hr"]
    var_scales = TARGET_VARS == :both ? [3600f0, 3600f0 * 1000f0] :
                 TARGET_VARS == :T   ? [3600f0] : [3600f0 * 1000f0]

    fig = Figure(size = (300 * n_cols, 300 * nv))
    Label(fig[0, :], "Epoch $epoch: predicted vs target tendency"; fontsize = 14)

    for (p, ci) in enumerate(1:n_cols)
        for v in 1:nv
            std_v  = Y_norm_stats.std[(v-1)*nz + 1]
            mean_v = Y_norm_stats.mean[(v-1)*nz + 1]
            scale  = var_scales[v]

            pred_phys = (pred_norms[v][:, ci] .* std_v .+ mean_v) .* scale
            tgt_phys  = (Yb_cpu[(v-1)*nz+1 : v*nz, ci] .* std_v .+ mean_v) .* scale

            ax = Axis(fig[v, p]; xlabel = var_units[v],
                       ylabel = p == 1 ? "Height (km)" : "",
                       title = "Col $ci — $(var_labels[v])")
            lines!(ax, tgt_phys, z_km; label = "Target", linewidth = 1.5)
            lines!(ax, pred_phys, z_km; label = "Predicted", linewidth = 1.5, linestyle = :dash)
            v == 1 && p == 1 && axislegend(ax; position = :rt, labelsize = 9)
        end
    end

    save(path, fig; px_per_unit = 2)
end

function get_z_km(nz)
    grid = load_model_grid()
    z_mask = grid.z_ref .<= Z_MAX
    z_idx = findall(z_mask)
    return Float32.(grid.z_ref[z_idx] ./ 1000.0)
end

# ─────────────────────────────────────────────────────────────────────────────
# Training loop
# ─────────────────────────────────────────────────────────────────────────────

function train()
    rng = Random.MersenneTwister(42)

    println("=" ^ 60)
    println("Building dataset...")
    println("=" ^ 60)

    cached = load_dataset_cache()
    X_raw, Y_raw, dz_all, nz = if cached !== nothing
        cached
    else
        result = build_dataset()
        save_dataset_cache(result...)
        result
    end
    n_samples = size(X_raw, 2)
    input_dim = size(X_raw, 1)
    println("Dataset: $n_samples columns, input_dim=$input_dim, nz=$nz")

    println("\nComputing normalisation statistics...")
    X_norm_stats = compute_norm(X_raw)
    Y_norm_stats = compute_target_norm(Y_raw, nz)

    X = normalise(X_raw, X_norm_stats)
    Y = normalise(Y_raw, Y_norm_stats)

    nv = n_target_vars()
    var_names = TARGET_VARS == :both ? ["T", "q"] : TARGET_VARS == :T ? ["T"] : ["q"]
    var_unit  = TARGET_VARS == :both ? ["K/s", "kg/kg/s"] : TARGET_VARS == :T ? ["K/s"] : ["kg/kg/s"]
    for v in 1:nv
        @printf("  target %s: mean=%.4e  std=%.4e (%s)\n",
                var_names[v], Y_norm_stats.mean[(v-1)*nz+1], Y_norm_stats.std[(v-1)*nz+1], var_unit[v])
    end
    println("  Target variables normalised to O(1)")

    n_train = round(Int, TRAIN_FRACTION * n_samples)
    perm = randperm(rng, n_samples)
    train_idx = perm[1:n_train]
    val_idx   = perm[n_train+1:end]

    X_train, Y_train, dz_train = X[:, train_idx], Y[:, train_idx], dz_all[:, train_idx]
    X_val,   Y_val,   dz_val   = X[:, val_idx],   Y[:, val_idx],   dz_all[:, val_idx]

    println("Train: $(size(X_train, 2)) columns, Val: $(size(X_val, 2)) columns")

    println("\nBuilding model...")
    model = build_model(input_dim, nz)
    ps, st = Lux.setup(rng, model)

    ps = ps |> DEV
    st = st |> DEV
    opt_state = Optimisers.setup(AdamW(LEARNING_RATE, (0.9f0, 0.999f0), WEIGHT_DECAY), ps)

    cosine_lr(epoch, T) = LR_MIN + 0.5f0 * (LEARNING_RATE - LR_MIN) * (1 + cos(Float32(π) * epoch / T))

    # Move ALL data to device once — dataset is small enough to fit in VRAM
    X_train_d  = X_train  |> DEV
    Y_train_d  = Y_train  |> DEV
    dz_train_d = dz_train |> DEV
    X_val_d    = X_val    |> DEV
    Y_val_d    = Y_val    |> DEV
    dz_val_d   = dz_val   |> DEV

    # Create run directory now so we can save progress during training
    run_dir = make_run_dir()
    plots_dir = joinpath(run_dir, "plots")
    mkpath(plots_dir)

    z_km = get_z_km(nz)

    # How often to save progress plots (every ~10% of training)
    plot_interval = max(1, N_EPOCHS ÷ 10)

    best_val_loss = Inf32
    best_ps = deepcopy(ps)
    epochs_without_improvement = 0
    train_losses = Float32[]
    val_losses = Float32[]

    # Initialise the CSV log
    log_path = joinpath(run_dir, "training_log.csv")
    open(log_path, "w") do io
        println(io, "epoch,lr,train_loss,val_loss,best_val_loss")
    end

    println("\nTraining for $N_EPOCHS epochs, batch_size=$BATCH_SIZE, device=$(USE_GPU ? "GPU" : "CPU")")
    println("Progress plots saved to: $plots_dir")
    println("-" ^ 60)

    n_train_samples = size(X_train_d, 2)
    for epoch in 1:N_EPOCHS
        lr = cosine_lr(epoch - 1, N_EPOCHS)
        Optimisers.adjust!(opt_state; eta = lr)

        shuffle_perm = randperm(rng, n_train_samples)
        n_batches = cld(n_train_samples, BATCH_SIZE)
        epoch_loss = 0.0f0

        for b in 1:n_batches
            i_start = (b - 1) * BATCH_SIZE + 1
            i_end   = min(b * BATCH_SIZE, n_train_samples)
            idx = shuffle_perm[i_start:i_end]
            Xb  = X_train_d[:, idx]
            Yb  = Y_train_d[:, idx]
            dzb = dz_train_d[:, idx]

            loss_val, grads = Zygote.withgradient(ps) do p
                pred, _ = Lux.apply(model, Xb, p, st)
                compute_loss(pred, Yb, dzb, nz)
            end

            epoch_loss += loss_val
            opt_state, ps = Optimisers.update(opt_state, ps, grads[1])
        end

        epoch_loss /= n_batches

        # Validation loss (disable dropout)
        st_test = Lux.testmode(st)
        pred_val, _ = Lux.apply(model, X_val_d, ps, st_test)
        val_loss = compute_loss(pred_val, Y_val_d, dz_val_d, nz)

        push!(train_losses, epoch_loss)
        push!(val_losses, val_loss)

        if val_loss < best_val_loss
            best_val_loss = val_loss
            best_ps = deepcopy(ps)
            epochs_without_improvement = 0
            marker = " *"
        else
            epochs_without_improvement += 1
            marker = ""
        end

        @printf("Epoch %3d/%d  lr=%.2e  train_loss=%.4e  val_loss=%.4e%s\n",
                epoch, N_EPOCHS, lr, epoch_loss, val_loss, marker)
        flush(stdout)

        # Append to CSV log (live-readable)
        open(log_path, "a") do io
            @printf(io, "%d,%.4e,%.6e,%.6e,%.6e\n", epoch, lr, epoch_loss, val_loss, best_val_loss)
        end

        # Save progress plots periodically
        if epoch % plot_interval == 0 || epoch == 1 || epoch == N_EPOCHS
            try
                save_loss_plot(train_losses, val_losses,
                               joinpath(plots_dir, "training_loss.png"))
                save_profile_snapshot(model, ps, st, X_val_d, Y_val_d, dz_val_d,
                                     nz, Y_norm_stats, z_km, epoch,
                                     joinpath(plots_dir, "profiles_epoch_$(@sprintf("%03d", epoch)).png"))
            catch e
                @warn "Failed to save progress plot" exception=e
            end
        end

        if epochs_without_improvement >= PATIENCE
            println("Early stopping: no improvement for $PATIENCE epochs")
            break
        end
    end

    println("-" ^ 60)
    @printf("Best validation loss: %.4e\n", best_val_loss)

    # Final loss plot
    save_loss_plot(train_losses, val_losses, joinpath(plots_dir, "training_loss.png"))

    # Move best params back to CPU for saving
    best_ps_cpu = best_ps |> cpu_device()
    st_cpu = st |> cpu_device()

    save_path = joinpath(run_dir, "flux_correction_model.bson")
    arch = ARCHITECTURE
    temporal_averaging = TEMPORAL_AVERAGING
    predict_mode = PREDICT_MODE
    target_vars = TARGET_VARS
    BSON.@save save_path ps=best_ps_cpu st=st_cpu X_norm_stats Y_norm_stats nz input_dim arch temporal_averaging predict_mode target_vars
    println("Model saved to $save_path")

    # Write a human-readable summary
    open(joinpath(run_dir, "config.txt"), "w") do io
        @printf(io, "architecture  = %s\n", ARCHITECTURE)
        @printf(io, "temporal_avg  = %s\n", TEMPORAL_AVERAGING)
        @printf(io, "predict_mode  = %s\n", PREDICT_MODE)
        @printf(io, "target_vars   = %s\n", TARGET_VARS)
        @printf(io, "hidden_dim    = %d\n", HIDDEN_DIM)
        @printf(io, "learning_rate = %.1e\n", LEARNING_RATE)
        @printf(io, "lr_min        = %.1e\n", LR_MIN)
        @printf(io, "lr_schedule   = cosine_annealing\n")
        @printf(io, "batch_size    = %d\n", BATCH_SIZE)
        @printf(io, "n_epochs      = %d\n", N_EPOCHS)
        @printf(io, "patience      = %d\n", PATIENCE)
        @printf(io, "train_frac    = %.2f\n", TRAIN_FRACTION)
        @printf(io, "z_max         = %.0f\n", Z_MAX)
        @printf(io, "features      = %s\n", join([f[1] for f in INPUT_FEATURES], ", "))
        @printf(io, "col_subsample = %.3f\n", COLUMN_SUBSAMPLE_FRACTION)
        @printf(io, "time_days     = %d-%d\n", TIME_DAY_START, TIME_DAY_END)
        @printf(io, "n_train       = %d\n", n_train)
        @printf(io, "n_val         = %d\n", n_samples - n_train)
        @printf(io, "best_val_loss = %.6e\n", best_val_loss)
        @printf(io, "weight_decay  = %.1e\n", WEIGHT_DECAY)
        @printf(io, "flux_smooth_w = %.1e\n", FLUX_SMOOTH_WEIGHT)
        @printf(io, "device        = %s\n", USE_GPU ? "GPU" : "CPU")
    end

    return model, best_ps, st, run_dir
end

# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

if abspath(PROGRAM_FILE) == @__FILE__
    _, _, _, run_dir = train()
    # Write a pointer so evaluate_flux_model.jl can find the latest run
    open(joinpath(@__DIR__, "latest_run.txt"), "w") do io
        println(io, run_dir)
    end
    println("\nRun directory: $run_dir")
end
