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

# GPU support — set USE_GPU=true via environment variable to enable.
# On SLURM GPU jobs, the sbatch script sets this before launching.
const USE_GPU = get(ENV, "USE_GPU", "false") == "true" && let
    local gpu_avail = false
    try
        @eval using CUDA
        gpu_avail = CUDA.functional()
    catch e
        @warn "CUDA loading failed, falling back to CPU" exception=e
    end
    gpu_avail
end

if USE_GPU
    println("GPU detected: ", CUDA.device())
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
const TRAIN_FRACTION = 0.8
const HIDDEN_DIM     = 128
const PATIENCE       = 15     # early stopping: stop after this many epochs without improvement

const LOSS_WEIGHT_T  = 1.0f0
# LOSS_WEIGHT_Q is computed from data: var(target_T) / var(target_q)
# so both terms contribute equally to the total loss.

const ARCHITECTURE   = :cnn    # :unet, :mlp, or :cnn

# ── Input features ───────────────────────────────────────────────────────────
# Each entry is (netcdf_varname, transform) where transform is applied per-element
# after loading. Add/remove entries here to change what the NN sees.
const INPUT_FEATURES = [
    ("ta",    identity),       # temperature (K) — the variable we're correcting
    ("hus",   identity),       # specific humidity (kg/kg) — the other target variable
    ("pfull", x -> log(max(x, 1f0))),  # log-pressure — encodes altitude / stratification
    ("hur",   identity),       # relative humidity — key for cloud/convection regimes
    ("tke",   identity),       # turbulent kinetic energy — characterises SGS turbulence
]
const N_FEATURES = length(INPUT_FEATURES)

# ── Data subsetting ──────────────────────────────────────────────────────────
# Time window: only use model timesteps within this day range (1-indexed from
# the start of the run). Set to (1, Inf) to use all available timesteps.
const TIME_DAY_START = 3       # skip day 1 (spinup)
const TIME_DAY_END   = 4      # through end of day 8

# Spatial subsampling: randomly keep this fraction of columns per timestep.
# Set to 1.0 to use all columns (no subsampling).
const COLUMN_SUBSAMPLE_FRACTION = 0.02   # 5% of 192×96 ≈ 922 columns per step

const DATASET_CACHE_DIR = joinpath(@__DIR__, "cached_datasets")

# ─────────────────────────────────────────────────────────────────────────────
# Dataset caching — skip the expensive ERA5 regridding pipeline on repeat runs
# ─────────────────────────────────────────────────────────────────────────────

function dataset_cache_key()
    parts = string(
        MODEL_DIR, "|", ERA5_DIR, "|",
        Z_MAX, "|", TIME_DAY_START, "-", TIME_DAY_END, "|",
        COLUMN_SUBSAMPLE_FRACTION, "|",
        join([f[1] for f in INPUT_FEATURES], ","),
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

"""
Build the full training dataset:
  - features X: (n_features × n_levels, n_columns)
  - targets  Y: (2 × n_levels, n_columns)   [T tendency, q tendency]
  - dz:         (n_levels, n_columns)
"""
function build_dataset()
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

        X_step = zeros(Float32, N_FEATURES * nz, nc)
        Y_step = zeros(Float32, 2 * nz, nc)
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

            Y_step[1:nz, col]     = dT_corr[i, j, :]
            Y_step[nz+1:2*nz, col] = dq_corr[i, j, :]
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

# ─────────────────────────────────────────────────────────────────────────────
# Normalisation
# ─────────────────────────────────────────────────────────────────────────────

struct NormStats
    mean::Vector{Float32}
    std::Vector{Float32}
end

function compute_norm(X::Matrix{Float32})
    m = vec(mean(X, dims = 2))
    s = vec(std(X, dims = 2))
    s[s .< 1f-8] .= 1f0
    NormStats(m, s)
end

normalise(X, ns::NormStats)   = (X .- ns.mean) ./ ns.std
denormalise(X, ns::NormStats) = X .* ns.std .+ ns.mean

# ─────────────────────────────────────────────────────────────────────────────
# Neural network architectures
# ─────────────────────────────────────────────────────────────────────────────

"""
Build the flux-correction NN. Dispatches on ARCHITECTURE.

Input:  (n_features × nz, batch)
Output: (2*(nz-1), batch) — interior face fluxes for T and q

Boundary fluxes F[bottom]=0, F[top]=0 are enforced in the loss.
"""
function build_model(input_dim::Int, nz::Int)
    if ARCHITECTURE == :unet
        return build_unet_model(input_dim, nz)
    elseif ARCHITECTURE == :mlp
        return build_mlp_model(input_dim, nz)
    elseif ARCHITECTURE == :cnn
        return build_cnn_model(input_dim, nz)
    else
        error("Unknown ARCHITECTURE: $ARCHITECTURE. Use :unet, :mlp, or :cnn")
    end
end

function build_mlp_model(input_dim::Int, nz::Int)
    n_interior_faces = nz - 1
    output_dim = 2 * n_interior_faces

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
function build_cnn_model(input_dim::Int, nz::Int)
    n_feat = N_FEATURES
    @assert input_dim == n_feat * nz "input_dim ($input_dim) != N_FEATURES ($n_feat) * nz ($nz)"

    nf = nz - 1
    ch = max(HIDDEN_DIM ÷ 4, 16)

    dilations = [1, 1, 2, 4, 1, 1]

    model = @compact(
        conv1 = Conv((3,), n_feat => ch, gelu; pad = 1 * dilations[1], dilation = dilations[1]),
        conv2 = Conv((3,), ch => ch, gelu; pad = 1 * dilations[2], dilation = dilations[2]),
        conv3 = Conv((3,), ch => 2ch, gelu; pad = 1 * dilations[3], dilation = dilations[3]),
        conv4 = Conv((3,), 2ch => 2ch, gelu; pad = 1 * dilations[4], dilation = dilations[4]),
        conv5 = Conv((3,), 2ch => ch, gelu; pad = 1 * dilations[5], dilation = dilations[5]),
        conv6 = Conv((3,), ch => ch, gelu; pad = 1 * dilations[6], dilation = dilations[6]),
        out_conv = Conv((1,), ch => 2),
    ) do x
        B = size(x, 2)

        # (input_dim, B) → (nz, n_feat, B)
        h = reshape(x, nz, n_feat, B)

        h = conv1(h)
        h = conv2(h)
        h = conv3(h)
        h = conv4(h)
        h = conv5(h)
        h = conv6(h)

        out = out_conv(h)  # (nz, 2, B)

        # Interpolate cell centres → interior faces
        faces = (out[1:nf, :, :] .+ out[2:nz, :, :]) ./ 2f0  # (nf, 2, B)

        @return vcat(reshape(faces[:, 1:1, :], nf, B),
                     reshape(faces[:, 2:2, :], nf, B))
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
function build_unet_model(input_dim::Int, nz::Int)
    n_feat = N_FEATURES
    @assert input_dim == n_feat * nz "input_dim ($input_dim) != N_FEATURES ($n_feat) * nz ($nz)"

    nf = nz - 1
    nz_pad = nz % 4 == 0 ? nz : nz + (4 - nz % 4)
    ch = max(HIDDEN_DIM ÷ 4, 16)  # base channels (32 for HIDDEN_DIM=128)

    model = @compact(
        # ── Encoder ──
        enc1 = Chain(Conv((3,), n_feat => ch, gelu; pad = 1),
                     Conv((3,), ch => ch, gelu; pad = 1)),
        enc2 = Chain(Conv((3,), ch => 2ch, gelu; pad = 1),
                     Conv((3,), 2ch => 2ch, gelu; pad = 1)),
        # ── Bottleneck ──
        bneck = Chain(Conv((3,), 2ch => 4ch, gelu; pad = 1),
                      Conv((3,), 4ch => 4ch, gelu; pad = 1)),
        # ── Decoder ──
        up2     = ConvTranspose((2,), 4ch => 2ch; stride = 2),
        dec2    = Chain(Conv((3,), 4ch => 2ch, gelu; pad = 1),
                        Conv((3,), 2ch => 2ch, gelu; pad = 1)),
        up1     = ConvTranspose((2,), 2ch => ch; stride = 2),
        dec1    = Chain(Conv((3,), 2ch => ch, gelu; pad = 1),
                        Conv((3,), ch => ch, gelu; pad = 1)),
        # ── Output head ──
        out_conv = Conv((1,), ch => 2),
        pool     = MaxPool((2,)),
    ) do x
        B = size(x, 2)

        # (input_dim, B) → (nz, n_feat, B)
        x3d = reshape(x, nz, n_feat, B)

        # Pad vertical dim to multiple of 4
        if nz_pad > nz
            x3d = cat(x3d, x3d[1:nz_pad - nz, :, :] .* 0f0; dims = 1)
        end

        # Encoder
        e1 = enc1(x3d)        # (nz_pad,   ch,  B)
        e2 = enc2(pool(e1))   # (nz_pad/2, 2ch, B)
        b  = bneck(pool(e2))  # (nz_pad/4, 4ch, B)

        # Decoder — upsample + cat skip + conv
        d2 = dec2(cat(up2(b), e2; dims = 2))   # (nz_pad/2, 2ch, B)
        d1 = dec1(cat(up1(d2), e1; dims = 2))  # (nz_pad,   ch,  B)

        # 1×1 conv → 2 channels (T, q) at cell centres
        out = out_conv(d1)                      # (nz_pad, 2, B)

        # Crop padding, interpolate centres → interior faces
        c = out[1:nz, :, :]                     # (nz, 2, B)
        faces = (c[1:nf, :, :] .+ c[2:nz, :, :]) ./ 2f0  # (nf, 2, B)

        # Stack as [T_faces; q_faces]  →  (2*nf, B)
        @return vcat(reshape(faces[:, 1:1, :], nf, B),
                     reshape(faces[:, 2:2, :], nf, B))
    end

    return model
end

# ─────────────────────────────────────────────────────────────────────────────
# Loss function: flux divergence must match target tendency
# ─────────────────────────────────────────────────────────────────────────────

"""
Given interior face fluxes (nz-1 per variable), reconstruct full face fluxes
(nz+1) with zero BCs, compute the divergence, and compare to target.

flux_interior: (2*(nz-1), batch)
target:        (2*nz, batch) — target correction tendencies
dz:            (nz, batch)
loss_weight_q: scalar weight for moisture term (computed from data variance ratio)
"""
function flux_divergence_loss(flux_interior, target, dz, nz, loss_weight_q)
    batch = size(flux_interior, 2)
    nf = nz - 1

    F_T_int = flux_interior[1:nf, :]
    F_q_int = flux_interior[nf+1:2*nf, :]

    z_row = flux_interior[1:1, :] .* 0f0
    F_T = vcat(z_row, F_T_int, z_row)
    F_q = vcat(z_row, F_q_int, z_row)

    dFdz_T = -(F_T[2:nz+1, :] .- F_T[1:nz, :]) ./ dz
    dFdz_q = -(F_q[2:nz+1, :] .- F_q[1:nz, :]) ./ dz

    target_T = target[1:nz, :]
    target_q = target[nz+1:2*nz, :]

    loss_T = mean((dFdz_T .- target_T) .^ 2)
    loss_q = mean((dFdz_q .- target_q) .^ 2)

    return LOSS_WEIGHT_T * loss_T + loss_weight_q * loss_q
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
    Y_norm_stats = compute_norm(Y_raw)

    X = normalise(X_raw, X_norm_stats)
    Y = Y_raw

    # Compute loss weight for q from the variance ratio so both terms
    # contribute equally to the total gradient signal.
    var_T = var(Y_raw[1:nz, :])
    var_q = var(Y_raw[nz+1:2*nz, :])
    loss_weight_q = Float32(var_T / var_q)
    @printf("  var(target_T)=%.4e  var(target_q)=%.4e  loss_weight_q=%.2f\n",
            var_T, var_q, loss_weight_q)

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
    opt_state = Optimisers.setup(Adam(LEARNING_RATE), ps)

    # Move ALL data to device once — dataset is small enough to fit in VRAM
    X_train_d  = X_train  |> DEV
    Y_train_d  = Y_train  |> DEV
    dz_train_d = dz_train |> DEV
    X_val_d    = X_val    |> DEV
    Y_val_d    = Y_val    |> DEV
    dz_val_d   = dz_val   |> DEV

    best_val_loss = Inf32
    best_ps = deepcopy(ps)
    epochs_without_improvement = 0

    println("\nTraining for $N_EPOCHS epochs, batch_size=$BATCH_SIZE, device=$(USE_GPU ? "GPU" : "CPU")")
    println("-" ^ 60)

    n_train_samples = size(X_train_d, 2)
    for epoch in 1:N_EPOCHS
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
                flux_pred, _ = Lux.apply(model, Xb, p, st)
                flux_divergence_loss(flux_pred, Yb, dzb, nz, loss_weight_q)
            end

            epoch_loss += loss_val
            opt_state, ps = Optimisers.update(opt_state, ps, grads[1])
        end

        epoch_loss /= n_batches

        # Validation loss
        flux_val, _ = Lux.apply(model, X_val_d, ps, st)
        val_loss = flux_divergence_loss(flux_val, Y_val_d, dz_val_d, nz, loss_weight_q)

        if val_loss < best_val_loss
            best_val_loss = val_loss
            best_ps = deepcopy(ps)
            epochs_without_improvement = 0
            marker = " *"
        else
            epochs_without_improvement += 1
            marker = ""
        end

        @printf("Epoch %3d/%d  train_loss=%.4e  val_loss=%.4e%s\n",
                epoch, N_EPOCHS, epoch_loss, val_loss, marker)
        flush(stdout)

        if epochs_without_improvement >= PATIENCE
            println("Early stopping: no improvement for $PATIENCE epochs")
            break
        end
    end

    println("-" ^ 60)
    @printf("Best validation loss: %.4e\n", best_val_loss)

    # Move best params back to CPU for saving
    best_ps_cpu = best_ps |> cpu_device()
    st_cpu = st |> cpu_device()

    save_path = joinpath(@__DIR__, "flux_correction_model.bson")
    arch = ARCHITECTURE
    BSON.@save save_path ps=best_ps_cpu st=st_cpu X_norm_stats nz input_dim arch
    println("Model saved to $save_path")

    return model, best_ps, st
end

# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

if abspath(PROGRAM_FILE) == @__FILE__
    train()
end
