#!/usr/bin/env julia
#=

julia -t 24 --project=. sanity_check_compute_tendencies.jl

Sanity checks for the offline *computed tendencies* used in flux-correction
training (ERA5 − model hourly finite differences).

Reuses the same loaders and time/level window as `train_flux_correction.jl`
(MODEL_DIR, ERA5_DIR, Z_MAX, TIME_DAY_START/END, DT, TEMPORAL_AVERAGING).

Supports both :hourly and :daily modes (reads TEMPORAL_AVERAGING from the
training config).  In :daily mode the hourly corrections within each calendar
day are averaged before accumulating — matching what the NN actually sees.

Checks
  1) Telescoping identity (hourly mode only): sum over valid hourly pairs of
     dT_corr*DT  ≈  (T_ERA5[last]−T_ERA5[first]) − (T_model[last]−T_model[first])

Plots (written to output dir, default `ml_flux/sanity_tendency_plots/`):
  01  Mean vertical profiles of dT/dt from ERA5, model, and correction (q too)
  02  Telescoping identity RMS vs height (T) — hourly mode only
  03  Telescoping identity RMS vs height (q) — hourly mode only
  04  Sample column profiles
  05  Spatial dT_corr at selected levels (single snapshot)
  06  Spatial dq_corr at selected levels
  07  Zonal-mean latitude-height cross-section of correction
  08  RMS vs |mean| vertical profiles
  09  Temporal evolution of global-mean correction
  10  Histograms of dT_corr at representative levels

Usage:
    julia -t 24 --project=. sanity_check_compute_tendencies.jl
    julia -t 24 --project=. sanity_check_compute_tendencies.jl runs/007_cnn_h128_day
=#

include("train_flux_correction.jl")

using CairoMakie
using Printf
using Random

function time_indices(model_dates, nt)
    if DATA_SOURCE == :daily_era5
        ti_start = max(1, TIME_DAY_START)
        ti_end   = min(nt, TIME_DAY_END == Inf ? nt : Int(TIME_DAY_END))
        ti_end   = min(ti_end, nt) - 1
        return ti_start, ti_end
    end
    t0 = model_dates[1]
    ti_start = max(
        1,
        something(
            findfirst(d -> (d - t0).value / 3600_000 >= 24 * (TIME_DAY_START - 1), model_dates),
            1,
        ),
    )
    ti_end = something(
        findlast(d -> (d - t0).value / 3600_000 < 24 * TIME_DAY_END, model_dates),
        nt,
    )
    ti_end = min(ti_end, nt) - 1
    return ti_start, ti_end
end

"""Pick model level indices for spatial maps: near-surface, lower/mid/upper troposphere."""
function level_indices_for_maps(nz::Int)
    ks = Int[1, max(1, nz ÷ 4), max(1, nz ÷ 2), max(1, 3 * nz ÷ 4), nz]
    unique!(sort!(ks))
    return ks
end

function regrid_era5_to_model_z(era5, grid, z_idx)
    t_rg = regrid_horizontal(era5.t, era5.lon, era5.lat, grid.lon, grid.lat)
    q_rg = regrid_horizontal(era5.q, era5.lon, era5.lat, grid.lon, grid.lat)
    z_rg = regrid_horizontal(era5.z, era5.lon, era5.lat, grid.lon, grid.lat)
    t_on_z = interp_vertical_3d(t_rg, z_rg, grid.z_phys)[:, :, z_idx]
    q_on_z = interp_vertical_3d(q_rg, z_rg, grid.z_phys)[:, :, z_idx]
    return t_on_z, q_on_z
end

# ─────────────────────────────────────────────────────────────────────────────
# Result container — populated by both hourly and daily loops
# ─────────────────────────────────────────────────────────────────────────────
struct SanityResults
    z_km::Vector{Float32}
    nlon::Int
    nlat::Int
    nz::Int
    lon::Vector{Float64}
    lat::Vector{Float64}
    mean_dT_era5::Vector{Float32}
    mean_dT_model::Vector{Float32}
    mean_dT_corr::Vector{Float32}
    mean_dq_era5::Vector{Float32}
    mean_dq_model::Vector{Float32}
    mean_dq_corr::Vector{Float32}
    # Telescoping (hourly only — empty vectors otherwise)
    rms_diff_T_z::Vector{Float64}
    rms_diff_q_z::Vector{Float64}
    max_abs_diff_T::Float64
    max_abs_diff_q::Float64
    n_valid::Int
    n_skipped::Int
    has_telescoping::Bool
    # Zonal mean, RMS, samples, snapshots, temporal, histogram
    acc_dT_corr_zonal::Matrix{Float64}
    acc_dq_corr_zonal::Matrix{Float64}
    acc_dT_corr_sq::Vector{Float64}
    acc_dq_corr_sq::Vector{Float64}
    n_prof::Int
    sample_dT_era5::Vector{Vector{Float32}}
    sample_dT_mod::Vector{Vector{Float32}}
    sample_dT_corr::Vector{Vector{Float32}}
    sample_meta::Vector{String}
    dT_corr_snapshot::Array{Float32, 3}
    dq_corr_snapshot::Array{Float32, 3}
    snapshot_label::String
    ts_dT_corr_mean::Vector{Float64}
    ts_dq_corr_mean::Vector{Float64}
    ts_dates::Vector{DateTime}
    hist_levels::Vector{Int}
    hist_dT_data::Dict{Int, Vector{Float32}}
    mode_label::String   # "hourly" or "daily-mean"
end

# ─────────────────────────────────────────────────────────────────────────────
# Hourly processing loop
# ─────────────────────────────────────────────────────────────────────────────
function process_hourly(grid, ta_model, hus_model, z_idx, nz, ti_start, ti_end, model_dates)
    nlon, nlat = length(grid.lon), length(grid.lat)
    rng = MersenneTwister(42)

    acc_dT_era5 = zeros(Float64, nz);  acc_dT_model = zeros(Float64, nz);  acc_dT_corr = zeros(Float64, nz)
    acc_dq_era5 = zeros(Float64, nz);  acc_dq_model = zeros(Float64, nz);  acc_dq_corr = zeros(Float64, nz)
    n_prof = 0

    sum_dT_corr_dt = zeros(Float64, nlon, nlat, nz)
    sum_dq_corr_dt = zeros(Float64, nlon, nlat, nz)
    n_skipped = 0;  n_valid = 0

    era5_t0 = zeros(Float32, nlon, nlat, nz); era5_t1 = similar(era5_t0)
    q_era5_t0 = similar(era5_t0);  q_era5_t1 = similar(era5_t0)
    model_t0 = similar(era5_t0);   model_t1 = similar(era5_t0)
    hus_t0 = similar(era5_t0);     hus_t1 = similar(era5_t0)
    chain_started = false

    acc_dT_corr_zonal = zeros(Float64, nlat, nz)
    acc_dq_corr_zonal = zeros(Float64, nlat, nz)
    acc_dT_corr_sq = zeros(Float64, nz)
    acc_dq_corr_sq = zeros(Float64, nz)

    ts_dT = Float64[];  ts_dq = Float64[];  ts_dt = DateTime[]

    hist_levels = level_indices_for_maps(nz)
    hist_dT = Dict{Int, Vector{Float32}}(k => Float32[] for k in hist_levels)

    n_samples_wanted = 12
    s_era5 = Vector{Vector{Float32}}();  s_mod = Vector{Vector{Float32}}()
    s_corr = Vector{Vector{Float32}}();  s_meta = Vector{String}()

    ti_mid = (ti_start + ti_end) ÷ 2
    dT_snap = zeros(Float32, nlon, nlat, nz)
    dq_snap = zeros(Float32, nlon, nlat, nz)
    snap_label = ""
    best_snap_dist = typemax(Int)

    prev_dt_next = nothing
    cache_t = Array{Float32}(undef, 0, 0, 0)
    cache_q = Array{Float32}(undef, 0, 0, 0)

    for ti in ti_start:ti_end
        dt_now = model_dates[ti];  dt_next = model_dates[ti + 1]
        era5_now = load_era5_pressure(dt_now)
        era5_next = load_era5_pressure(dt_next)
        if era5_now === nothing || era5_next === nothing
            n_skipped += 1;  prev_dt_next = nothing;  continue
        end
        n_valid += 1

        if ti == ti_start || ti % 24 == 0
            @printf("  [hourly] ti %d / %d  %s\n", ti, ti_end, dt_now);  flush(stdout)
        end

        if prev_dt_next == dt_now
            t_now_z, q_now_z = cache_t, cache_q
        else
            t_now_z, q_now_z = regrid_era5_to_model_z(era5_now, grid, z_idx)
        end
        t_next_z, q_next_z = regrid_era5_to_model_z(era5_next, grid, z_idx)
        cache_t, cache_q, prev_dt_next = t_next_z, q_next_z, dt_next

        ta_now_z  = ta_model[ti, :, :, z_idx];     ta_next_z  = ta_model[ti+1, :, :, z_idx]
        hus_now_z = hus_model[ti, :, :, z_idx];     hus_next_z = hus_model[ti+1, :, :, z_idx]

        dT_era5  = (t_next_z  .- t_now_z)  ./ DT
        dq_era5  = (q_next_z  .- q_now_z)  ./ DT
        dT_model = (ta_next_z .- ta_now_z) ./ DT
        dq_model = (hus_next_z .- hus_now_z) ./ DT
        dT_corr  = dT_era5 .- dT_model
        dq_corr  = dq_era5 .- dq_model

        n_prof += nlon * nlat
        for k in 1:nz
            acc_dT_era5[k]  += sum(Float64, @view dT_era5[:, :, k])
            acc_dT_model[k] += sum(Float64, @view dT_model[:, :, k])
            acc_dT_corr[k]  += sum(Float64, @view dT_corr[:, :, k])
            acc_dq_era5[k]  += sum(Float64, @view dq_era5[:, :, k])
            acc_dq_model[k] += sum(Float64, @view dq_model[:, :, k])
            acc_dq_corr[k]  += sum(Float64, @view dq_corr[:, :, k])
        end
        for k in 1:nz, j in 1:nlat
            acc_dT_corr_zonal[j, k] += sum(Float64, @view dT_corr[:, j, k])
            acc_dq_corr_zonal[j, k] += sum(Float64, @view dq_corr[:, j, k])
        end
        for k in 1:nz
            acc_dT_corr_sq[k] += mapreduce(x -> Float64(x)^2, +, @view dT_corr[:, :, k])
            acc_dq_corr_sq[k] += mapreduce(x -> Float64(x)^2, +, @view dq_corr[:, :, k])
        end

        push!(ts_dT, Float64(mean(dT_corr)));  push!(ts_dq, Float64(mean(dq_corr)));  push!(ts_dt, dt_now)

        for k in hist_levels
            slab = @view dT_corr[:, :, k]
            stride = max(1, length(slab) ÷ 500)
            append!(hist_dT[k], slab[1:stride:end])
        end

        sum_dT_corr_dt .+= Float64.(dT_corr) .* Float64(DT)
        sum_dq_corr_dt .+= Float64.(dq_corr) .* Float64(DT)

        if !chain_started
            era5_t0 .= t_now_z;  q_era5_t0 .= q_now_z
            model_t0 .= ta_now_z;  hus_t0 .= hus_now_z
            chain_started = true
        end
        era5_t1 .= t_next_z;  q_era5_t1 .= q_next_z
        model_t1 .= ta_next_z;  hus_t1 .= hus_next_z

        d = abs(ti - ti_mid)
        if d < best_snap_dist
            best_snap_dist = d
            dT_snap .= dT_corr;  dq_snap .= dq_corr
            snap_label = "hourly ti=$ti  $dt_now"
        end

        if length(s_corr) < n_samples_wanted
            for _ in 1:min(2, n_samples_wanted - length(s_corr))
                i, j = rand(rng, 1:nlon), rand(rng, 1:nlat)
                push!(s_era5, collect(Float32, dT_era5[i, j, :]))
                push!(s_mod,  collect(Float32, dT_model[i, j, :]))
                push!(s_corr, collect(Float32, dT_corr[i, j, :]))
                push!(s_meta, @sprintf("ti=%d i=%d j=%d %s", ti, i, j, dt_now))
            end
        end
    end

    chain_started || error("No valid ERA5 pairs — check ERA5_DIR and dates.")
    z_km = Float32.(grid.z_ref[z_idx] ./ 1000.0)

    inv_n = 1.0 / max(n_prof, 1)
    mn(a) = Float32.(a .* inv_n)

    direct_dT = era5_t1 .- era5_t0 .- (model_t1 .- model_t0)
    direct_dq = q_era5_t1 .- q_era5_t0 .- (hus_t1 .- hus_t0)
    diff_T = sum_dT_corr_dt .- Float64.(direct_dT)
    diff_q = sum_dq_corr_dt .- Float64.(direct_dq)
    rms_T = [sqrt(mean(diff_T[:, :, k] .^ 2)) for k in 1:nz]
    rms_q = [sqrt(mean(diff_q[:, :, k] .^ 2)) for k in 1:nz]

    return SanityResults(
        z_km, nlon, nlat, nz, grid.lon, grid.lat,
        mn(acc_dT_era5), mn(acc_dT_model), mn(acc_dT_corr),
        mn(acc_dq_era5), mn(acc_dq_model), mn(acc_dq_corr),
        rms_T, rms_q, maximum(abs, diff_T), maximum(abs, diff_q),
        n_valid, n_skipped, true,
        acc_dT_corr_zonal, acc_dq_corr_zonal,
        acc_dT_corr_sq, acc_dq_corr_sq, n_prof,
        s_era5, s_mod, s_corr, s_meta,
        dT_snap, dq_snap, snap_label,
        ts_dT, ts_dq, ts_dt,
        hist_levels, hist_dT,
        "hourly",
    )
end

# ─────────────────────────────────────────────────────────────────────────────
# Daily-mean processing loop
# ─────────────────────────────────────────────────────────────────────────────
function process_daily(grid, ta_model, hus_model, z_idx, nz, ti_start, ti_end, model_dates)
    nlon, nlat = length(grid.lon), length(grid.lat)
    rng = MersenneTwister(42)

    # Group hourly timestep indices by calendar day
    day_groups = Dict{Date, Vector{Int}}()
    for ti in ti_start:ti_end
        d = Date(model_dates[ti])
        push!(get!(day_groups, d, Int[]), ti)
    end
    days_sorted = sort(collect(keys(day_groups)))
    println("  Daily-mean mode: $(length(days_sorted)) days ($(days_sorted[1]) to $(days_sorted[end]))")

    acc_dT_era5 = zeros(Float64, nz);  acc_dT_model = zeros(Float64, nz);  acc_dT_corr = zeros(Float64, nz)
    acc_dq_era5 = zeros(Float64, nz);  acc_dq_model = zeros(Float64, nz);  acc_dq_corr = zeros(Float64, nz)
    n_prof = 0;  n_valid_days = 0;  n_skipped_hours = 0

    acc_dT_corr_zonal = zeros(Float64, nlat, nz)
    acc_dq_corr_zonal = zeros(Float64, nlat, nz)
    acc_dT_corr_sq = zeros(Float64, nz)
    acc_dq_corr_sq = zeros(Float64, nz)

    ts_dT = Float64[];  ts_dq = Float64[];  ts_dt = DateTime[]

    hist_levels = level_indices_for_maps(nz)
    hist_dT = Dict{Int, Vector{Float32}}(k => Float32[] for k in hist_levels)

    n_samples_wanted = 12
    s_era5 = Vector{Vector{Float32}}();  s_mod = Vector{Vector{Float32}}()
    s_corr = Vector{Vector{Float32}}();  s_meta = Vector{String}()

    mid_day_idx = max(1, length(days_sorted) ÷ 2)
    dT_snap = zeros(Float32, nlon, nlat, nz)
    dq_snap = zeros(Float32, nlon, nlat, nz)
    snap_label = ""

    for (di, day) in enumerate(days_sorted)
        hours = day_groups[day]
        @printf("  [daily] day %s  (%d hourly pairs)\n", day, length(hours))
        flush(stdout)

        dT_era5_acc = zeros(Float64, nlon, nlat, nz)
        dT_mod_acc  = zeros(Float64, nlon, nlat, nz)
        dq_era5_acc = zeros(Float64, nlon, nlat, nz)
        dq_mod_acc  = zeros(Float64, nlon, nlat, nz)
        n_hrs = 0

        prev_dt_next_d = nothing
        cache_t_d = Array{Float32}(undef, 0, 0, 0)
        cache_q_d = Array{Float32}(undef, 0, 0, 0)

        for ti in hours
            dt_now = model_dates[ti];  dt_next = model_dates[ti + 1]
            era5_now = load_era5_pressure(dt_now)
            era5_next = load_era5_pressure(dt_next)
            if era5_now === nothing || era5_next === nothing
                n_skipped_hours += 1;  prev_dt_next_d = nothing;  continue
            end

            if prev_dt_next_d == dt_now
                t_now_z, q_now_z = cache_t_d, cache_q_d
            else
                t_now_z, q_now_z = regrid_era5_to_model_z(era5_now, grid, z_idx)
            end
            t_next_z, q_next_z = regrid_era5_to_model_z(era5_next, grid, z_idx)
            cache_t_d, cache_q_d, prev_dt_next_d = t_next_z, q_next_z, dt_next

            ta_now_z  = ta_model[ti, :, :, z_idx];     ta_next_z  = ta_model[ti+1, :, :, z_idx]
            hus_now_z = hus_model[ti, :, :, z_idx];     hus_next_z = hus_model[ti+1, :, :, z_idx]

            dT_era5_acc .+= Float64.((t_next_z .- t_now_z) ./ DT)
            dq_era5_acc .+= Float64.((q_next_z .- q_now_z) ./ DT)
            dT_mod_acc  .+= Float64.((ta_next_z .- ta_now_z) ./ DT)
            dq_mod_acc  .+= Float64.((hus_next_z .- hus_now_z) ./ DT)
            n_hrs += 1
        end

        n_hrs == 0 && continue
        n_valid_days += 1
        inv_h = 1.0 / n_hrs
        dT_era5_day = Float32.(dT_era5_acc .* inv_h)
        dT_mod_day  = Float32.(dT_mod_acc .* inv_h)
        dq_era5_day = Float32.(dq_era5_acc .* inv_h)
        dq_mod_day  = Float32.(dq_mod_acc .* inv_h)
        dT_corr_day = dT_era5_day .- dT_mod_day
        dq_corr_day = dq_era5_day .- dq_mod_day

        n_prof += nlon * nlat
        for k in 1:nz
            acc_dT_era5[k]  += sum(Float64, @view dT_era5_day[:, :, k])
            acc_dT_model[k] += sum(Float64, @view dT_mod_day[:, :, k])
            acc_dT_corr[k]  += sum(Float64, @view dT_corr_day[:, :, k])
            acc_dq_era5[k]  += sum(Float64, @view dq_era5_day[:, :, k])
            acc_dq_model[k] += sum(Float64, @view dq_mod_day[:, :, k])
            acc_dq_corr[k]  += sum(Float64, @view dq_corr_day[:, :, k])
        end
        for k in 1:nz, j in 1:nlat
            acc_dT_corr_zonal[j, k] += sum(Float64, @view dT_corr_day[:, j, k])
            acc_dq_corr_zonal[j, k] += sum(Float64, @view dq_corr_day[:, j, k])
        end
        for k in 1:nz
            acc_dT_corr_sq[k] += mapreduce(x -> Float64(x)^2, +, @view dT_corr_day[:, :, k])
            acc_dq_corr_sq[k] += mapreduce(x -> Float64(x)^2, +, @view dq_corr_day[:, :, k])
        end

        push!(ts_dT, Float64(mean(dT_corr_day)));  push!(ts_dq, Float64(mean(dq_corr_day)))
        push!(ts_dt, DateTime(day))

        for k in hist_levels
            slab = @view dT_corr_day[:, :, k]
            stride = max(1, length(slab) ÷ 500)
            append!(hist_dT[k], slab[1:stride:end])
        end

        if di == mid_day_idx
            dT_snap .= dT_corr_day;  dq_snap .= dq_corr_day
            snap_label = "daily-mean $day ($n_hrs hrs)"
        end

        if length(s_corr) < n_samples_wanted
            for _ in 1:min(2, n_samples_wanted - length(s_corr))
                i, j = rand(rng, 1:nlon), rand(rng, 1:nlat)
                push!(s_era5, collect(Float32, dT_era5_day[i, j, :]))
                push!(s_mod,  collect(Float32, dT_mod_day[i, j, :]))
                push!(s_corr, collect(Float32, dT_corr_day[i, j, :]))
                push!(s_meta, @sprintf("day=%s i=%d j=%d", day, i, j))
            end
        end
    end

    n_valid_days > 0 || error("No valid days — check ERA5_DIR and dates.")
    z_km = Float32.(grid.z_ref[z_idx] ./ 1000.0)
    inv_n = 1.0 / max(n_prof, 1)
    mn(a) = Float32.(a .* inv_n)

    return SanityResults(
        z_km, nlon, nlat, nz, grid.lon, grid.lat,
        mn(acc_dT_era5), mn(acc_dT_model), mn(acc_dT_corr),
        mn(acc_dq_era5), mn(acc_dq_model), mn(acc_dq_corr),
        Float64[], Float64[], 0.0, 0.0,
        n_valid_days, n_skipped_hours, false,
        acc_dT_corr_zonal, acc_dq_corr_zonal,
        acc_dT_corr_sq, acc_dq_corr_sq, n_prof,
        s_era5, s_mod, s_corr, s_meta,
        dT_snap, dq_snap, snap_label,
        ts_dT, ts_dq, ts_dt,
        hist_levels, hist_dT,
        "daily-mean",
    )
end

# ─────────────────────────────────────────────────────────────────────────────
# Shared plotting
# ─────────────────────────────────────────────────────────────────────────────
function make_plots(r::SanityResults, out_dir::String)
    mkpath(out_dir)
    z_km = r.z_km

    # ── 01: mean vertical profiles ──────────────────────────────────────
    fig = Figure(size = (1000, 500))
    ax1 = Axis(fig[1, 1]; xlabel = "dT/dt (K/hr)", ylabel = "Height (km)",
               title = "Mean T tendencies ($(r.mode_label))")
    lines!(ax1, r.mean_dT_era5 .* 3600, z_km; label = "ERA5", linewidth = 2)
    lines!(ax1, r.mean_dT_model .* 3600, z_km; label = "Model", linewidth = 2)
    lines!(ax1, r.mean_dT_corr .* 3600, z_km; label = "Correction", linewidth = 2, linestyle = :dash)
    axislegend(ax1; position = :rt)
    ax2 = Axis(fig[1, 2]; xlabel = "dq/dt (g/kg/hr)", ylabel = "Height (km)",
               title = "Mean q tendencies ($(r.mode_label))")
    lines!(ax2, r.mean_dq_era5 .* 3600 .* 1000, z_km; label = "ERA5", linewidth = 2)
    lines!(ax2, r.mean_dq_model .* 3600 .* 1000, z_km; label = "Model", linewidth = 2)
    lines!(ax2, r.mean_dq_corr .* 3600 .* 1000, z_km; label = "Correction", linewidth = 2, linestyle = :dash)
    axislegend(ax2; position = :rt)
    save(joinpath(out_dir, "01_mean_tendency_profiles.png"), fig; px_per_unit = 2)
    println("Wrote 01_mean_tendency_profiles.png")

    # ── 02/03: telescoping (hourly only) ─────────────────────────────────
    if r.has_telescoping
        for (idx, rms_z, unit, var, color) in [
            ("02", r.rms_diff_T_z, "K", "T", :darkred),
            ("03", r.rms_diff_q_z, "kg/kg", "q", :steelblue)]
            f = Figure(size = (700, 500))
            ax = Axis(f[1, 1]; xlabel = "RMS column error ($unit)", ylabel = "Height (km)",
                       title = "Telescoping mismatch ($var) vs z")
            lines!(ax, rms_z, z_km; linewidth = 2, color = color)
            save(joinpath(out_dir, "$(idx)_integral_identity_rms_$(var).png"), f; px_per_unit = 2)
            println("Wrote $(idx)_integral_identity_rms_$(var).png")
        end
    end

    # ── 04: sample column profiles ─────────────────────────────────────
    n_s = length(r.sample_dT_corr)
    if n_s > 0
        ncols = min(4, max(1, n_s));  nrows = cld(n_s, ncols)
        f = Figure(size = (280 * ncols, 320 * nrows))
        for p in 1:n_s
            rc, cc = cld(p, ncols), mod1(p, ncols)
            ax = Axis(f[rc, cc]; xlabel = "dT/dt (K/hr)",
                       ylabel = cc == 1 ? "z (km)" : "", title = r.sample_meta[p])
            lines!(ax, r.sample_dT_era5[p] .* 3600, z_km; label = "ERA5", linewidth = 1.5)
            lines!(ax, r.sample_dT_mod[p] .* 3600, z_km; label = "Model", linewidth = 1.5)
            lines!(ax, r.sample_dT_corr[p] .* 3600, z_km; label = "Corr", linewidth = 1.5, linestyle = :dash)
            p == 1 && axislegend(ax; position = :rt, labelsize = 8)
        end
        Label(f[0, :], "Sample dT/dt profiles ($(r.mode_label))"; fontsize = 14)
        save(joinpath(out_dir, "04_sample_column_profiles.png"), f; px_per_unit = 2)
        println("Wrote 04_sample_column_profiles.png")
    end

    # ── 05/06: spatial maps (per-level colorscale) ─────────────────────
    ks = level_indices_for_maps(r.nz);  nmaps = length(ks)
    for (pnum, snap, unit, scale, label) in [
        ("05", r.dT_corr_snapshot, "K/hr", 3600f0, "dT_corr"),
        ("06", r.dq_corr_snapshot, "g/kg/hr", 3600f0 * 1000f0, "dq_corr")]
        f = Figure(size = (320 * nmaps, 320))
        Label(f[0, :], "$label ($unit)  $(r.snapshot_label)"; fontsize = 12)
        slab = snap .* scale
        for (ip, k) in enumerate(ks)
            ax = Axis(f[1, ip]; aspect = DataAspect(), xlabel = "lon",
                       ylabel = ip == 1 ? "lat" : "",
                       title = @sprintf("z ≈ %.1f km", z_km[k]))
            ld = slab[:, :, k]
            cl = Float64(maximum(abs, ld; init = 1f-6))
            hm = heatmap!(ax, r.lon, r.lat, ld; colormap = :balance, colorrange = (-cl, cl))
            Colorbar(f[2, ip], hm; label = unit, vertical = false, flipaxis = false)
        end
        save(joinpath(out_dir, "$(pnum)_spatial_$(label)_levels.png"), f; px_per_unit = 2)
        println("Wrote $(pnum)_spatial_$(label)_levels.png")
    end

    # ── 07: zonal-mean lat-height ──────────────────────────────────────
    nd = max(r.n_prof ÷ r.nlat, 1)
    zT = Float32.(r.acc_dT_corr_zonal ./ nd) .* 3600
    zq = Float32.(r.acc_dq_corr_zonal ./ nd) .* 3600 .* 1000
    f = Figure(size = (1000, 450))
    for (col, zd, unit, ttl) in [(1, zT, "K/hr", "T"), (3, zq, "g/kg/hr", "q")]
        ax = Axis(f[1, col]; xlabel = "Latitude", ylabel = "Height (km)",
                   title = "Zonal-mean d$(ttl)_corr ($unit, $(r.mode_label))")
        cl = Float64(maximum(abs, zd; init = 1f-6))
        hm = heatmap!(ax, Float32.(r.lat), z_km, zd; colormap = :balance, colorrange = (-cl, cl))
        Colorbar(f[1, col+1], hm; label = unit)
    end
    save(joinpath(out_dir, "07_zonal_mean_correction.png"), f; px_per_unit = 2)
    println("Wrote 07_zonal_mean_correction.png")

    # ── 08: RMS vs |mean| profiles ────────────────────────────────────
    inv_np = 1.0 / max(r.n_prof, 1)
    rms_dT = Float32.(sqrt.(r.acc_dT_corr_sq .* inv_np)) .* 3600
    rms_dq = Float32.(sqrt.(r.acc_dq_corr_sq .* inv_np)) .* 3600 .* 1000
    abs_mT = abs.(r.mean_dT_corr) .* 3600;  abs_mq = abs.(r.mean_dq_corr) .* 3600 .* 1000
    f = Figure(size = (1000, 500))
    for (col, rms, abm, unit, ttl, clr) in [
        (1, rms_dT, abs_mT, "K/hr", "T", :firebrick),
        (2, rms_dq, abs_mq, "g/kg/hr", "q", :steelblue)]
        ax = Axis(f[1, col]; xlabel = unit, ylabel = "Height (km)",
                   title = "$ttl correction: RMS vs |mean|")
        lines!(ax, rms, z_km; label = "RMS", linewidth = 2, color = clr)
        lines!(ax, abm, z_km; label = "|Mean|", linewidth = 2, color = clr, linestyle = :dash)
        axislegend(ax; position = :rt)
    end
    save(joinpath(out_dir, "08_rms_vs_mean_profiles.png"), f; px_per_unit = 2)
    println("Wrote 08_rms_vs_mean_profiles.png")

    # ── 09: temporal evolution ─────────────────────────────────────────
    n_ts = length(r.ts_dates)
    if n_ts > 1
        x_label = r.mode_label == "daily-mean" ? "Day index" : "Hours since start"
        xs = r.mode_label == "daily-mean" ?
             Float64.(1:n_ts) :
             Float64[(r.ts_dates[i] - r.ts_dates[1]).value / 3600_000 for i in 1:n_ts]
        f = Figure(size = (1000, 450))
        ax1 = Axis(f[1, 1]; xlabel = x_label, ylabel = "K/hr",
                    title = "Global-mean dT_corr vs time ($(r.mode_label))")
        lines!(ax1, xs, r.ts_dT_corr_mean .* 3600; linewidth = 1.5, color = :firebrick)
        hlines!(ax1, [0]; color = :gray, linestyle = :dot)
        ax2 = Axis(f[1, 2]; xlabel = x_label, ylabel = "g/kg/hr",
                    title = "Global-mean dq_corr vs time")
        lines!(ax2, xs, r.ts_dq_corr_mean .* 3600 .* 1000; linewidth = 1.5, color = :steelblue)
        hlines!(ax2, [0]; color = :gray, linestyle = :dot)
        save(joinpath(out_dir, "09_temporal_evolution.png"), f; px_per_unit = 2)
        println("Wrote 09_temporal_evolution.png")
    end

    # ── 10: histograms ────────────────────────────────────────────────
    nhist = length(r.hist_levels)
    f = Figure(size = (300 * min(nhist, 5), 350))
    Label(f[0, :], "Distribution of dT_corr (K/hr, $(r.mode_label))"; fontsize = 13)
    for (ip, k) in enumerate(r.hist_levels)
        vals = r.hist_dT_data[k] .* 3600
        ax = Axis(f[1, ip]; xlabel = "K/hr",
                   ylabel = ip == 1 ? "Count" : "",
                   title = @sprintf("z ≈ %.1f km", z_km[k]))
        hist!(ax, vals; bins = 80, color = (:firebrick, 0.6))
        vlines!(ax, [0]; color = :black, linestyle = :dot)
    end
    save(joinpath(out_dir, "10_histogram_dT_corr.png"), f; px_per_unit = 2)
    println("Wrote 10_histogram_dT_corr.png")
end

# ─────────────────────────────────────────────────────────────────────────────
# Daily-ERA5 processing loop (pre-computed daily-mean ERA5 + model 1d_average)
# ─────────────────────────────────────────────────────────────────────────────
function process_daily_era5(grid, ta_model, hus_model, z_idx, nz, ti_start, ti_end, model_dates)
    nlon, nlat = length(grid.lon), length(grid.lat)
    rng = MersenneTwister(42)

    start_date = Date(model_dates[1])

    acc_dT_era5 = zeros(Float64, nz);  acc_dT_model = zeros(Float64, nz);  acc_dT_corr = zeros(Float64, nz)
    acc_dq_era5 = zeros(Float64, nz);  acc_dq_model = zeros(Float64, nz);  acc_dq_corr = zeros(Float64, nz)
    n_prof = 0;  n_valid = 0;  n_skipped = 0

    acc_dT_corr_zonal = zeros(Float64, nlat, nz)
    acc_dq_corr_zonal = zeros(Float64, nlat, nz)
    acc_dT_corr_sq = zeros(Float64, nz)
    acc_dq_corr_sq = zeros(Float64, nz)

    ts_dT = Float64[];  ts_dq = Float64[];  ts_dt = DateTime[]

    hist_levels = level_indices_for_maps(nz)
    hist_dT = Dict{Int, Vector{Float32}}(k => Float32[] for k in hist_levels)

    n_samples_wanted = 12
    s_era5 = Vector{Vector{Float32}}();  s_mod = Vector{Vector{Float32}}()
    s_corr = Vector{Vector{Float32}}();  s_meta = Vector{String}()

    ti_mid = (ti_start + ti_end) ÷ 2
    dT_snap = zeros(Float32, nlon, nlat, nz)
    dq_snap = zeros(Float32, nlon, nlat, nz)
    snap_label = ""

    era5_cache = Dict{Date, Any}()

    for ti in ti_start:ti_end
        date_now  = start_date + Day(ti - 1)
        date_next = start_date + Day(ti)

        era5_now  = get!(era5_cache, date_now)  do; load_era5_daily(date_now);  end
        era5_next = get!(era5_cache, date_next) do; load_era5_daily(date_next); end

        if era5_now === nothing || era5_next === nothing
            n_skipped += 1;  continue
        end
        n_valid += 1

        @printf("  [daily_era5] pair %d: %s → %s\n", ti, date_now, date_next);  flush(stdout)

        t_now_z, q_now_z = regrid_era5_to_model_z(era5_now, grid, z_idx)
        t_next_z, q_next_z = regrid_era5_to_model_z(era5_next, grid, z_idx)

        ta_now_z  = ta_model[ti, :, :, z_idx];     ta_next_z  = ta_model[ti+1, :, :, z_idx]
        hus_now_z = hus_model[ti, :, :, z_idx];     hus_next_z = hus_model[ti+1, :, :, z_idx]

        dT_era5  = (t_next_z  .- t_now_z)  ./ DT
        dq_era5  = (q_next_z  .- q_now_z)  ./ DT
        dT_model = (ta_next_z .- ta_now_z) ./ DT
        dq_model = (hus_next_z .- hus_now_z) ./ DT
        dT_corr  = dT_era5 .- dT_model
        dq_corr  = dq_era5 .- dq_model

        n_prof += nlon * nlat
        for k in 1:nz
            acc_dT_era5[k]  += sum(Float64, @view dT_era5[:, :, k])
            acc_dT_model[k] += sum(Float64, @view dT_model[:, :, k])
            acc_dT_corr[k]  += sum(Float64, @view dT_corr[:, :, k])
            acc_dq_era5[k]  += sum(Float64, @view dq_era5[:, :, k])
            acc_dq_model[k] += sum(Float64, @view dq_model[:, :, k])
            acc_dq_corr[k]  += sum(Float64, @view dq_corr[:, :, k])
        end
        for k in 1:nz, j in 1:nlat
            acc_dT_corr_zonal[j, k] += sum(Float64, @view dT_corr[:, j, k])
            acc_dq_corr_zonal[j, k] += sum(Float64, @view dq_corr[:, j, k])
        end
        for k in 1:nz
            acc_dT_corr_sq[k] += mapreduce(x -> Float64(x)^2, +, @view dT_corr[:, :, k])
            acc_dq_corr_sq[k] += mapreduce(x -> Float64(x)^2, +, @view dq_corr[:, :, k])
        end

        push!(ts_dT, Float64(mean(dT_corr)));  push!(ts_dq, Float64(mean(dq_corr)))
        push!(ts_dt, DateTime(date_now))

        for k in hist_levels
            slab = @view dT_corr[:, :, k]
            stride = max(1, length(slab) ÷ 500)
            append!(hist_dT[k], slab[1:stride:end])
        end

        if ti == ti_mid
            dT_snap .= dT_corr;  dq_snap .= dq_corr
            snap_label = "daily-era5 $date_now → $date_next"
        end

        if length(s_corr) < n_samples_wanted
            for _ in 1:min(2, n_samples_wanted - length(s_corr))
                i, j = rand(rng, 1:nlon), rand(rng, 1:nlat)
                push!(s_era5, collect(Float32, dT_era5[i, j, :]))
                push!(s_mod,  collect(Float32, dT_model[i, j, :]))
                push!(s_corr, collect(Float32, dT_corr[i, j, :]))
                push!(s_meta, @sprintf("%s i=%d j=%d", date_now, i, j))
            end
        end

        delete!(era5_cache, date_now)
    end

    n_valid > 0 || error("No valid day pairs — check ERA5_DIR_DAILY and dates.")
    z_km = Float32.(grid.z_ref[z_idx] ./ 1000.0)
    inv_n = 1.0 / max(n_prof, 1)
    mn(a) = Float32.(a .* inv_n)

    return SanityResults(
        z_km, nlon, nlat, nz, grid.lon, grid.lat,
        mn(acc_dT_era5), mn(acc_dT_model), mn(acc_dT_corr),
        mn(acc_dq_era5), mn(acc_dq_model), mn(acc_dq_corr),
        Float64[], Float64[], 0.0, 0.0,
        n_valid, n_skipped, false,
        acc_dT_corr_zonal, acc_dq_corr_zonal,
        acc_dT_corr_sq, acc_dq_corr_sq, n_prof,
        s_era5, s_mod, s_corr, s_meta,
        dT_snap, dq_snap, snap_label,
        ts_dT, ts_dq, ts_dt,
        hist_levels, hist_dT,
        "daily-era5",
    )
end

# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────
function next_sanity_index(base_dir)
    isdir(base_dir) || return 1
    existing = filter(d -> occursin(r"^\d{3}_", d), readdir(base_dir))
    isempty(existing) && return 1
    max_idx = maximum(parse(Int, m.match) for d in existing for m in eachmatch(r"^(\d{3})", d))
    return max_idx + 1
end

function make_sanity_dir()
    base = joinpath(@__DIR__, "sanity_tendency_plots")
    idx = next_sanity_index(base)
    src = DATA_SOURCE == :daily_era5 ? "daily_era5" :
          TEMPORAL_AVERAGING == :daily ? "daily" : "hourly"
    name = @sprintf("%03d_%s_d%d-%d", idx, src, TIME_DAY_START, Int(min(TIME_DAY_END, 999)))
    path = joinpath(base, name)
    mkpath(path)
    return path
end

function main()
    out_dir = length(ARGS) >= 1 ? ARGS[1] : make_sanity_dir()
    mkpath(out_dir)
    println("Output directory: ", out_dir)
    println("Data source:      ", DATA_SOURCE)
    println("Temporal mode:    ", TEMPORAL_AVERAGING)

    println("Loading model grid...")
    grid = load_model_grid()
    nlon, nlat = length(grid.lon), length(grid.lat)
    nz_full = length(grid.z_ref)
    model_dates = grid.dates
    nt = length(model_dates)

    z_mask = grid.z_ref .<= Z_MAX
    z_idx = findall(z_mask)
    nz = length(z_idx)

    ti_start, ti_end = time_indices(model_dates, nt)
    println("  Levels below Z_MAX: $nz / $nz_full")
    println("  Time pairs: ti = $ti_start : $ti_end ($(model_dates[ti_start]) → $(model_dates[ti_end+1]))")

    println("Loading ta, hus...")
    ta_mod = load_model_var("ta")
    hus_mod = load_model_var("hus")

    r = if DATA_SOURCE == :daily_era5
        process_daily_era5(grid, ta_mod, hus_mod, z_idx, nz, ti_start, ti_end, model_dates)
    elseif TEMPORAL_AVERAGING == :daily
        process_daily(grid, ta_mod, hus_mod, z_idx, nz, ti_start, ti_end, model_dates)
    else
        process_hourly(grid, ta_mod, hus_mod, z_idx, nz, ti_start, ti_end, model_dates)
    end

    # Print summary
    println()
    println("=== $(r.mode_label) statistics ===")
    @printf("  Valid %s: %d   Skipped ERA5 hours: %d\n",
            r.mode_label == "daily-mean" ? "days" : "pairs", r.n_valid, r.n_skipped)

    if r.has_telescoping
        if r.n_skipped > 0
            @warn "Skipped hours break the telescoping identity — expect nonzero mismatch."
        end
        println()
        println("=== Telescoping identity ===")
        @printf("  max |sum(dT_corr*dt) - direct_delta_T| : %.3e K\n", r.max_abs_diff_T)
        @printf("  max |sum(dq_corr*dt) - direct_delta_q| : %.3e kg/kg\n", r.max_abs_diff_q)
    end

    open(joinpath(out_dir, "sanity_summary.txt"), "w") do io
        println(io, "sanity_check_compute_tendencies.jl")
        println(io, "mode=", r.mode_label)
        println(io, "MODEL_DIR=", MODEL_DIR)
        println(io, "ERA5_DIR=", ERA5_DIR)
        println(io, "TIME_DAY_START=", TIME_DAY_START, " TIME_DAY_END=", TIME_DAY_END)
        println(io, "Z_MAX=", Z_MAX, " nz=", r.nz)
        println(io, "valid=", r.n_valid, " skipped_hours=", r.n_skipped)
        if r.has_telescoping
            @printf(io, "max_abs_diff_T_K %.6e\n", r.max_abs_diff_T)
            @printf(io, "max_abs_diff_q_kgkg %.6e\n", r.max_abs_diff_q)
        end
    end

    make_plots(r, out_dir)
    println("\nDone.")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
