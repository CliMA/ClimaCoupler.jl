#!/usr/bin/env julia
#=
Sanity checks for the offline *computed tendencies* used in flux-correction
training (ERA5 − model hourly finite differences).

Reuses the same loaders and time/level window as `train_flux_correction.jl`
(MODEL_DIR, ERA5_DIR, Z_MAX, TIME_DAY_START/END, DT).

Checks
  1) Telescoping identity: sum over valid hourly pairs of
     dT_corr*DT  should match  (T_ERA5[last]−T_ERA5[first]) − (T_model[last]−T_model[first])
     (and similarly for humidity), where “first/last” are the endpoints of the
     chain of valid pairs (skips break the calendar identity — we report skips).

Plots (written to output dir, default `ml_flux/sanity_tendency_plots/`):
  - Mean vertical profiles of dT/dt from ERA5, model, and correction (q too)
  - Several sample column profiles (random columns / timesteps)
  - Profile of RMS(|sum(dT_corr*DT) − direct_delta|) vs height (column-wise check)
  - Spatial maps of dT_corr at selected levels for one snapshot hour

Usage:
    julia --project=. sanity_check_compute_tendencies.jl
    julia --project=. sanity_check_compute_tendencies.jl /path/to/output_dir
=#

include("train_flux_correction.jl")

using CairoMakie
using Printf
using Random

function time_indices(model_dates, nt)
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

function main()
    out_dir = length(ARGS) >= 1 ? ARGS[1] : joinpath(@__DIR__, "sanity_tendency_plots")
    mkpath(out_dir)
    println("Output directory: ", out_dir)

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
    ta_model = load_model_var("ta")
    hus_model = load_model_var("hus")

    # Accumulators for mean profiles (all columns, valid pairs only)
    acc_dT_era5 = zeros(Float64, nz)
    acc_dT_model = zeros(Float64, nz)
    acc_dT_corr = zeros(Float64, nz)
    acc_dq_era5 = zeros(Float64, nz)
    acc_dq_model = zeros(Float64, nz)
    acc_dq_corr = zeros(Float64, nz)
    n_prof = Ref(0)

    # Integral identity: sum of d*corr * DT over valid pairs
    sum_dT_corr_dt = zeros(Float64, nlon, nlat, nz)
    sum_dq_corr_dt = zeros(Float64, nlon, nlat, nz)
    n_skipped = Ref(0)
    n_valid_pairs = Ref(0)

    # Chain endpoints for direct delta (updated each valid pair)
    era5_t0 = zeros(Float32, nlon, nlat, nz)
    era5_t1 = zeros(Float32, nlon, nlat, nz)
    q_era5_t0 = zeros(Float32, nlon, nlat, nz)
    q_era5_t1 = zeros(Float32, nlon, nlat, nz)
    model_t0 = zeros(Float32, nlon, nlat, nz)
    model_t1 = zeros(Float32, nlon, nlat, nz)
    hus_t0 = zeros(Float32, nlon, nlat, nz)
    hus_t1 = zeros(Float32, nlon, nlat, nz)
    chain_started = Ref(false)

    # Sample profiles: store (dT_era5, dT_model, dT_corr) columns for a few draws
    rng = MersenneTwister(42)
    n_samples_wanted = 12
    sample_dT_era5 = Vector{Vector{Float32}}()
    sample_dT_mod = Vector{Vector{Float32}}()
    sample_dT_corr = Vector{Vector{Float32}}()
    sample_meta = Vector{String}()

    # Zonal-mean accumulator: (nlat, nz) for latitude-height cross-section
    acc_dT_corr_zonal = zeros(Float64, nlat, nz)
    acc_dq_corr_zonal = zeros(Float64, nlat, nz)

    # Per-level RMS accumulator
    acc_dT_corr_sq = zeros(Float64, nz)
    acc_dq_corr_sq = zeros(Float64, nz)

    # Temporal evolution: mean correction per valid timestep
    ts_dT_corr_mean = Float64[]
    ts_dq_corr_mean = Float64[]
    ts_dates = DateTime[]

    # Histogram data: store all correction values at a few representative levels
    hist_levels = level_indices_for_maps(nz)
    hist_dT_data = Dict{Int, Vector{Float32}}(k => Float32[] for k in hist_levels)

    # Snapshot for spatial maps: use valid timestep closest to ti_mid (ti_mid itself may be skipped)
    ti_mid = (ti_start + ti_end) ÷ 2
    dT_corr_snapshot = zeros(Float32, nlon, nlat, nz)
    dq_corr_snapshot = zeros(Float32, nlon, nlat, nz)
    snapshot_time = Ref{Union{Nothing, DateTime}}(nothing)
    snapshot_ti = Ref{Int}(ti_mid)
    best_snap_dist = Ref(typemax(Int))

    # Cache regridded/interpolated ERA5 fields: the "next" of step ti becomes
    # the "now" of step ti+1, cutting ERA5 I/O roughly in half.
    prev_dt_next = nothing                        # DateTime of cached "next"
    cache_t_next_z = Array{Float32}(undef, 0, 0, 0)
    cache_q_next_z = Array{Float32}(undef, 0, 0, 0)

    function regrid_era5_to_model_z(era5, grid, z_idx)
        t_rg = regrid_horizontal(era5.t, era5.lon, era5.lat, grid.lon, grid.lat)
        q_rg = regrid_horizontal(era5.q, era5.lon, era5.lat, grid.lon, grid.lat)
        z_rg = regrid_horizontal(era5.z, era5.lon, era5.lat, grid.lon, grid.lat)
        t_on_z = interp_vertical_3d(t_rg, z_rg, grid.z_phys)[:, :, z_idx]
        q_on_z = interp_vertical_3d(q_rg, z_rg, grid.z_phys)[:, :, z_idx]
        return t_on_z, q_on_z
    end

    for ti in ti_start:ti_end
        dt_now = model_dates[ti]
        dt_next = model_dates[ti + 1]

        era5_now = load_era5_pressure(dt_now)
        era5_next = load_era5_pressure(dt_next)
        if era5_now === nothing || era5_next === nothing
            n_skipped[] += 1
            prev_dt_next = nothing   # invalidate cache on skip
            continue
        end
        n_valid_pairs[] += 1

        if ti == ti_start || ti % 24 == 0
            @printf("  timestep %d / %d  %s\n", ti, ti_end, dt_now)
            flush(stdout)
        end

        # Reuse cached "next" from previous iteration if available
        if prev_dt_next == dt_now
            t_now_z = cache_t_next_z
            q_now_z = cache_q_next_z
        else
            t_now_z, q_now_z = regrid_era5_to_model_z(era5_now, grid, z_idx)
        end

        t_next_z, q_next_z = regrid_era5_to_model_z(era5_next, grid, z_idx)

        # Cache this "next" for the following iteration
        cache_t_next_z = t_next_z
        cache_q_next_z = q_next_z
        prev_dt_next = dt_next

        ta_now_z = ta_model[ti, :, :, z_idx]
        ta_next_z = ta_model[ti + 1, :, :, z_idx]
        hus_now_z = hus_model[ti, :, :, z_idx]
        hus_next_z = hus_model[ti + 1, :, :, z_idx]

        dT_era5 = (t_next_z .- t_now_z) ./ DT
        dq_era5 = (q_next_z .- q_now_z) ./ DT
        dT_model = (ta_next_z .- ta_now_z) ./ DT
        dq_model = (hus_next_z .- hus_now_z) ./ DT
        dT_corr = dT_era5 .- dT_model
        dq_corr = dq_era5 .- dq_model

        # Mean profiles
        n_prof[] += nlon * nlat
        for k in 1:nz
            acc_dT_era5[k] += sum(Float64, @view dT_era5[:, :, k])
            acc_dT_model[k] += sum(Float64, @view dT_model[:, :, k])
            acc_dT_corr[k] += sum(Float64, @view dT_corr[:, :, k])
            acc_dq_era5[k] += sum(Float64, @view dq_era5[:, :, k])
            acc_dq_model[k] += sum(Float64, @view dq_model[:, :, k])
            acc_dq_corr[k] += sum(Float64, @view dq_corr[:, :, k])
        end

        # Zonal-mean cross-section (average over longitudes)
        for k in 1:nz, j in 1:nlat
            acc_dT_corr_zonal[j, k] += sum(Float64, @view dT_corr[:, j, k])
            acc_dq_corr_zonal[j, k] += sum(Float64, @view dq_corr[:, j, k])
        end

        # Per-level sum-of-squares for RMS
        for k in 1:nz
            acc_dT_corr_sq[k] += mapreduce(x -> Float64(x)^2, +, @view dT_corr[:, :, k])
            acc_dq_corr_sq[k] += mapreduce(x -> Float64(x)^2, +, @view dq_corr[:, :, k])
        end

        # Temporal evolution
        push!(ts_dT_corr_mean, Float64(mean(dT_corr)))
        push!(ts_dq_corr_mean, Float64(mean(dq_corr)))
        push!(ts_dates, dt_now)

        # Histogram samples at representative levels (subsample to keep memory bounded)
        for k in hist_levels
            slab = @view dT_corr[:, :, k]
            stride = max(1, length(slab) ÷ 500)
            append!(hist_dT_data[k], slab[1:stride:end])
        end

        sum_dT_corr_dt .+= Float64.(dT_corr) .* Float64(DT)
        sum_dq_corr_dt .+= Float64.(dq_corr) .* Float64(DT)

        if !chain_started[]
            era5_t0 .= t_now_z
            q_era5_t0 .= q_now_z
            model_t0 .= ta_now_z
            hus_t0 .= hus_now_z
            chain_started[] = true
        end
        era5_t1 .= t_next_z
        q_era5_t1 .= q_next_z
        model_t1 .= ta_next_z
        hus_t1 .= hus_next_z

        dist_snap = abs(ti - ti_mid)
        if dist_snap < best_snap_dist[]
            best_snap_dist[] = dist_snap
            dT_corr_snapshot .= dT_corr
            dq_corr_snapshot .= dq_corr
            snapshot_time[] = dt_now
            snapshot_ti[] = ti
        end

        # Random sample columns for profile plots (spread across timesteps)
        if length(sample_dT_corr) < n_samples_wanted
            n_draw = min(2, n_samples_wanted - length(sample_dT_corr))
            for _ in 1:n_draw
                i = rand(rng, 1:nlon)
                j = rand(rng, 1:nlat)
                push!(sample_dT_era5, collect(Float32, dT_era5[i, j, :]))
                push!(sample_dT_mod, collect(Float32, dT_model[i, j, :]))
                push!(sample_dT_corr, collect(Float32, dT_corr[i, j, :]))
                push!(sample_meta, @sprintf("ti=%d i=%d j=%d %s", ti, i, j, dt_now))
            end
        end
    end

    if !chain_started[]
        error("No valid ERA5 pairs in the requested time window — check ERA5_DIR and dates.")
    end

    z_km = Float32.(grid.z_ref[z_idx] ./ 1000.0)

    mean_dT_era5 = Float32.(acc_dT_era5 ./ max(n_prof[], 1))
    mean_dT_model = Float32.(acc_dT_model ./ max(n_prof[], 1))
    mean_dT_corr = Float32.(acc_dT_corr ./ max(n_prof[], 1))
    mean_dq_era5 = Float32.(acc_dq_era5 ./ max(n_prof[], 1))
    mean_dq_model = Float32.(acc_dq_model ./ max(n_prof[], 1))
    mean_dq_corr = Float32.(acc_dq_corr ./ max(n_prof[], 1))

    # Direct residual over the integrated chain
    direct_dT = era5_t1 .- era5_t0 .- (model_t1 .- model_t0)
    direct_dq = q_era5_t1 .- q_era5_t0 .- (hus_t1 .- hus_t0)

    diff_T = sum_dT_corr_dt .- Float64.(direct_dT)
    diff_q = sum_dq_corr_dt .- Float64.(direct_dq)

    rms_diff_T_z = [sqrt(mean(diff_T[:, :, k] .^ 2)) for k in 1:nz]
    rms_diff_q_z = [sqrt(mean(diff_q[:, :, k] .^ 2)) for k in 1:nz]
    max_abs_diff_T = maximum(abs, diff_T)
    max_abs_diff_q = maximum(abs, diff_q)

    println()
    println("=== Hourly pair statistics ===")
    @printf("  Valid hourly pairs: %d   Skipped (missing ERA5): %d\n", n_valid_pairs[], n_skipped[])
    if n_skipped[] > 0
        @warn "Skipped model hours break the telescoping check: sum(d*corr*dt) uses only processed pairs, while direct delta uses ERA5/model at the first and last valid pair endpoints. Expect nonzero mismatch if gaps exist."
    end

    println()
    println("=== Telescoping identity (full horizontal column, all levels) ===")
    @printf("  max |sum(dT_corr*dt) - direct_delta_T| : %.3e K\n", max_abs_diff_T)
    @printf("  max |sum(dq_corr*dt) - direct_delta_q| : %.3e kg/kg\n", max_abs_diff_q)
    @printf("  RMS of column mismatch, T (all lev): %.3e K\n", sqrt(mean(diff_T .^ 2)))
    @printf("  RMS of column mismatch, q (all lev): %.3e kg/kg\n", sqrt(mean(diff_q .^ 2)))

    open(joinpath(out_dir, "sanity_summary.txt"), "w") do io
        println(io, "sanity_check_compute_tendencies.jl")
        println(io, "MODEL_DIR=", MODEL_DIR)
        println(io, "ERA5_DIR=", ERA5_DIR)
        println(io, "TIME_DAY_START=", TIME_DAY_START, " TIME_DAY_END=", TIME_DAY_END)
        println(io, "Z_MAX=", Z_MAX, " nz=", nz)
        println(io, "ti_start=", ti_start, " ti_end=", ti_end)
        println(io, "valid_hourly_pairs=", n_valid_pairs[])
        println(io, "skipped_era5_pairs=", n_skipped[])
        println(io, "snapshot_ti=", snapshot_ti[])
        println(io, "note_telescoping=Exact only if every hour in [ti_start,ti_end] has ERA5 pairs (no skips).")
        @printf(io, "max_abs_diff_T_K %.6e\n", max_abs_diff_T)
        @printf(io, "max_abs_diff_q_kgkg %.6e\n", max_abs_diff_q)
        @printf(io, "rms_diff_T_K %.6e\n", sqrt(mean(diff_T .^ 2)))
        @printf(io, "rms_diff_q_kgkg %.6e\n", sqrt(mean(diff_q .^ 2)))
    end

    # ── Plot: mean vertical profiles (K/hr, g/kg/hr) ─────────────────────
    fig_mean = Figure(size = (1000, 500))
    ax1 = Axis(
        fig_mean[1, 1];
        xlabel = "dT/dt (K/hr)",
        ylabel = "Height (km)",
        title = "Mean vertical profiles (computed hourly)",
    )
    lines!(ax1, mean_dT_era5 .* 3600, z_km; label = "ERA5 dT/dt", linewidth = 2)
    lines!(ax1, mean_dT_model .* 3600, z_km; label = "Model dT/dt", linewidth = 2)
    lines!(ax1, mean_dT_corr .* 3600, z_km; label = "Correction (ERA5−model)", linewidth = 2, linestyle = :dash)
    axislegend(ax1; position = :rt)

    ax2 = Axis(
        fig_mean[1, 2];
        xlabel = "dq/dt (g/kg/hr)",
        ylabel = "Height (km)",
        title = "Mean moisture tendencies",
    )
    lines!(ax2, mean_dq_era5 .* 3600 .* 1000, z_km; label = "ERA5", linewidth = 2)
    lines!(ax2, mean_dq_model .* 3600 .* 1000, z_km; label = "Model", linewidth = 2)
    lines!(ax2, mean_dq_corr .* 3600 .* 1000, z_km; label = "Correction", linewidth = 2, linestyle = :dash)
    axislegend(ax2; position = :rt)

    save(joinpath(out_dir, "01_mean_tendency_profiles.png"), fig_mean; px_per_unit = 2)
    println("Wrote ", joinpath(out_dir, "01_mean_tendency_profiles.png"))

    # ── Plot: telescoping error vs height ──────────────────────────────────
    fig_err = Figure(size = (700, 500))
    axe = Axis(
        fig_err[1, 1];
        xlabel = "RMS column error (K)",
        ylabel = "Height (km)",
        title = "RMS(|Σ dT_corr Δt − ΔT_residual|) vs z",
    )
    lines!(axe, rms_diff_T_z, z_km; linewidth = 2, color = :darkred)
    save(joinpath(out_dir, "02_integral_identity_rms_T.png"), fig_err; px_per_unit = 2)
    println("Wrote ", joinpath(out_dir, "02_integral_identity_rms_T.png"))

    fig_errq = Figure(size = (700, 500))
    axq = Axis(
        fig_errq[1, 1];
        xlabel = "RMS column error (kg/kg)",
        ylabel = "Height (km)",
        title = "RMS(|Σ dq_corr Δt − Δq_residual|) vs z",
    )
    lines!(axq, rms_diff_q_z, z_km; linewidth = 2, color = :steelblue)
    save(joinpath(out_dir, "03_integral_identity_rms_q.png"), fig_errq; px_per_unit = 2)
    println("Wrote ", joinpath(out_dir, "03_integral_identity_rms_q.png"))

    # ── Plot: sample columns ────────────────────────────────────────────────
    n_s = length(sample_dT_corr)
    if n_s > 0
        ncols = min(4, max(1, n_s))
        nrows = cld(n_s, ncols)
        fig_samp = Figure(size = (280 * ncols, 320 * nrows))
        for p in 1:n_s
            r = cld(p, ncols)
            c = mod1(p, ncols)
            ax = Axis(
                fig_samp[r, c];
                xlabel = "dT/dt (K/hr)",
                ylabel = c == 1 ? "z (km)" : "",
                title = sample_meta[p],
            )
            lines!(ax, sample_dT_era5[p] .* 3600, z_km; label = "ERA5", linewidth = 1.5)
            lines!(ax, sample_dT_mod[p] .* 3600, z_km; label = "Model", linewidth = 1.5)
            lines!(ax, sample_dT_corr[p] .* 3600, z_km; label = "Corr", linewidth = 1.5, linestyle = :dash)
            p == 1 && axislegend(ax; position = :rt, labelsize = 8)
        end
        Label(fig_samp[0, :], "Sample hourly dT/dt profiles"; fontsize = 14)
        save(joinpath(out_dir, "04_sample_column_profiles.png"), fig_samp; px_per_unit = 2)
        println("Wrote ", joinpath(out_dir, "04_sample_column_profiles.png"))
    else
        @warn "No sample columns collected; skipping 04_sample_column_profiles.png"
    end

    # ── Spatial maps with PER-LEVEL colorscale ──────────────────────────────
    ks = level_indices_for_maps(nz)
    nmaps = length(ks)
    lon = grid.lon
    lat = grid.lat
    snap_str = something(snapshot_time[], model_dates[ti_mid])
    snap_ti = snapshot_ti[]

    # T correction
    fig_map = Figure(size = (320 * nmaps, 320))
    Label(fig_map[0, :], "dT_corr (K/hr) snapshot: $snap_str (ti=$snap_ti)"; fontsize = 12)
    slab_Khr = dT_corr_snapshot .* 3600
    for (ip, k) in enumerate(ks)
        ax = Axis(
            fig_map[1, ip];
            aspect = DataAspect(),
            xlabel = "lon",
            ylabel = ip == 1 ? "lat" : "",
            title = @sprintf("z ≈ %.1f km (k=%d)", z_km[k], k),
        )
        level_data = slab_Khr[:, :, k]'
        clev = Float64(maximum(abs, level_data; init = 1f-6))
        hm = heatmap!(ax, lon, lat, level_data; colormap = :balance,
                       colorrange = (-clev, clev))
        Colorbar(fig_map[2, ip], hm; label = "K/hr", vertical = false,
                 flipaxis = false)
    end
    save(joinpath(out_dir, "05_spatial_dT_corr_levels.png"), fig_map; px_per_unit = 2)
    println("Wrote ", joinpath(out_dir, "05_spatial_dT_corr_levels.png"))

    # q correction
    fig_mapq = Figure(size = (320 * nmaps, 320))
    Label(fig_mapq[0, :], "dq_corr (g/kg/hr) snapshot: $snap_str (ti=$snap_ti)"; fontsize = 12)
    slab_qhr = dq_corr_snapshot .* 3600 .* 1000
    for (ip, k) in enumerate(ks)
        ax = Axis(
            fig_mapq[1, ip];
            aspect = DataAspect(),
            xlabel = "lon",
            ylabel = ip == 1 ? "lat" : "",
            title = @sprintf("z ≈ %.1f km (k=%d)", z_km[k], k),
        )
        level_data = slab_qhr[:, :, k]'
        clev = Float64(maximum(abs, level_data; init = 1f-6))
        hm = heatmap!(ax, lon, lat, level_data; colormap = :balance,
                       colorrange = (-clev, clev))
        Colorbar(fig_mapq[2, ip], hm; label = "g/kg/hr", vertical = false,
                 flipaxis = false)
    end
    save(joinpath(out_dir, "06_spatial_dq_corr_levels.png"), fig_mapq; px_per_unit = 2)
    println("Wrote ", joinpath(out_dir, "06_spatial_dq_corr_levels.png"))

    # ── Zonal-mean latitude-height cross-section ─────────────────────────
    n_zonal_denom = max(n_prof[] ÷ nlat, 1)  # nlon * n_valid_pairs
    zonal_dT = Float32.(acc_dT_corr_zonal ./ (n_zonal_denom)) .* 3600
    zonal_dq = Float32.(acc_dq_corr_zonal ./ (n_zonal_denom)) .* 3600 .* 1000

    fig_zonal = Figure(size = (1000, 450))
    ax_z1 = Axis(fig_zonal[1, 1]; xlabel = "Latitude", ylabel = "Height (km)",
                  title = "Zonal-mean dT_corr (K/hr)")
    clim_z1 = Float64(maximum(abs, zonal_dT; init = 1f-6))
    hm_z1 = heatmap!(ax_z1, Float32.(grid.lat), z_km, zonal_dT;
                      colormap = :balance, colorrange = (-clim_z1, clim_z1))
    Colorbar(fig_zonal[1, 2], hm_z1; label = "K/hr")

    ax_z2 = Axis(fig_zonal[1, 3]; xlabel = "Latitude", ylabel = "Height (km)",
                  title = "Zonal-mean dq_corr (g/kg/hr)")
    clim_z2 = Float64(maximum(abs, zonal_dq; init = 1f-6))
    hm_z2 = heatmap!(ax_z2, Float32.(grid.lat), z_km, zonal_dq;
                      colormap = :balance, colorrange = (-clim_z2, clim_z2))
    Colorbar(fig_zonal[1, 4], hm_z2; label = "g/kg/hr")
    save(joinpath(out_dir, "07_zonal_mean_correction.png"), fig_zonal; px_per_unit = 2)
    println("Wrote ", joinpath(out_dir, "07_zonal_mean_correction.png"))

    # ── RMS and std vertical profiles of correction ──────────────────────
    rms_dT_corr = Float32.(sqrt.(acc_dT_corr_sq ./ max(n_prof[], 1))) .* 3600
    rms_dq_corr = Float32.(sqrt.(acc_dq_corr_sq ./ max(n_prof[], 1))) .* 3600 .* 1000
    abs_mean_dT = abs.(mean_dT_corr) .* 3600
    abs_mean_dq = abs.(mean_dq_corr) .* 3600 .* 1000

    fig_rms = Figure(size = (1000, 500))
    ax_r1 = Axis(fig_rms[1, 1]; xlabel = "K/hr", ylabel = "Height (km)",
                  title = "T correction: RMS vs |mean|")
    lines!(ax_r1, rms_dT_corr, z_km; label = "RMS", linewidth = 2, color = :firebrick)
    lines!(ax_r1, abs_mean_dT, z_km; label = "|Mean|", linewidth = 2, color = :firebrick,
           linestyle = :dash)
    axislegend(ax_r1; position = :rt)

    ax_r2 = Axis(fig_rms[1, 2]; xlabel = "g/kg/hr", ylabel = "Height (km)",
                  title = "q correction: RMS vs |mean|")
    lines!(ax_r2, rms_dq_corr, z_km; label = "RMS", linewidth = 2, color = :steelblue)
    lines!(ax_r2, abs_mean_dq, z_km; label = "|Mean|", linewidth = 2, color = :steelblue,
           linestyle = :dash)
    axislegend(ax_r2; position = :rt)
    save(joinpath(out_dir, "08_rms_vs_mean_profiles.png"), fig_rms; px_per_unit = 2)
    println("Wrote ", joinpath(out_dir, "08_rms_vs_mean_profiles.png"))

    # ── Temporal evolution of mean correction ─────────────────────────────
    n_ts = length(ts_dates)
    if n_ts > 1
        hours = Float64[(ts_dates[i] - ts_dates[1]).value / 3600_000 for i in 1:n_ts]
        fig_ts = Figure(size = (1000, 450))
        ax_t1 = Axis(fig_ts[1, 1]; xlabel = "Hours since start", ylabel = "K/hr",
                      title = "Global-mean dT_corr vs time")
        lines!(ax_t1, hours, ts_dT_corr_mean .* 3600; linewidth = 1.5, color = :firebrick)
        hlines!(ax_t1, [0]; color = :gray, linestyle = :dot)

        ax_t2 = Axis(fig_ts[1, 2]; xlabel = "Hours since start", ylabel = "g/kg/hr",
                      title = "Global-mean dq_corr vs time")
        lines!(ax_t2, hours, ts_dq_corr_mean .* 3600 .* 1000; linewidth = 1.5, color = :steelblue)
        hlines!(ax_t2, [0]; color = :gray, linestyle = :dot)
        save(joinpath(out_dir, "09_temporal_evolution.png"), fig_ts; px_per_unit = 2)
        println("Wrote ", joinpath(out_dir, "09_temporal_evolution.png"))
    end

    # ── Histograms of correction at representative levels ─────────────────
    n_hist = length(hist_levels)
    fig_hist = Figure(size = (300 * min(n_hist, 5), 350))
    Label(fig_hist[0, :], "Distribution of dT_corr (K/hr) by level"; fontsize = 13)
    for (ip, k) in enumerate(hist_levels)
        vals = hist_dT_data[k] .* 3600
        ax = Axis(fig_hist[1, ip]; xlabel = "K/hr",
                   ylabel = ip == 1 ? "Count" : "",
                   title = @sprintf("z ≈ %.1f km", z_km[k]))
        hist!(ax, vals; bins = 80, color = (:firebrick, 0.6))
        vlines!(ax, [0]; color = :black, linestyle = :dot)
    end
    save(joinpath(out_dir, "10_histogram_dT_corr.png"), fig_hist; px_per_unit = 2)
    println("Wrote ", joinpath(out_dir, "10_histogram_dT_corr.png"))

    println("\nDone.")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
