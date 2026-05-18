#!/usr/bin/env julia
#=
Evaluate and visualise a trained flux correction model.

Loads a saved .bson model, rebuilds the dataset (or a subset), runs inference,
and produces diagnostic plots comparing predicted vs target tendencies, flux
profiles, and spatial error maps.

Usage:
    julia --project=. evaluate_flux_model.jl                          # latest run
    julia --project=. evaluate_flux_model.jl runs/001_cnn_h128_lr1e-3 # specific run dir
    julia --project=. evaluate_flux_model.jl path/to/model.bson       # direct model path
=#

include("train_flux_correction.jl")

using CairoMakie

# ─────────────────────────────────────────────────────────────────────────────
# Configuration
# ─────────────────────────────────────────────────────────────────────────────

function resolve_run_dir()
    if length(ARGS) >= 1
        arg = ARGS[1]
        if isdir(arg)
            return arg
        elseif endswith(arg, ".bson") && isfile(arg)
            return dirname(abspath(arg))
        end
    end
    # Default: read latest_run.txt written by training
    latest_path = joinpath(@__DIR__, "latest_run.txt")
    if isfile(latest_path)
        d = strip(read(latest_path, String))
        if isdir(d)
            return d
        end
    end
    # Fallback: highest-numbered run directory
    if isdir(RUNS_DIR)
        dirs = sort(filter(d -> occursin(r"^\d{3}_", d), readdir(RUNS_DIR)))
        if !isempty(dirs)
            return joinpath(RUNS_DIR, dirs[end])
        end
    end
    return @__DIR__
end

const RUN_DIR    = resolve_run_dir()
const MODEL_PATH = joinpath(RUN_DIR, "flux_correction_model.bson")
const PLOT_DIR   = joinpath(RUN_DIR, "plots")
const N_EXAMPLE_COLS = 8

# ─────────────────────────────────────────────────────────────────────────────
# Load model
# ─────────────────────────────────────────────────────────────────────────────

function load_trained_model(path)
    println("Loading model from $path ...")
    data = BSON.load(path)
    ps = data[:ps]
    st = data[:st]
    X_norm_stats = data[:X_norm_stats]
    Y_norm_stats = get(data, :Y_norm_stats, nothing)
    nz = data[:nz]
    input_dim = data[:input_dim]
    arch = get(data, :arch, ARCHITECTURE)
    temporal_avg = get(data, :temporal_averaging, TEMPORAL_AVERAGING)
    pm = get(data, :predict_mode, :flux)
    tv = get(data, :target_vars, :both)
    model = if arch == :unet
        build_unet_model(input_dim, nz; predict_mode = pm)
    elseif arch == :cnn
        build_cnn_model(input_dim, nz; predict_mode = pm)
    else
        build_mlp_model(input_dim, nz; predict_mode = pm)
    end
    println("  Architecture: $arch, predict_mode=$pm, target_vars=$tv, input_dim=$input_dim, nz=$nz")
    return model, ps, st, X_norm_stats, Y_norm_stats, nz, temporal_avg, pm, tv
end

# ─────────────────────────────────────────────────────────────────────────────
# Compute predicted tendencies from model output
# ─────────────────────────────────────────────────────────────────────────────

"""
Run the model on input X and produce cell-centre tendency predictions.
Returns a Dict mapping variable symbol (:T and/or :q) to (nz, n_samples) arrays
in physical units, plus a Dict of face fluxes (nothing entries in direct mode).
"""
function predict_tendencies(model, ps, st, X, dz, nz, tv)
    st_test = Lux.testmode(st)
    X_norm = normalise(X, X_norm_stats_global)
    raw_pred, _ = Lux.apply(model, X_norm, ps, st_test)

    nv = tv == :both ? 2 : 1

    if predict_mode_global == :direct
        pred_norm = raw_pred
        fluxes = Dict{Symbol, Any}()
    else
        nf = nz - 1
        z_row = raw_pred[1:1, :] .* 0f0
        pred_parts = Matrix{Float32}[]
        for v in 1:nv
            F_int = raw_pred[(v-1)*nf+1 : v*nf, :]
            F_full = vcat(z_row, F_int, z_row)
            push!(pred_parts, -(F_full[2:nz+1, :] .- F_full[1:nz, :]) ./ dz)
        end
        pred_norm = vcat(pred_parts...)
        fluxes = Dict{Symbol, Any}()
    end

    if Y_norm_stats_global !== nothing
        pred_phys = denormalise(pred_norm, Y_norm_stats_global)
    else
        pred_phys = pred_norm
    end

    results = Dict{Symbol, Matrix{Float32}}()
    var_order = tv == :both ? [:T, :q] : tv == :T ? [:T] : [:q]
    for (vi, sym) in enumerate(var_order)
        results[sym] = pred_phys[(vi-1)*nz+1 : vi*nz, :]
    end

    return results
end

# ─────────────────────────────────────────────────────────────────────────────
# Plotting functions
# ─────────────────────────────────────────────────────────────────────────────

function get_z_centres(nz)
    grid = load_model_grid()
    z_mask = grid.z_ref .<= Z_MAX
    z_idx = findall(z_mask)
    z_km = grid.z_ref[z_idx] ./ 1000.0
    return z_km
end

function get_z_faces(z_centres_km)
    nz = length(z_centres_km)
    z_faces = zeros(nz + 1)
    z_faces[1] = z_centres_km[1] - (z_centres_km[2] - z_centres_km[1]) / 2
    z_faces[end] = z_centres_km[end] + (z_centres_km[end] - z_centres_km[end-1]) / 2
    for k in 2:nz
        z_faces[k] = (z_centres_km[k-1] + z_centres_km[k]) / 2
    end
    return z_faces
end

const VAR_INFO = Dict(
    :T => (label = "Temperature", xlabel_tendency = "dT/dt (K/hr)", unit = "K/hr",
           color = :firebrick, scale = 3600.0f0),
    :q => (label = "Moisture", xlabel_tendency = "dq/dt (g/kg/hr)", unit = "g/kg/hr",
           color = :steelblue, scale = 3600.0f0 * 1000.0f0),
)

"""Plot 1: Mean vertical profiles of target vs predicted tendency."""
function plot_mean_profiles(z_km, targets, preds, vars; temporal_avg = :hourly)
    nv = length(vars)
    fig = Figure(size = (500 * nv, 500))
    avg_label = temporal_avg == :daily ? " (daily mean)" : ""

    for (i, sym) in enumerate(vars)
        info = VAR_INFO[sym]
        ax = Axis(fig[1, i]; xlabel = info.xlabel_tendency, ylabel = "Height (km)",
                  title = "$(info.label) tendency correction" * avg_label)
        tgt_mean = vec(mean(targets[sym], dims = 2)) .* info.scale
        prd_mean = vec(mean(preds[sym], dims = 2)) .* info.scale
        lines!(ax, tgt_mean, z_km; label = "Target (ERA5 - model)", linewidth = 2)
        lines!(ax, prd_mean, z_km; label = "Predicted (-dF/dz)", linewidth = 2, linestyle = :dash)
        axislegend(ax; position = :rt)
    end
    return fig
end

"""Plot 3: Example individual column profiles (target vs predicted)."""
function plot_example_columns(z_km, targets, preds, vars, n_examples; temporal_avg = :hourly)
    rng = Random.MersenneTwister(99)
    first_var = first(vars)
    n_samples = size(targets[first_var], 2)
    idx = randperm(rng, n_samples)[1:min(n_examples, n_samples)]
    avg_label = temporal_avg == :daily ? " (daily mean)" : ""

    ncols = min(4, length(idx))
    nrows = cld(length(idx), ncols)
    fig = Figure(size = (300 * ncols, 350 * nrows))

    for (panel, ci) in enumerate(idx)
        r = cld(panel, ncols)
        c = mod1(panel, ncols)
        info = VAR_INFO[first_var]
        ax = Axis(fig[r, c]; xlabel = info.unit, ylabel = "Height (km)",
                  title = "Column $ci")
        lines!(ax, targets[first_var][:, ci] .* info.scale, z_km; label = "Target", linewidth = 1.5)
        lines!(ax, preds[first_var][:, ci] .* info.scale, z_km; label = "Predicted",
               linewidth = 1.5, linestyle = :dash)
        panel == 1 && axislegend(ax; position = :rt, labelsize = 10)
    end

    Label(fig[0, :], "$(VAR_INFO[first_var].label) tendency: individual columns" * avg_label; fontsize = 16)
    return fig
end

"""Plot 4: Scatter plot of predicted vs target tendencies."""
function plot_scatter(targets, preds, vars; temporal_avg = :hourly)
    nv = length(vars)
    fig = Figure(size = (500 * nv, 500))
    avg_label = temporal_avg == :daily ? " (daily mean)" : ""

    for (i, sym) in enumerate(vars)
        info = VAR_INFO[sym]
        n = length(targets[sym])
        max_pts = 50_000
        stride = max(1, n ÷ max_pts)

        ax = Axis(fig[1, i]; xlabel = "Target " * info.xlabel_tendency,
                  ylabel = "Predicted " * info.xlabel_tendency,
                  title = "$(info.label) tendency" * avg_label, aspect = 1)
        tgt = vec(targets[sym])[1:stride:end] .* info.scale
        prd = vec(preds[sym])[1:stride:end] .* info.scale
        scatter!(ax, tgt, prd; markersize = 1, color = (info.color, 0.15))
        lims = extrema(vcat(tgt, prd))
        lines!(ax, [lims...], [lims...]; color = :black, linewidth = 1)
    end
    return fig
end

"""Plot 5: Vertical profile of RMSE and bias."""
function plot_error_profiles(z_km, targets, preds, vars; temporal_avg = :hourly)
    nv = length(vars)
    fig = Figure(size = (500 * nv, 500))
    avg_label = temporal_avg == :daily ? " (daily mean)" : ""

    for (i, sym) in enumerate(vars)
        info = VAR_INFO[sym]
        err = preds[sym] .- targets[sym]
        rmse = sqrt.(vec(mean(err .^ 2, dims = 2))) .* info.scale
        bias = vec(mean(err, dims = 2)) .* info.scale

        ax = Axis(fig[1, i]; xlabel = info.unit, ylabel = "Height (km)",
                  title = "$(info.label) error profile" * avg_label)
        lines!(ax, rmse, z_km; label = "RMSE", linewidth = 2, color = info.color)
        lines!(ax, bias, z_km; label = "Bias", linewidth = 2, color = info.color, linestyle = :dash)
        vlines!(ax, [0]; color = :gray, linestyle = :dot)
        axislegend(ax; position = :rt)
    end
    return fig
end

"""Plot 6: R² skill score as a function of height."""
function plot_r2_profile(z_km, targets, preds, vars; temporal_avg = :hourly)
    fig = Figure(size = (600, 500))
    avg_label = temporal_avg == :daily ? " (daily mean)" : ""
    ax = Axis(fig[1, 1]; xlabel = "R² score", ylabel = "Height (km)",
              title = "Skill score by level" * avg_label)

    for sym in vars
        info = VAR_INFO[sym]
        r2 = zeros(length(z_km))
        for k in eachindex(z_km)
            ss_res = sum((preds[sym][k, :] .- targets[sym][k, :]) .^ 2)
            ss_tot = sum((targets[sym][k, :] .- mean(targets[sym][k, :])) .^ 2)
            r2[k] = ss_tot > 0 ? 1 - ss_res / ss_tot : 0.0
        end
        lines!(ax, r2, z_km; label = info.label, linewidth = 2, color = info.color)
    end

    vlines!(ax, [0]; color = :gray, linestyle = :dot)
    vlines!(ax, [1]; color = :gray, linestyle = :dot)
    xlims!(ax, -0.5, 1.05)
    axislegend(ax; position = :lb)
    return fig
end

"""Plot 7: Magnitude of correction — RMS target, predicted, and residual."""
function plot_tendency_magnitudes(z_km, targets, preds, vars; temporal_avg = :hourly)
    nv = length(vars)
    fig = Figure(size = (500 * nv, 500))
    avg_label = temporal_avg == :daily ? " (daily mean)" : ""

    for (i, sym) in enumerate(vars)
        info = VAR_INFO[sym]
        rms_tgt = sqrt.(vec(mean(targets[sym] .^ 2, dims = 2))) .* info.scale
        rms_prd = sqrt.(vec(mean(preds[sym] .^ 2, dims = 2))) .* info.scale
        rms_res = sqrt.(vec(mean((preds[sym] .- targets[sym]) .^ 2, dims = 2))) .* info.scale

        ax = Axis(fig[1, i]; xlabel = "RMS tendency (" * info.unit * ")", ylabel = "Height (km)",
                  title = "$(info.label): correction magnitude" * avg_label)
        lines!(ax, rms_tgt, z_km; label = "Target correction", linewidth = 2)
        lines!(ax, rms_prd, z_km; label = "Predicted correction", linewidth = 2, linestyle = :dash)
        lines!(ax, rms_res, z_km; label = "Residual error", linewidth = 2, linestyle = :dot)
        axislegend(ax; position = :rt)
    end
    return fig
end

# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

function evaluate()
    println("=" ^ 60)
    println("Evaluating run: $RUN_DIR")
    println("=" ^ 60)
    mkpath(PLOT_DIR)

    model, ps, st, X_norm_stats, Y_norm_stats, nz, temporal_avg, pm, tv = load_trained_model(MODEL_PATH)
    global X_norm_stats_global = X_norm_stats
    global Y_norm_stats_global = Y_norm_stats
    global predict_mode_global = pm

    vars = tv == :both ? [:T, :q] : tv == :T ? [:T] : [:q]

    println("\nBuilding evaluation dataset (temporal_averaging=$temporal_avg)...")
    cached = load_dataset_cache()
    X_raw, Y_raw, dz_all, nz_data = if cached !== nothing
        cached
    else
        result = build_dataset()
        save_dataset_cache(result...)
        result
    end
    @assert nz == nz_data "Model nz ($nz) != data nz ($nz_data)"
    println("Evaluation set: $(size(X_raw, 2)) columns")

    println("\nRunning inference...")
    preds = predict_tendencies(model, ps, st, X_raw, dz_all, nz, tv)

    nv = length(vars)
    targets = Dict{Symbol, Matrix{Float32}}()
    for (vi, sym) in enumerate(vars)
        targets[sym] = Y_raw[(vi-1)*nz+1 : vi*nz, :]
    end

    # Overall metrics
    open(joinpath(RUN_DIR, "metrics.txt"), "w") do io
        for sym in vars
            info = VAR_INFO[sym]
            rmse = sqrt(mean((preds[sym] .- targets[sym]) .^ 2)) * info.scale
            corr_val = cor(vec(preds[sym]), vec(targets[sym]))
            @printf("  RMSE %s: %.4f %s   Corr: %.4f\n", sym, rmse, info.unit, corr_val)
            @printf(io, "rmse_%s = %.6f\n", sym, rmse)
            @printf(io, "corr_%s = %.6f\n", sym, corr_val)
        end
        @printf(io, "temporal_avg = %s\n", temporal_avg)
        @printf(io, "target_vars  = %s\n", tv)
    end

    z_km = get_z_centres(nz)

    println("\nGenerating plots...")

    fig1 = plot_mean_profiles(z_km, targets, preds, vars; temporal_avg)
    save(joinpath(PLOT_DIR, "01_mean_tendency_profiles.png"), fig1; px_per_unit = 2)
    println("  01_mean_tendency_profiles.png")

    fig3 = plot_example_columns(z_km, targets, preds, vars, N_EXAMPLE_COLS; temporal_avg)
    save(joinpath(PLOT_DIR, "03_example_columns.png"), fig3; px_per_unit = 2)
    println("  03_example_columns.png")

    fig4 = plot_scatter(targets, preds, vars; temporal_avg)
    save(joinpath(PLOT_DIR, "04_scatter_pred_vs_target.png"), fig4; px_per_unit = 2)
    println("  04_scatter_pred_vs_target.png")

    fig5 = plot_error_profiles(z_km, targets, preds, vars; temporal_avg)
    save(joinpath(PLOT_DIR, "05_error_profiles.png"), fig5; px_per_unit = 2)
    println("  05_error_profiles.png")

    fig6 = plot_r2_profile(z_km, targets, preds, vars; temporal_avg)
    save(joinpath(PLOT_DIR, "06_r2_skill_score.png"), fig6; px_per_unit = 2)
    println("  06_r2_skill_score.png")

    fig7 = plot_tendency_magnitudes(z_km, targets, preds, vars; temporal_avg)
    save(joinpath(PLOT_DIR, "07_tendency_magnitudes.png"), fig7; px_per_unit = 2)
    println("  07_tendency_magnitudes.png")

    println("\nAll plots saved to $(PLOT_DIR)/")
    println("Done.")
end

if abspath(PROGRAM_FILE) == @__FILE__
    evaluate()
end
