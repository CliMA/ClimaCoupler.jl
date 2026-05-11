#!/usr/bin/env julia
#=
Evaluate and visualise a trained flux correction model.

Loads a saved .bson model, rebuilds the dataset (or a subset), runs inference,
and produces diagnostic plots comparing predicted vs target tendencies, flux
profiles, and spatial error maps.

Usage:
    julia --project=. evaluate_flux_model.jl                          # defaults
    julia --project=. evaluate_flux_model.jl path/to/model.bson       # custom model
=#

# Include shared code from the training script (config, data loaders, regridding,
# normalisation, model builders). The train() function won't run because we
# don't call it.
include("train_flux_correction.jl")

using CairoMakie

# ─────────────────────────────────────────────────────────────────────────────
# Configuration
# ─────────────────────────────────────────────────────────────────────────────
const MODEL_PATH = length(ARGS) >= 1 ? ARGS[1] : "flux_correction_model.bson"
const PLOT_DIR   = "plots"
const N_EXAMPLE_COLS = 8   # number of example columns for profile plots

# ─────────────────────────────────────────────────────────────────────────────
# Load model
# ─────────────────────────────────────────────────────────────────────────────

function load_trained_model(path)
    println("Loading model from $path ...")
    data = BSON.load(path)
    ps = data[:ps]
    st = data[:st]
    X_norm_stats = data[:X_norm_stats]
    nz = data[:nz]
    input_dim = data[:input_dim]
    arch = get(data, :arch, ARCHITECTURE)
    model = if arch == :unet
        build_unet_model(input_dim, nz)
    elseif arch == :cnn
        build_cnn_model(input_dim, nz)
    else
        build_mlp_model(input_dim, nz)
    end
    println("  Architecture: $arch, input_dim=$input_dim, nz=$nz")
    return model, ps, st, X_norm_stats, nz
end

# ─────────────────────────────────────────────────────────────────────────────
# Compute predicted tendencies from model output
# ─────────────────────────────────────────────────────────────────────────────

"""
Run the model on input X and convert face fluxes → cell-centre tendencies.
Returns (dT_pred, dq_pred) each of shape (nz, n_samples).
"""
function predict_tendencies(model, ps, st, X, dz, nz)
    X_norm = normalise(X, X_norm_stats_global)
    flux_pred, _ = Lux.apply(model, X_norm, ps, st)

    nf = nz - 1
    F_T_int = flux_pred[1:nf, :]
    F_q_int = flux_pred[nf+1:2*nf, :]

    z_row = F_T_int[1:1, :] .* 0f0
    F_T = vcat(z_row, F_T_int, z_row)
    F_q = vcat(z_row, F_q_int, z_row)

    dT_pred = -(F_T[2:nz+1, :] .- F_T[1:nz, :]) ./ dz
    dq_pred = -(F_q[2:nz+1, :] .- F_q[1:nz, :]) ./ dz

    return dT_pred, dq_pred, F_T, F_q
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

"""Plot 1: Mean vertical profiles of target vs predicted tendency."""
function plot_mean_profiles(z_km, dT_target, dT_pred, dq_target, dq_pred)
    fig = Figure(size = (1000, 500))

    ax1 = Axis(fig[1, 1]; xlabel = "dT/dt tendency (K/hr)", ylabel = "Height (km)",
               title = "Temperature tendency correction")
    dT_tgt_mean = vec(mean(dT_target, dims = 2)) .* 3600
    dT_prd_mean = vec(mean(dT_pred, dims = 2)) .* 3600
    lines!(ax1, dT_tgt_mean, z_km; label = "Target (ERA5 - model)", linewidth = 2)
    lines!(ax1, dT_prd_mean, z_km; label = "Predicted (-dF/dz)", linewidth = 2, linestyle = :dash)
    axislegend(ax1; position = :rt)

    ax2 = Axis(fig[1, 2]; xlabel = "dq/dt tendency (g/kg/hr)", ylabel = "Height (km)",
               title = "Moisture tendency correction")
    dq_tgt_mean = vec(mean(dq_target, dims = 2)) .* 3600 .* 1000
    dq_prd_mean = vec(mean(dq_pred, dims = 2)) .* 3600 .* 1000
    lines!(ax2, dq_tgt_mean, z_km; label = "Target", linewidth = 2)
    lines!(ax2, dq_prd_mean, z_km; label = "Predicted", linewidth = 2, linestyle = :dash)
    axislegend(ax2; position = :rt)

    return fig
end

"""Plot 2: Mean predicted flux profiles F_T(z) and F_q(z)."""
function plot_mean_fluxes(z_faces_km, F_T, F_q)
    fig = Figure(size = (1000, 500))

    ax1 = Axis(fig[1, 1]; xlabel = "F_T flux (K·m/s)", ylabel = "Height (km)",
               title = "Mean temperature flux correction")
    F_T_mean = vec(mean(F_T, dims = 2))
    lines!(ax1, F_T_mean, z_faces_km; linewidth = 2, color = :firebrick)
    vlines!(ax1, [0]; color = :gray, linestyle = :dot)

    ax2 = Axis(fig[1, 2]; xlabel = "F_q flux (g/kg·m/s)", ylabel = "Height (km)",
               title = "Mean moisture flux correction")
    F_q_mean = vec(mean(F_q, dims = 2)) .* 1000
    lines!(ax2, F_q_mean, z_faces_km; linewidth = 2, color = :steelblue)
    vlines!(ax2, [0]; color = :gray, linestyle = :dot)

    return fig
end

"""Plot 3: Example individual column profiles (target vs predicted)."""
function plot_example_columns(z_km, dT_target, dT_pred, dq_target, dq_pred, n_examples)
    rng = Random.MersenneTwister(99)
    n_samples = size(dT_target, 2)
    idx = randperm(rng, n_samples)[1:min(n_examples, n_samples)]

    ncols = min(4, length(idx))
    nrows = cld(length(idx), ncols)
    fig = Figure(size = (300 * ncols, 350 * nrows))

    for (panel, ci) in enumerate(idx)
        r = cld(panel, ncols)
        c = mod1(panel, ncols)
        ax = Axis(fig[r, c]; xlabel = "dT/dt (K/hr)", ylabel = "Height (km)",
                  title = "Column $ci")
        lines!(ax, dT_target[:, ci] .* 3600, z_km; label = "Target", linewidth = 1.5)
        lines!(ax, dT_pred[:, ci] .* 3600, z_km; label = "Predicted",
               linewidth = 1.5, linestyle = :dash)
        if panel == 1
            axislegend(ax; position = :rt, labelsize = 10)
        end
    end

    Label(fig[0, :], "Temperature tendency: individual columns"; fontsize = 16)
    return fig
end

"""Plot 4: Scatter plot of predicted vs target tendencies."""
function plot_scatter(dT_target, dT_pred, dq_target, dq_pred)
    fig = Figure(size = (1000, 500))

    # Subsample for plotting speed
    n = size(dT_target, 2) * size(dT_target, 1)
    max_pts = 50_000
    stride = max(1, n ÷ max_pts)

    ax1 = Axis(fig[1, 1]; xlabel = "Target dT/dt (K/hr)", ylabel = "Predicted dT/dt (K/hr)",
               title = "Temperature tendency", aspect = 1)
    tgt = vec(dT_target)[1:stride:end] .* 3600
    prd = vec(dT_pred)[1:stride:end] .* 3600
    scatter!(ax1, tgt, prd; markersize = 1, color = (:firebrick, 0.15))
    lims = extrema(vcat(tgt, prd))
    lines!(ax1, [lims...], [lims...]; color = :black, linewidth = 1)

    ax2 = Axis(fig[1, 2]; xlabel = "Target dq/dt (g/kg/hr)", ylabel = "Predicted dq/dt (g/kg/hr)",
               title = "Moisture tendency", aspect = 1)
    tgt_q = vec(dq_target)[1:stride:end] .* 3600 .* 1000
    prd_q = vec(dq_pred)[1:stride:end] .* 3600 .* 1000
    scatter!(ax2, tgt_q, prd_q; markersize = 1, color = (:steelblue, 0.15))
    lims_q = extrema(vcat(tgt_q, prd_q))
    lines!(ax2, [lims_q...], [lims_q...]; color = :black, linewidth = 1)

    return fig
end

"""Plot 5: Vertical profile of RMSE and bias."""
function plot_error_profiles(z_km, dT_target, dT_pred, dq_target, dq_pred)
    fig = Figure(size = (1000, 500))

    dT_err = dT_pred .- dT_target
    dq_err = dq_pred .- dq_target

    dT_rmse = sqrt.(vec(mean(dT_err .^ 2, dims = 2))) .* 3600
    dT_bias = vec(mean(dT_err, dims = 2)) .* 3600
    dq_rmse = sqrt.(vec(mean(dq_err .^ 2, dims = 2))) .* 3600 .* 1000
    dq_bias = vec(mean(dq_err, dims = 2)) .* 3600 .* 1000

    ax1 = Axis(fig[1, 1]; xlabel = "K/hr", ylabel = "Height (km)",
               title = "Temperature error profile")
    lines!(ax1, dT_rmse, z_km; label = "RMSE", linewidth = 2, color = :firebrick)
    lines!(ax1, dT_bias, z_km; label = "Bias", linewidth = 2, color = :firebrick,
           linestyle = :dash)
    vlines!(ax1, [0]; color = :gray, linestyle = :dot)
    axislegend(ax1; position = :rt)

    ax2 = Axis(fig[1, 2]; xlabel = "g/kg/hr", ylabel = "Height (km)",
               title = "Moisture error profile")
    lines!(ax2, dq_rmse, z_km; label = "RMSE", linewidth = 2, color = :steelblue)
    lines!(ax2, dq_bias, z_km; label = "Bias", linewidth = 2, color = :steelblue,
           linestyle = :dash)
    vlines!(ax2, [0]; color = :gray, linestyle = :dot)
    axislegend(ax2; position = :rt)

    return fig
end

"""Plot 6: R² skill score as a function of height."""
function plot_r2_profile(z_km, dT_target, dT_pred, dq_target, dq_pred)
    fig = Figure(size = (600, 500))

    r2_T = zeros(length(z_km))
    r2_q = zeros(length(z_km))
    for k in eachindex(z_km)
        ss_res_T = sum((dT_pred[k, :] .- dT_target[k, :]) .^ 2)
        ss_tot_T = sum((dT_target[k, :] .- mean(dT_target[k, :])) .^ 2)
        r2_T[k] = ss_tot_T > 0 ? 1 - ss_res_T / ss_tot_T : 0.0

        ss_res_q = sum((dq_pred[k, :] .- dq_target[k, :]) .^ 2)
        ss_tot_q = sum((dq_target[k, :] .- mean(dq_target[k, :])) .^ 2)
        r2_q[k] = ss_tot_q > 0 ? 1 - ss_res_q / ss_tot_q : 0.0
    end

    ax = Axis(fig[1, 1]; xlabel = "R² score", ylabel = "Height (km)",
              title = "Skill score by level")
    lines!(ax, r2_T, z_km; label = "Temperature", linewidth = 2, color = :firebrick)
    lines!(ax, r2_q, z_km; label = "Moisture", linewidth = 2, color = :steelblue)
    vlines!(ax, [0]; color = :gray, linestyle = :dot)
    vlines!(ax, [1]; color = :gray, linestyle = :dot)
    xlims!(ax, -0.5, 1.05)
    axislegend(ax; position = :lb)

    return fig
end

"""Plot 7: Magnitude of correction vs model/ERA5 tendencies."""
function plot_tendency_magnitudes(z_km, dT_pred, dq_pred, Y, nz)
    fig = Figure(size = (1000, 500))

    target_T = Y[1:nz, :]
    target_q = Y[nz+1:2*nz, :]

    # RMS of target tendency (= ERA5 - model) gives the correction magnitude
    rms_correction_T = sqrt.(vec(mean(target_T .^ 2, dims = 2))) .* 3600
    rms_pred_T = sqrt.(vec(mean(dT_pred .^ 2, dims = 2))) .* 3600
    rms_residual_T = sqrt.(vec(mean((dT_pred .- target_T) .^ 2, dims = 2))) .* 3600

    ax1 = Axis(fig[1, 1]; xlabel = "RMS tendency (K/hr)", ylabel = "Height (km)",
               title = "Temperature: correction magnitude")
    lines!(ax1, rms_correction_T, z_km; label = "Target correction", linewidth = 2)
    lines!(ax1, rms_pred_T, z_km; label = "Predicted correction", linewidth = 2, linestyle = :dash)
    lines!(ax1, rms_residual_T, z_km; label = "Residual error", linewidth = 2, linestyle = :dot)
    axislegend(ax1; position = :rt)

    rms_correction_q = sqrt.(vec(mean(target_q .^ 2, dims = 2))) .* 3600 .* 1000
    rms_pred_q = sqrt.(vec(mean(dq_pred .^ 2, dims = 2))) .* 3600 .* 1000
    rms_residual_q = sqrt.(vec(mean((dq_pred .- target_q) .^ 2, dims = 2))) .* 3600 .* 1000

    ax2 = Axis(fig[1, 2]; xlabel = "RMS tendency (g/kg/hr)", ylabel = "Height (km)",
               title = "Moisture: correction magnitude")
    lines!(ax2, rms_correction_q, z_km; label = "Target correction", linewidth = 2)
    lines!(ax2, rms_pred_q, z_km; label = "Predicted correction", linewidth = 2, linestyle = :dash)
    lines!(ax2, rms_residual_q, z_km; label = "Residual error", linewidth = 2, linestyle = :dot)
    axislegend(ax2; position = :rt)

    return fig
end

# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

function evaluate()
    mkpath(PLOT_DIR)

    model, ps, st, X_norm_stats, nz = load_trained_model(MODEL_PATH)
    global X_norm_stats_global = X_norm_stats

    println("\nBuilding evaluation dataset...")
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
    dT_pred, dq_pred, F_T, F_q = predict_tendencies(model, ps, st, X_raw, dz_all, nz)

    dT_target = Y_raw[1:nz, :]
    dq_target = Y_raw[nz+1:2*nz, :]

    # Overall metrics
    rmse_T = sqrt(mean((dT_pred .- dT_target) .^ 2)) * 3600
    rmse_q = sqrt(mean((dq_pred .- dq_target) .^ 2)) * 3600 * 1000
    corr_T = cor(vec(dT_pred), vec(dT_target))
    corr_q = cor(vec(dq_pred), vec(dq_target))
    @printf("  RMSE  T: %.4f K/hr    q: %.4f g/kg/hr\n", rmse_T, rmse_q)
    @printf("  Corr  T: %.4f         q: %.4f\n", corr_T, corr_q)

    z_km = get_z_centres(nz)
    z_faces_km = get_z_faces(z_km)

    println("\nGenerating plots...")

    fig1 = plot_mean_profiles(z_km, dT_target, dT_pred, dq_target, dq_pred)
    save(joinpath(PLOT_DIR, "01_mean_tendency_profiles.png"), fig1; px_per_unit = 2)
    println("  01_mean_tendency_profiles.png")

    fig2 = plot_mean_fluxes(z_faces_km, F_T, F_q)
    save(joinpath(PLOT_DIR, "02_mean_flux_profiles.png"), fig2; px_per_unit = 2)
    println("  02_mean_flux_profiles.png")

    fig3 = plot_example_columns(z_km, dT_target, dT_pred, dq_target, dq_pred, N_EXAMPLE_COLS)
    save(joinpath(PLOT_DIR, "03_example_columns.png"), fig3; px_per_unit = 2)
    println("  03_example_columns.png")

    fig4 = plot_scatter(dT_target, dT_pred, dq_target, dq_pred)
    save(joinpath(PLOT_DIR, "04_scatter_pred_vs_target.png"), fig4; px_per_unit = 2)
    println("  04_scatter_pred_vs_target.png")

    fig5 = plot_error_profiles(z_km, dT_target, dT_pred, dq_target, dq_pred)
    save(joinpath(PLOT_DIR, "05_error_profiles.png"), fig5; px_per_unit = 2)
    println("  05_error_profiles.png")

    fig6 = plot_r2_profile(z_km, dT_target, dT_pred, dq_target, dq_pred)
    save(joinpath(PLOT_DIR, "06_r2_skill_score.png"), fig6; px_per_unit = 2)
    println("  06_r2_skill_score.png")

    fig7 = plot_tendency_magnitudes(z_km, dT_pred, dq_pred, Y_raw, nz)
    save(joinpath(PLOT_DIR, "07_tendency_magnitudes.png"), fig7; px_per_unit = 2)
    println("  07_tendency_magnitudes.png")

    println("\nAll plots saved to $(PLOT_DIR)/")
    println("Done.")
end

if abspath(PROGRAM_FILE) == @__FILE__
    evaluate()
end
