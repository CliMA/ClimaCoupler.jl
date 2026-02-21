#!/usr/bin/env julia
#=
analyze_calibration.jl - Post-calibration parameter analysis and export

Analyzes calibration parameters at a given iteration by reading directly from 
ensemble member output directories (ground truth parameters.toml files).

Features:
- Parameter statistics: mean, std, min, max, spread, quartiles
- Extremes detection and ensemble spread analysis
- Multiple selection strategies for final parameters:
  - mean_nearest: ensemble member nearest to the mean (default)
  - lowest_error: member with lowest EKP error
  - median_nearest: member nearest to the median
  - ensemble_mean: use actual ensemble mean values (not from a specific member)
- Merges calibrated parameters with base TOML overrides

Usage:
  julia --project=experiments/ClimaEarth experiments/calibration/subseasonal/analyze_calibration.jl

Configuration is at the top of the script.
=#

using Printf
using Statistics
using LinearAlgebra
import TOML
import YAML
import JLD2
import EnsembleKalmanProcesses as EKP
import EnsembleKalmanProcesses.ParameterDistributions as PD
import ClimaCoupler

# Override JLD2's default_iotype to avoid Lustre issues
JLD2.default_iotype() = IOStream

# =============================================================================
# CONFIGURATION - Edit these as needed
# =============================================================================
# Calibration output directory (where iteration_XXX folders are)
const CALIBRATION_OUTPUT_DIR = "/glade/derecho/scratch/zhaoyi/calibration/"

# Iteration to analyze (0-indexed). Set to -1 to use the latest available iteration.
const ITERATION_TO_ANALYZE = 1

# Config file used for calibration (yml). The base TOML is inferred from this.
# This should match the config_file used in run_calibration.jl's CALIBRATE_CONFIG
const CALIBRATION_CONFIG_FILE = joinpath(
    pkgdir(ClimaCoupler),
    "config/subseasonal_configs/wxquest_diagedmf.yml"
)

# Base output directory for generated files (created in current working directory)
# Subfolders are auto-created: calibration_analysis/<exp_name>/iter<XXX>/
const LOCAL_OUTPUT_BASE = "calibration_analysis"

# Selection strategy for default output:
# Options: "mean_nearest", "lowest_error", "median_nearest", "ensemble_mean"
const DEFAULT_SELECTION_STRATEGY = "mean_nearest"

# =============================================================================
# DATA STRUCTURES
# =============================================================================
struct ParameterStats
    names::Vector{String}
    values::Matrix{Float64}  # (n_params, n_members) - raw values from each member
    mean::Vector{Float64}
    std::Vector{Float64}
    min::Vector{Float64}
    max::Vector{Float64}
    median::Vector{Float64}
    q25::Vector{Float64}
    q75::Vector{Float64}
    spread::Vector{Float64}  # max - min
    cv::Vector{Float64}      # coefficient of variation (std/mean)
end

struct MemberError
    member_id::Int
    error::Float64
end

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

"""
    get_experiment_name(output_dir)

Extract the experiment name from the calibration output directory.
E.g., "/glade/derecho/scratch/user/calibration/exp25" -> "exp25"
"""
function get_experiment_name(output_dir::String)
    return basename(rstrip(output_dir, '/'))
end

"""
    get_output_dir(base_dir, exp_name, iteration)

Construct the output directory path with experiment name and iteration.
Creates: <base_dir>/<exp_name>/iter<XXX>/
"""
function get_output_dir(base_dir::String, exp_name::String, iteration::Int)
    return joinpath(base_dir, exp_name, @sprintf("iter%03d", iteration))
end

"""
    get_base_toml_paths(config_file)

Extract the base TOML file path(s) from the calibration config file (yml).
Returns a vector of absolute paths to the TOML files.
"""
function get_base_toml_paths(config_file::String)
    if !isfile(config_file)
        error("Config file not found: $config_file")
    end
    
    # Load the YAML config file
    config_dict = YAML.load_file(config_file)
    
    # Get coupler_toml - it may be a single string or a vector
    coupler_toml = get(config_dict, "coupler_toml", String[])
    if coupler_toml isa String
        coupler_toml = [coupler_toml]
    end
    
    if isempty(coupler_toml)
        @warn "No coupler_toml found in config file $config_file"
        return String[]
    end
    
    # Convert to absolute paths (relative to ClimaCoupler package directory)
    toml_paths = String[]
    for toml_file in coupler_toml
        if isfile(toml_file)
            push!(toml_paths, abspath(toml_file))
        else
            # Try relative to ClimaCoupler package
            pkg_path = joinpath(pkgdir(ClimaCoupler), toml_file)
            if isfile(pkg_path)
                push!(toml_paths, pkg_path)
            else
                @warn "TOML file not found: $toml_file (also checked $pkg_path)"
            end
        end
    end
    
    return toml_paths
end

"""
    find_latest_iteration(output_dir)

Find the latest iteration directory number in the output directory.
"""
function find_latest_iteration(output_dir)
    dirs = filter(d -> startswith(d, "iteration_"), readdir(output_dir))
    if isempty(dirs)
        error("No iteration directories found in $output_dir")
    end
    # Extract iteration numbers and find max
    iter_nums = [parse(Int, replace(d, "iteration_" => "")) for d in dirs]
    return maximum(iter_nums)
end

"""
    get_iteration_path(output_dir, iteration)

Get the path to the iteration directory.
"""
function get_iteration_path(output_dir, iteration)
    return joinpath(output_dir, @sprintf("iteration_%03d", iteration))
end

"""
    load_member_parameters(output_dir, iteration, member)

Load parameters.toml from a specific ensemble member directory.
Returns Dict with parameter names and their values.
"""
function load_member_parameters(output_dir, iteration, member)
    member_path = joinpath(
        get_iteration_path(output_dir, iteration),
        @sprintf("member_%03d", member)
    )
    param_file = joinpath(member_path, "parameters.toml")
    
    if !isfile(param_file)
        @warn "Parameter file not found: $param_file"
        return nothing
    end
    
    return TOML.parsefile(param_file)
end

"""
    count_ensemble_members(output_dir, iteration)

Count the number of ensemble members in an iteration directory.
"""
function count_ensemble_members(output_dir, iteration)
    iter_path = get_iteration_path(output_dir, iteration)
    member_dirs = filter(d -> startswith(d, "member_"), readdir(iter_path))
    return length(member_dirs)
end

"""
    extract_param_names_and_values(param_dict)

Extract parameter names and values from a TOML dict.
Returns (names, values) where names is a vector of strings and values is a vector of floats.
"""
function extract_param_names_and_values(param_dict)
    names = String[]
    values = Float64[]
    for (name, entry) in param_dict
        if haskey(entry, "value")
            push!(names, name)
            push!(values, Float64(entry["value"]))
        end
    end
    # Sort by name for consistent ordering
    perm = sortperm(names)
    return names[perm], values[perm]
end

"""
    load_all_member_parameters(output_dir, iteration)

Load parameters from all ensemble members for a given iteration.
Returns a ParameterStats struct with comprehensive statistics.
"""
function load_all_member_parameters(output_dir, iteration)
    n_members = count_ensemble_members(output_dir, iteration)
    @info "Found $n_members ensemble members"
    
    # Load first member to get parameter names
    first_params = load_member_parameters(output_dir, iteration, 1)
    if isnothing(first_params)
        error("Could not load first member parameters")
    end
    param_names, _ = extract_param_names_and_values(first_params)
    n_params = length(param_names)
    
    # Load all members
    values = zeros(n_params, n_members)
    for m in 1:n_members
        member_params = load_member_parameters(output_dir, iteration, m)
        if isnothing(member_params)
            @warn "Skipping member $m - no parameters found"
            values[:, m] .= NaN
            continue
        end
        names, vals = extract_param_names_and_values(member_params)
        
        # Verify parameter names match
        if names != param_names
            @warn "Parameter names mismatch in member $m"
        end
        values[:, m] = vals
    end
    
    # Compute statistics (ignoring NaN)
    mean_vals = [mean(filter(!isnan, values[i, :])) for i in 1:n_params]
    std_vals = [std(filter(!isnan, values[i, :])) for i in 1:n_params]
    min_vals = [minimum(filter(!isnan, values[i, :])) for i in 1:n_params]
    max_vals = [maximum(filter(!isnan, values[i, :])) for i in 1:n_params]
    median_vals = [median(filter(!isnan, values[i, :])) for i in 1:n_params]
    q25_vals = [quantile(filter(!isnan, values[i, :]), 0.25) for i in 1:n_params]
    q75_vals = [quantile(filter(!isnan, values[i, :]), 0.75) for i in 1:n_params]
    spread = max_vals .- min_vals
    cv = std_vals ./ abs.(mean_vals)  # Coefficient of variation
    
    return ParameterStats(
        param_names, values, mean_vals, std_vals, min_vals, max_vals,
        median_vals, q25_vals, q75_vals, spread, cv
    )
end

"""
    load_ekp_data(output_dir, iteration)

Load the EKP object and compute per-member errors.
Returns vector of MemberError structs sorted by error.
"""
function load_ekp_data(output_dir, iteration)
    ekp_path = joinpath(get_iteration_path(output_dir, iteration), "eki_file.jld2")
    g_path = joinpath(get_iteration_path(output_dir, iteration), "G_ensemble.jld2")
    
    if !isfile(ekp_path)
        @warn "EKP file not found: $ekp_path"
        return nothing
    end
    
    ekp = JLD2.load(ekp_path)["single_stored_object"]
    
    # Load G ensemble if available for per-member error computation
    member_errors = MemberError[]
    if isfile(g_path)
        G = JLD2.load(g_path)["single_stored_object"]
        
        # Get observations
        y = EKP.get_obs(ekp)
        
        # Compute per-member squared error (using all iterations up to this one)
        n_members = size(G, 2)
        for m in 1:n_members
            g_m = G[:, m]
            if any(isnan, g_m)
                push!(member_errors, MemberError(m, Inf))
            else
                mse = mean((g_m .- y).^2)
                push!(member_errors, MemberError(m, sqrt(mse)))  # RMSE
            end
        end
    else
        @warn "G_ensemble.jld2 not found - cannot compute per-member errors"
    end
    
    return sort(member_errors, by=x -> x.error)
end

"""
    load_prior(output_dir)

Load the prior distribution from iteration_000.
"""
function load_prior(output_dir)
    prior_path = joinpath(get_iteration_path(output_dir, 0), "prior.jld2")
    if !isfile(prior_path)
        @warn "Prior file not found: $prior_path"
        return nothing
    end
    return JLD2.load(prior_path)["single_stored_object"]
end

# =============================================================================
# ANALYSIS & REPORTING
# =============================================================================

"""
    print_parameter_statistics(stats::ParameterStats, member_errors=nothing)

Print a formatted table of parameter statistics.
"""
function print_parameter_statistics(stats::ParameterStats, member_errors=nothing)
    n_params = length(stats.names)
    n_members = size(stats.values, 2)
    
    println("\n" * "="^110)
    println("PARAMETER STATISTICS")
    println("="^110)
    println(@sprintf("%-45s %12s %12s %12s %12s %10s", 
                     "Parameter", "Mean", "Std", "Min", "Max", "CV%"))
    println("-"^110)
    
    for i in 1:n_params
        cv_pct = stats.cv[i] * 100
        println(@sprintf("%-45s %12.6g %12.6g %12.6g %12.6g %10.1f",
                         stats.names[i], stats.mean[i], stats.std[i],
                         stats.min[i], stats.max[i], cv_pct))
    end
    println("="^110)
    
    # Print quartile information
    println("\nQUARTILE SUMMARY")
    println("-"^90)
    println(@sprintf("%-45s %12s %12s %12s", "Parameter", "Q25", "Median", "Q75"))
    println("-"^90)
    for i in 1:n_params
        println(@sprintf("%-45s %12.6g %12.6g %12.6g",
                         stats.names[i], stats.q25[i], stats.median[i], stats.q75[i]))
    end
    println("-"^90)
    
    # Print extreme values (members with min/max for each parameter)
    println("\nEXTREME VALUES BY MEMBER")
    println("-"^80)
    for i in 1:n_params
        valid_vals = [(m, stats.values[i, m]) for m in 1:n_members if !isnan(stats.values[i, m])]
        if isempty(valid_vals)
            continue
        end
        min_member = argmin(v -> v[2], valid_vals)[1]
        max_member = argmax(v -> v[2], valid_vals)[1]
        println(@sprintf("%-45s min=member_%03d (%.6g), max=member_%03d (%.6g)",
                         stats.names[i], min_member, stats.min[i], max_member, stats.max[i]))
    end
    
    # Print per-member errors if available
    if !isnothing(member_errors) && !isempty(member_errors)
        println("\n" * "="^60)
        println("MEMBER ERRORS (sorted by RMSE)")
        println("="^60)
        for (rank, me) in enumerate(member_errors[1:min(10, length(member_errors))])
            println(@sprintf("  Rank %2d: member_%03d  RMSE = %.6g", rank, me.member_id, me.error))
        end
        if length(member_errors) > 10
            println("  ... ($(length(member_errors) - 10) more members)")
        end
    end
end

"""
    select_member(stats::ParameterStats, strategy::String, member_errors=nothing)

Select an ensemble member based on the given strategy.
Returns (member_id, parameter_values).
"""
function select_member(stats::ParameterStats, strategy::String, member_errors=nothing)
    n_members = size(stats.values, 2)
    n_params = length(stats.names)
    
    if strategy == "mean_nearest"
        # Find member nearest to ensemble mean (L2 distance in normalized space)
        distances = Float64[]
        for m in 1:n_members
            if any(isnan, stats.values[:, m])
                push!(distances, Inf)
            else
                # Normalize by std to give equal weight to each parameter
                normed_diff = (stats.values[:, m] .- stats.mean) ./ max.(stats.std, 1e-10)
                push!(distances, norm(normed_diff))
            end
        end
        member_id = argmin(distances)
        return member_id, stats.values[:, member_id]
        
    elseif strategy == "median_nearest"
        # Find member nearest to median values
        distances = Float64[]
        for m in 1:n_members
            if any(isnan, stats.values[:, m])
                push!(distances, Inf)
            else
                normed_diff = (stats.values[:, m] .- stats.median) ./ max.(stats.std, 1e-10)
                push!(distances, norm(normed_diff))
            end
        end
        member_id = argmin(distances)
        return member_id, stats.values[:, member_id]
        
    elseif strategy == "lowest_error"
        if isnothing(member_errors) || isempty(member_errors)
            error("Cannot use 'lowest_error' strategy: no error data available")
        end
        member_id = first(member_errors).member_id
        return member_id, stats.values[:, member_id]
        
    elseif strategy == "ensemble_mean"
        # Return the ensemble mean (not from a specific member)
        return 0, stats.mean
        
    else
        error("Unknown selection strategy: $strategy")
    end
end

"""
    merge_with_base_tomls(calibrated_params::Dict, base_toml_paths::Vector{String})

Merge calibrated parameters with base TOML file(s).
Calibrated parameters override base parameters where they overlap.
Non-calibrated parameters from base TOML files are preserved.

When multiple TOML files are provided, they are merged in order (later files override earlier ones),
and then calibrated parameters override everything.
"""
function merge_with_base_tomls(calibrated_params::Dict{String,Any}, base_toml_paths::Vector{String})
    merged = Dict{String,Any}()
    
    # Load and merge base TOML files in order
    for base_toml_path in base_toml_paths
        if isfile(base_toml_path)
            base_params = TOML.parsefile(base_toml_path)
            # Merge into existing (later files override)
            for (name, entry) in base_params
                merged[name] = deepcopy(entry)
            end
            @info "Loaded $(length(base_params)) parameters from $(basename(base_toml_path))"
        else
            @warn "Base TOML not found: $base_toml_path"
        end
    end
    
    if isempty(base_toml_paths)
        @warn "No base TOML files provided - output will only contain calibrated parameters"
    end
    
    # Override with calibrated parameters
    n_overwritten = 0
    n_added = 0
    for (name, entry) in calibrated_params
        if haskey(merged, name)
            n_overwritten += 1
        else
            n_added += 1
        end
        merged[name] = entry
    end
    
    @info "Merged: $n_overwritten overwritten, $n_added added (total: $(length(merged)) parameters)"
    return merged
end

"""
    create_toml_dict(param_names::Vector{String}, param_values::Vector{Float64})

Create a TOML-compatible dictionary from parameter names and values.
"""
function create_toml_dict(param_names::Vector{String}, param_values::Vector{Float64})
    result = Dict{String,Any}()
    for (name, value) in zip(param_names, param_values)
        result[name] = Dict{String,Any}(
            "value" => value,
            "type" => "float"
        )
    end
    return result
end

"""
    save_parameters_toml(params::Dict, filepath::String)

Save parameters to a TOML file.
"""
function save_parameters_toml(params::Dict, filepath::String)
    open(filepath, "w") do f
        TOML.print(f, params)
    end
    @info "Saved parameters to: $filepath"
end

"""
    save_statistics_csv(stats::ParameterStats, filepath::String)

Save statistics to a CSV file.
"""
function save_statistics_csv(stats::ParameterStats, filepath::String)
    open(filepath, "w") do f
        # Header
        write(f, "parameter,mean,std,min,max,median,q25,q75,spread,cv\n")
        for i in eachindex(stats.names)
            write(f, "$(stats.names[i]),")
            write(f, "$(stats.mean[i]),")
            write(f, "$(stats.std[i]),")
            write(f, "$(stats.min[i]),")
            write(f, "$(stats.max[i]),")
            write(f, "$(stats.median[i]),")
            write(f, "$(stats.q25[i]),")
            write(f, "$(stats.q75[i]),")
            write(f, "$(stats.spread[i]),")
            write(f, "$(stats.cv[i])\n")
        end
    end
    @info "Saved statistics CSV to: $filepath"
end

"""
    save_member_values_csv(stats::ParameterStats, filepath::String)

Save raw parameter values for each member to a CSV file.
"""
function save_member_values_csv(stats::ParameterStats, filepath::String)
    n_members = size(stats.values, 2)
    open(filepath, "w") do f
        # Header
        write(f, "member,")
        write(f, join(stats.names, ","))
        write(f, "\n")
        # Data rows
        for m in 1:n_members
            write(f, @sprintf("member_%03d,", m))
            write(f, join([@sprintf("%.10g", stats.values[i, m]) for i in eachindex(stats.names)], ","))
            write(f, "\n")
        end
    end
    @info "Saved member values CSV to: $filepath"
end

# =============================================================================
# MAIN
# =============================================================================

function main()
    # Determine iteration to analyze
    iteration = ITERATION_TO_ANALYZE
    if iteration < 0
        iteration = find_latest_iteration(CALIBRATION_OUTPUT_DIR)
        @info "Auto-detected latest iteration: $iteration"
    end
    
    # Get experiment name and construct output directory
    exp_name = get_experiment_name(CALIBRATION_OUTPUT_DIR)
    output_dir = get_output_dir(LOCAL_OUTPUT_BASE, exp_name, iteration)
    
    # Get base TOML paths from config file
    base_toml_paths = get_base_toml_paths(CALIBRATION_CONFIG_FILE)
    
    @info "="^70
    @info "CALIBRATION ANALYSIS"
    @info "="^70
    @info "Experiment: $exp_name"
    @info "Calibration dir: $CALIBRATION_OUTPUT_DIR"
    @info "Analyzing iteration: $iteration"
    @info "Config file: $(basename(CALIBRATION_CONFIG_FILE))"
    @info "Base TOML(s): $(join(basename.(base_toml_paths), ", "))"
    @info "Output dir: $output_dir"
    
    # Create output directory
    mkpath(output_dir)
    
    # Load parameters from all members
    @info "\nLoading parameters from ensemble members..."
    stats = load_all_member_parameters(CALIBRATION_OUTPUT_DIR, iteration)
    
    # Load EKP data for error information
    @info "Loading EKP data..."
    member_errors = load_ekp_data(CALIBRATION_OUTPUT_DIR, iteration)
    
    # Load prior for reference
    prior = load_prior(CALIBRATION_OUTPUT_DIR)
    if !isnothing(prior)
        @info "Prior loaded with $(PD.ndims(prior)) parameters"
    end
    
    # Print statistics
    print_parameter_statistics(stats, member_errors)
    
    # Save statistics files
    save_statistics_csv(stats, joinpath(output_dir, "statistics.csv"))
    save_member_values_csv(stats, joinpath(output_dir, "member_values.csv"))
    
    # Generate final parameter files for all strategies
    strategies = ["mean_nearest", "median_nearest", "lowest_error", "ensemble_mean"]
    
    println("\n" * "="^70)
    println("GENERATING FINAL PARAMETER FILES")
    println("="^70)
    
    for strategy in strategies
        if strategy == "lowest_error" && (isnothing(member_errors) || isempty(member_errors))
            @warn "Skipping '$strategy' strategy: no error data"
            continue
        end
        
        member_id, param_values = select_member(stats, strategy, member_errors)
        
        # Create calibrated params dict
        calibrated_dict = create_toml_dict(stats.names, param_values)
        
        # Merge with base TOML(s)
        merged_dict = merge_with_base_tomls(calibrated_dict, base_toml_paths)
        
        # Save
        output_file = joinpath(output_dir, "parameters_$(strategy).toml")
        save_parameters_toml(merged_dict, output_file)
        
        if member_id == 0
            println("  Strategy '$strategy': using ensemble mean values")
        else
            println("  Strategy '$strategy': selected member_$(lpad(member_id, 3, '0'))")
        end
    end
    
    # Create a "default" symlink or copy for the default strategy
    default_source = joinpath(output_dir, "parameters_$(DEFAULT_SELECTION_STRATEGY).toml")
    default_dest = joinpath(output_dir, "parameters_final.toml")
    if isfile(default_source)
        # Copy instead of symlink for portability
        cp(default_source, default_dest; force=true)
        @info "\nDefault parameters ($(DEFAULT_SELECTION_STRATEGY)) copied to: $default_dest"
    end
    
    # Summary
    println("\n" * "="^70)
    println("ANALYSIS COMPLETE")
    println("="^70)
    println("Output files in: $output_dir/")
    println("  - statistics.csv")
    println("  - member_values.csv")
    for strategy in strategies
        if strategy == "lowest_error" && (isnothing(member_errors) || isempty(member_errors))
            continue
        end
        println("  - parameters_$(strategy).toml")
    end
    println("  - parameters_final.toml (default = $(DEFAULT_SELECTION_STRATEGY))")
end

# Run if executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
