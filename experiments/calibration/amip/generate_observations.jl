# Prevent MPI from being used when generating observations.
ENV["CLIMACOMMS_CONTEXT"] = "SINGLETON"

using Dates
import Statistics
import ClimaAnalysis
import ClimaAnalysis: OutputVar
import ClimaCalibrate
import ClimaCalibrate: ObservationRecipe, SampleBuilder
import ClimaCoupler
import ClimaCoupler: CalibrationTools
import EnsembleKalmanProcesses as EKP
import EnsembleKalmanProcesses.ParameterDistributions as PD
import JLD2

include(joinpath(@__DIR__, "preprocessing.jl"))

"""
    make_observation_vector(
        var,
        sample_date_ranges;
        covariance_date_ranges = sample_date_ranges,
        beta = 0.05,
        rank = 2,
        sigma = 1e-3,
        use_latitude_weights = true,
        min_cosd_lat = 0.1,
        noise_sd = nothing,
        effective_entries = nothing,
        FT = Float32,
    )

Make a vector of `EKP.Observation`s for `var`, one for each sample corresponding to the
dates in `sample_date_ranges`. When `noise_sd` is given (for a one-value-per-date reduction
such as a global mean), the observation has independent noise of standard deviation
`noise_sd`. Otherwise it uses an `SVDplusD` covariance matrix estimated across
`covariance_date_ranges`.

The `SVDplusD` covariance is `ObservationRecipe.SVDplusModelErrorCovariance`: a truncated
SVD of the sample covariance with `rank` modes, plus the sample variance those modes leave
out, plus a model error term with relative scale `beta`. When the sample mean of `var`
crosses zero (such as an annual pattern or a seasonal anomaly), `model_error_covariance`
floors the squared mean at its spatial median so the model error does not vanish along zero
contours.

With `use_latitude_weights`, the covariance is estimated from samples weighted by
`SampleBuilder.LatitudeWeighting`, which scales the variance at each point by
`1 / max(cosd(lat), min_cosd_lat)`. When `effective_entries` is set, `covariance_scale`
multiplies the covariance by `n_entries / effective_entries`:

1. Multiplying by `n_entries` normalizes the block (whether 1 scalar or a multi-season map)
   to its average squared whitened error per entry (`(r' * inv(Gamma) * r) / n_entries`).
2. Dividing by `effective_entries` gives that per-entry average weight `effective_entries`
   in the loss and adds `effective_entries` to the total observation count in
   `EffectiveEntriesController`.

The observation itself is unweighted.

Each date range is shrunk to the dates `var` has inside it (`snap_date_ranges`).
"""
function covariance_scale(samples, effective_entries)
    isnothing(effective_entries) && return 1.0
    n_entries = size(SampleBuilder.get_samples(samples), 1)
    return n_entries / effective_entries
end

function make_observation_vector(
    var,
    sample_date_ranges;
    covariance_date_ranges = sample_date_ranges,
    beta = 0.05,
    rank = 2,
    sigma = 1e-3,
    use_latitude_weights = true,
    min_cosd_lat = 0.1,
    noise_sd = nothing,
    effective_entries = nothing,
    FT = Float32,
)
    sample_date_ranges = snap_date_ranges(sample_date_ranges, var)
    if !isnothing(noise_sd)
        return map(sample_date_ranges) do date_range
            samples = SampleBuilder.build_samples_by_times([var], [date_range]; FT)
            scale = covariance_scale(samples, effective_entries)
            cov = ObservationRecipe.DiagonalCovariance(
                ObservationRecipe.ScalarDiagonal(scale * noise_sd^2),
            )
            ObservationRecipe.observation(cov, samples, 1)
        end
    end
    covariance_date_ranges = snap_date_ranges(covariance_date_ranges, var)
    # The targets are samples of this collection, so every observation's covariance is
    # estimated from all of its samples.
    samples = SampleBuilder.build_samples_by_times([var], covariance_date_ranges; FT)
    @info "Built covariance samples" samples
    sample_mean = Statistics.mean(SampleBuilder.get_samples(samples); dims = 2)
    floor_mean = any(<(0), sample_mean) && any(>(0), sample_mean)
    @info "Using SVDplusD covariance matrix with" name = ClimaAnalysis.short_name(var) beta rank sigma floor_mean use_latitude_weights min_cosd_lat effective_entries
    covar_estimator = model_error_covariance(; rank, beta, sigma, floor_mean)

    weighted_samples = samples
    use_latitude_weights &&
        ClimaAnalysis.has_latitude(var) && (
            weighted_samples =
                weighted_samples |> SampleBuilder.LatitudeWeighting(; min_cosd_lat)
        )
    !isnothing(effective_entries) && (
        weighted_samples =
            weighted_samples |> SampleBuilder.PerCollectionWeighting(
                sqrt(covariance_scale(samples, effective_entries)),
            )
    )
    return map(sample_date_ranges) do date_range
        i = findfirst(==(date_range), covariance_date_ranges)
        isnothing(i) && error(
            "The target $date_range is not one of the covariance_date_ranges, which the \
            observation is taken from",
        )
        ObservationRecipe.observation(covar_estimator, weighted_samples, i)
    end
end

"""
    model_error_covariance(; rank, beta, sigma, floor_mean = false)

`ObservationRecipe.SVDplusModelErrorCovariance(; rank, beta, sigma)`, with the squared mean
in its model-error term raised to at least its median when `floor_mean` is true. An anomaly
map's mean crosses zero, so without the floor its model error vanishes there and those
points dominate the loss.
"""
function model_error_covariance(; rank, beta, sigma, floor_mean = false)
    floor_mean ||
        return ObservationRecipe.SVDplusModelErrorCovariance(; rank, beta, sigma)
    mean_squared = ObservationRecipe.ModelErrorScaleDiagonal(1.0)
    floor = max.(
        sigma^2 .*
        ObservationRecipe.QuantileDiagonal(0.5, ObservationRecipe.VarianceDiagonal()),
        ObservationRecipe.QuantileDiagonal(0.5, mean_squared),
    )
    diagonal =
        ObservationRecipe.ResidualVarianceDiagonal(rank) .+
        ObservationRecipe.ScalarDiagonal(float.(vcat(beta))) .^ 2 .*
        max.(mean_squared, floor)
    return ObservationRecipe.SVDplusDCovariance(diagonal; rank)
end

"""
    build_observation_vector(vars, config)

One `EKP.Observation` per target date, with one noise covariance block per observable
(`make_observation_vector`), scaled by each observable's `effective_entries` when set.
Building one block per observable gives each field and reduction its own `rank` leading
interannual SVD modes (instead of forcing all observables to share `rank` global modes) and
lets observables with different post-reduction time axes (`annual_average` versus
`seasonal_anomaly`) be combined in one `Observation`. `vars` are the reduced values, one per
observable, in the order of `config.observables`.
"""
function build_observation_vector(vars, config)
    per_observable = map(vars, config.observables) do var, o
        make_observation_vector(
            var,
            config.sample_date_ranges;
            covariance_date_ranges = config.covariance_date_ranges,
            beta = something(o.model_error_fraction, 0.0),
            noise_sd = o.noise_sd,
            effective_entries = o.effective_entries,
        )
    end
    return [
        EKP.combine_observations([obs[i] for obs in per_observable]) for
        i in eachindex(first(per_observable))
    ]
end

"""
    generate_observations(config; data_loader = CalibrationTools.default_data_loader())

Build and save the coverage masks and observation vector for `config` to
`config.output_dir`.
"""
function generate_observations(
    config;
    data_loader = CalibrationTools.default_data_loader(),
)
    (; observables, sample_date_ranges, covariance_date_ranges, output_dir) = config
    length(unique(covariance_date_ranges)) >= 3 ||
        error("The config needs at least 3 distinct covariance_date_ranges")
    mkpath(output_dir)

    vars = map(unique(o.field for o in observables)) do field
        source_data_loader = CalibrationTools.find_source_loader(data_loader, field)
        @info "Retrieving $field from $(typeof(source_data_loader))"
        return get(source_data_loader, field)
    end

    vars = select_pressure_levels.(vars, Ref(config.pressure_levels))
    lonlat_regridder = get_lonlat_regridder(config.config_file)
    vars = lonlat_regridder.(vars)

    # Record where each observation has data over the dates in use, while longitude still
    # exists. A globally covered product gets no mask. `observation_map.jl` loads the same
    # masks, and `mask_and_reduce` applies them on both sides.
    # A target is dated by its first month, and its reduced value is made from the months
    # up to the end of `extend`, so the mask covers all of them.
    mask_date_ranges = unique([
        (first(r), last(r) + config.extend - Dates.Day(1)) for
        r in vcat(sample_date_ranges, covariance_date_ranges)
    ])
    coverage_masks = Dict{String, Any}()
    for var in vars
        mask = coverage_mask(var, mask_date_ranges)
        isnothing(mask) || (coverage_masks[ClimaAnalysis.short_name(var)] = mask)
    end

    vars = mask_and_reduce(vars, coverage_masks, config)
    check_extend(vars, config)

    # Give every sample the same NaN mask, or `build_samples_by_times` rejects a product
    # whose coverage varies by year.
    vars = ClimaAnalysis.propagate_nans.(vars; dims = ("time",))

    # `observation_map.jl` reads these so the simulation is sampled at the same points as
    # the observation.
    JLD2.save_object(coverage_mask_path(output_dir), coverage_masks)

    # Create the observation vector, in physical units, with one noise covariance per
    # observable.
    observation_vec = build_observation_vector(vars, config)

    JLD2.save_object(observation_vector_path(output_dir), observation_vec)

    for (i, obs) in enumerate(observation_vec)
        @info "Observation $i"
        @info ObservationRecipe.reconstruct_vars(obs)
    end
    return observation_vec
end

if abspath(PROGRAM_FILE) == @__FILE__
    if !isdefined(Main, :CALIBRATE_CONFIG)
        config_dir =
            joinpath(pkgdir(ClimaCoupler), "experiments", "calibration", "amip", "config")
        default_config_path = joinpath(config_dir, "year2010_rsut_rlut_pr.jl")
        test_calibration_config_path = joinpath(config_dir, "pipeline_test.jl")
        test_calibration = haskey(ENV, "TEST_CALIBRATION")
        config_path =
            test_calibration ? test_calibration_config_path :
            get(ENV, "CALIBRATION_CONFIG", default_config_path)
        @info "Using calibration configuration in: $config_path"
        include(config_path)
    end
    generate_observations(CALIBRATE_CONFIG)
end
