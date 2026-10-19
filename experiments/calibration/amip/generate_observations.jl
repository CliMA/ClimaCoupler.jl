import ClimaAnalysis
import ClimaAnalysis: OutputVar
import ClimaCalibrate
import ClimaCalibrate: ObservationRecipe, SampleBuilder
import ClimaCoupler
import ClimaCoupler: CalibrationTools
import JLD2
import LinearAlgebra

include(
    joinpath(
        pkgdir(ClimaCoupler),
        "experiments",
        "calibration",
        "amip",
        "run_calibration.jl",
    ),
)
include(
    joinpath(
        pkgdir(ClimaCoupler),
        "experiments",
        "calibration",
        "amip",
        "preprocessing.jl",
    ),
)
include(
    joinpath(pkgdir(ClimaCoupler), "experiments", "calibration", "amip", "noise_model.jl"),
)

"""
    make_svdplusd_observation_vector(
        vars,
        sample_date_ranges;
        covariance_date_ranges = sample_date_ranges,
        beta = 0.05^2,
        rank = 2,
        sigma2 = 1e-6,
        use_latitude_weights = true,
        min_cosd_lat = 0.1,
        independent_noise_sd = Dict{String, Float64}(),
        FT = Float32,
    )

Make a vector of `EKP.Observation`s with an `SVDplusD` covariance matrix, one for each
sample corresponding to the dates in `sample_date_ranges`.

Variables named in `independent_noise_sd`, one value per target date such as global means,
are left out of the `SVDplusD` covariance and form a second block of the observation, with
independent noise of the given standard deviation. It sets how strongly the calibration holds
such a value, for example the energy balance, against the maps in the first block.

The covariance is the `Gamma` of `noise_model.jl`, estimated from one sample per date
range in `covariance_date_ranges`, so it is the interannual spread of the observation
across those dates. Give it enough date ranges to estimate that spread: `rank` modes need
appreciably more than `rank` samples, and identical date ranges produce a singular
covariance. The calibration targets in `sample_date_ranges` may repeat, and may be a
subset of the covariance dates, since the model only needs initial conditions for those.

`beta` accepts one value for all variables or one per variable, in the order of `vars`
without those in `independent_noise_sd`.
"""
function make_svdplusd_observation_vector(
    vars,
    sample_date_ranges;
    covariance_date_ranges = sample_date_ranges,
    beta = 0.05^2,
    rank = 2,
    sigma2 = 1e-6,
    use_latitude_weights = true,
    min_cosd_lat = 0.1,
    independent_noise_sd = Dict{String, Float64}(),
    FT = Float32,
)
    independent(v) = haskey(independent_noise_sd, ClimaAnalysis.short_name(v))
    independent_vars = filter(independent, vars)
    vars = filter(!independent, vars)
    @info "Using SVDplusD covariance matrix with" beta rank sigma2 use_latitude_weights min_cosd_lat
    covar_estimator =
        noise_covariance_estimator(; beta, rank, sigma2, use_latitude_weights, min_cosd_lat)

    covariance_samples =
        SampleBuilder.build_samples_by_times(vars, covariance_date_ranges; FT)
    @info "Built covariance samples" covariance_samples
    covar = ObservationRecipe.covariance(covar_estimator, covariance_samples)

    target_samples = SampleBuilder.build_samples_by_times(vars, sample_date_ranges; FT)
    @info "Built target samples" target_samples
    n_cov = size(SampleBuilder.get_samples(covariance_samples), 1)
    n_target = size(SampleBuilder.get_samples(target_samples), 1)
    n_cov == n_target || error(
        "The covariance samples ($n_cov values) and the targets ($n_target values) have different lengths",
    )

    # The same assembly as ObservationRecipe.observation, with the covariance supplied.
    obs_vec = map(1:SampleBuilder.num_samples(target_samples)) do i
        sample = collect(view(SampleBuilder.get_samples(target_samples), :, i))
        metadata = collect(view(SampleBuilder.get_metadata(target_samples), :, i))
        name = join(ClimaAnalysis.short_name.(metadata), ";")
        EKP.Observation(
            Dict(
                "samples" => sample,
                "covariances" => covar,
                "names" => name,
                "metadata" => metadata,
            ),
        )
    end
    isempty(independent_vars) && return obs_vec

    # A second block, one value per global mean, with independent noise.
    independent_samples =
        SampleBuilder.build_samples_by_times(independent_vars, sample_date_ranges; FT)
    @info "Global means with independent noise" independent_noise_sd
    return map(enumerate(obs_vec)) do (i, obs)
        metadata = collect(view(SampleBuilder.get_metadata(independent_samples), :, i))
        sample = collect(view(SampleBuilder.get_samples(independent_samples), :, i))
        sd = [FT(independent_noise_sd[ClimaAnalysis.short_name(m)]) for m in metadata]
        length(sd) == length(sample) ||
            error("A global mean must be one value per target date")
        EKP.Observation(
            Dict(
                "samples" => [EKP.get_obs(obs), sample],
                "covariances" => [covar, LinearAlgebra.Diagonal(sd .^ 2)],
                "names" => [
                    only(EKP.get_names(obs)),
                    join(ClimaAnalysis.short_name.(metadata), ";"),
                ],
                "metadata" => vcat(EKP.get_metadata(obs), metadata),
            ),
        )
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    # Prevent MPI from being used which is not needed for generating
    # observations
    ENV["CLIMACOMMS_CONTEXT"] = "SINGLETON"

    config = CALIBRATE_CONFIG
    (; observables, sample_date_ranges, covariance_date_ranges, output_dir) = config
    isempty(observables) && error("The config has no observables")
    length(unique(covariance_date_ranges)) >= 3 ||
        error("The config needs at least 3 distinct covariance_date_ranges")

    data_loader = CalibrationTools.default_data_loader()
    vars = map(observables) do o
        source_data_loader = CalibrationTools.find_source_loader(data_loader, o.name)
        @info "Retrieving $(o.name) from $(typeof(source_data_loader))"
        return get(source_data_loader, o.name)
    end

    # For now, we apply the preprocessing to all the variables if possible
    # If the preprocessing does not apply, then it is a no-op.
    # In the future, if we want to do specific preprocessing, this needs to
    # change
    vars = select_pressure_levels.(vars, Ref(config.pressure_levels))
    lonlat_regridder = get_lonlat_regridder(config.config_file)
    vars = lonlat_regridder.(vars)

    # Record where each observation has data over the dates in use, while longitude still
    # exists. A globally covered product gets no mask. `observation_map.jl` loads the same
    # masks, and `loss_values` applies them on both sides.
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

    vars = loss_values(vars, coverage_masks, config)
    check_extend(vars, config)
    for (var, o) in zip(vars, observables)
        haskey(o, :model_error_fraction) &&
            !ClimaAnalysis.has_latitude(var) &&
            error(
                "$(o.name) has no latitude after its reduction, so give it noise_sd instead of model_error_fraction",
            )
    end

    # Give every sample the same NaN mask, or `build_samples_by_times` rejects a product
    # whose coverage varies by year.
    vars = ClimaAnalysis.propagate_nans.(vars; dims = ("time",))

    # `observation_map.jl` reads these so the simulation is sampled at the same points as
    # the observation.
    JLD2.save_object(coverage_mask_path(output_dir), coverage_masks)

    # Create the observation vector, in physical units: the noise model's model-error term
    # is a fraction of the observed mean. `beta` is the model-error variance, the square of
    # each map's model_error_fraction.
    observation_vec = make_svdplusd_observation_vector(
        vars,
        sample_date_ranges;
        covariance_date_ranges,
        beta = [
            o.model_error_fraction^2 for
            o in observables if haskey(o, :model_error_fraction)
        ],
        rank = 2,
        sigma2 = 1e-6,
        use_latitude_weights = true,
        min_cosd_lat = 0.1,
        independent_noise_sd = Dict{String, Float64}(
            o.name => o.noise_sd for o in observables if haskey(o, :noise_sd)
        ),
    )

    JLD2.save_object(observation_vector_path(output_dir), observation_vec)

    # Reconstruct the variables from the observation and show them for debugging
    for (i, obs) in enumerate(observation_vec)
        @info "Observation $i"
        @info ObservationRecipe.reconstruct_vars(obs)
    end
end
