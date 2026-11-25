import Test: @test, @testset, @test_throws, @test_logs
import Artifacts
import Dates
import ClimaCoupler
import ClimaCoupler: CalibrationTools
import ClimaAnalysis
import ClimaAnalysis: Template

# TODO: Remove this once the name is added to ClimaAnalysis
push!(ClimaAnalysis.Var.ALTITUDE_NAMES, "height")

# The reductions observables compose, such as `annual_average` and `seasonal_anomaly`.
include(
    joinpath(
        pkgdir(ClimaCoupler),
        "experiments",
        "calibration",
        "amip",
        "preprocessing.jl",
    ),
)

@testset "Calibrate config" begin
    config_file = joinpath(pkgdir(ClimaCoupler), "experiments", "test", "amip_test.jl")
    valid_config_args = (;
        config_file,
        observables = [
            CalibrationTools.Observable(
                name = "pr",
                reduce = identity,
                model_error_fraction = 0.1,
            ),
            CalibrationTools.Observable(name = "rtmt", reduce = identity, noise_sd = 0.2),
        ],
        n_iterations = 3,
        sample_date_ranges = [["2010-01", "2010-01"], ["2011-01", "2011-01"]],
        extend = Dates.Week(1),
        spinup = Dates.Week(1),
        output_dir = ".",
        rng_seed = 42,
        # CalibrateConfig only checks that priors are given.
        priors = :test_priors,
    )
    config = CalibrationTools.CalibrateConfig(; valid_config_args...)

    # The targets are used in order, one per iteration plus one.
    @test config.sample_date_ranges ==
          [(Dates.DateTime(y, 1), Dates.DateTime(y, 1)) for y in (2010, 2011, 2010, 2011)]

    function make_calibrate_config(; kwargs...)
        CalibrationTools.CalibrateConfig(; valid_config_args..., kwargs...)
    end

    @test_throws ErrorException make_calibrate_config(
        config_file = "this_does_not_exist.yml",
    )

    @test_throws ErrorException make_calibrate_config(sample_date_ranges = [])
    @test_throws ErrorException make_calibrate_config(observables = [])
    @test_throws ErrorException make_calibrate_config(priors = nothing)

    @test_throws ErrorException make_calibrate_config(extend = -Dates.Week(1))

    @test_throws ErrorException make_calibrate_config(spinup = -Dates.Week(1))

    @test_throws ErrorException make_calibrate_config(
        sample_date_ranges = [["2010-02", "2010-01"], ["2011-02", "2011-01"]],
    )

    @test_throws ErrorException make_calibrate_config(n_iterations = -1)

    # Each observable needs a name, a reduce function and exactly one of
    # model_error_fraction and noise_sd, and names may not repeat.
    @test_throws UndefKeywordError CalibrationTools.Observable(
        name = "pr",
        model_error_fraction = 0.1,
    )
    @test_throws ErrorException CalibrationTools.Observable(name = "pr", reduce = identity)
    @test_throws ErrorException CalibrationTools.Observable(
        name = "pr",
        reduce = identity,
        model_error_fraction = 0.1,
        noise_sd = 0.2,
    )
    @test_throws ErrorException make_calibrate_config(
        observables = [
            CalibrationTools.Observable(name = "pr", reduce = identity, noise_sd = 0.2),
            CalibrationTools.Observable(name = "pr", reduce = identity, noise_sd = 0.1),
        ],
    )

    # An observable is computed from the field of its own name by default, or from `field`,
    # can set its own `effective_entries` or leave it as `nothing` to count every entry, and
    # derives `name` from `field` and `reduce` when `name` is omitted.
    pr = CalibrationTools.Observable(name = "pr", reduce = identity, noise_sd = 0.2)
    @test pr.field == "pr"
    @test isnothing(pr.effective_entries)
    @test CalibrationTools.Observable("pr", zonal_average; model_error_fraction = 0.1).name ==
          "pr_zonal_average"
    @test_throws ErrorException CalibrationTools.Observable(
        name = "pr",
        reduce = identity,
        noise_sd = 0.2,
        effective_entries = 0,
    )

    parts = mean_pattern_seasonal_observables(
        "pr";
        global_mean_sd = 0.85,
        annual_pattern_fraction = 0.1,
        seasonal_cycle_fraction = 0.3,
    )
    @test [o.name for o in parts] ==
          ["pr_global_mean", "pr_annual_pattern", "pr_seasonal_cycle"]
    @test all(o -> o.field == "pr", parts)
    @test [o.effective_entries for o in parts] == [29.0, 29.0, 29.0]

    custom_parts = mean_pattern_seasonal_observables(
        "pr";
        global_mean_sd = 0.85,
        annual_pattern_fraction = 0.1,
        seasonal_cycle_fraction = 0.3,
        effective_entries = (1, 15, 30),
    )
    @test [o.effective_entries for o in custom_parts] == [1.0, 15.0, 30.0]

    # Observables with different `effective_entries` (including `nothing` for an untempered
    # scalar) can be combined in one config, with commas or semicolons.
    config_mixed = make_calibrate_config(observables = [custom_parts, pr])
    @test length(config_mixed.observables) == 4
end

@testset "Three-part reductions" begin
    # Monthly means from December 2009 to November 2010 on a coarse grid.
    start_date = Dates.DateTime(2009, 12, 1)
    months = [start_date + Dates.Month(m) for m in 0:11]
    var =
        Template.TemplateVar() |>
        Template.add_attribs(short_name = "rsut", start_date = string(start_date)) |>
        Template.add_time_dim(
            dim = [Float64(Dates.value(Dates.Second(d - start_date))) for d in months],
        ) |>
        Template.add_lon_dim(dim = collect(0.0:60.0:300.0)) |>
        Template.add_lat_dim(dim = collect(-75.0:30.0:75.0)) |>
        Template.add_data(data = rand(12, 6, 6)) |>
        Template.initialize

    seasons = seasonal_average(var)
    annual = annual_average(var)
    global_part = annual_global_mean(var)
    pattern = annual_pattern(var)
    anomaly = seasonal_anomaly(var)

    @test ClimaAnalysis.dates(annual) == [start_date]
    @test ClimaAnalysis.dates(anomaly) == ClimaAnalysis.dates(seasons)
    # The pattern has no global mean, and a year's four anomalies sum to zero at every point.
    @test abs(only(global_mean(pattern).data)) < 1e-12
    @test maximum(abs, sum(anomaly.data; dims = 1)) < 1e-12
    # The three parts add up to the seasonal maps.
    @test anomaly.data .+ pattern.data .+ global_part.data ≈ seasons.data

    # snap_date_ranges shrinks a December-to-November range to the reduced dates inside it
    # and rejects empty or inconsistently covered ranges.
    year_range = [(start_date, start_date + Dates.Month(11))]
    @test snap_date_ranges(year_range, annual) == [(start_date, start_date)]
    @test snap_date_ranges(year_range, anomaly) ==
          [(start_date, start_date + Dates.Month(9))]
    @test_throws ErrorException snap_date_ranges(
        [(Dates.DateTime(2020, 1, 1), Dates.DateTime(2020, 12, 1))],
        annual,
    )
    @test_throws ErrorException snap_date_ranges(
        [
            (start_date, start_date + Dates.Month(11)),
            (start_date, start_date + Dates.Month(3)),
        ],
        anomaly,
    )

    # mask_and_reduce, coverage_mask, and check_extend on the three-part observables.
    config_file = joinpath(pkgdir(ClimaCoupler), "experiments", "test", "amip_test.jl")
    rsut_observables = mean_pattern_seasonal_observables(
        "rsut";
        global_mean_sd = 12.6,
        annual_pattern_fraction = 0.42,
        seasonal_cycle_fraction = 0.10,
    )
    cfg = CalibrationTools.CalibrateConfig(;
        config_file,
        observables = rsut_observables,
        n_iterations = 1,
        sample_date_ranges = year_range,
        extend = Dates.Month(1),
        spinup = Dates.Month(1),
        output_dir = ".",
        priors = :test_priors,
    )
    @test isnothing(coverage_mask(var, cfg.sample_date_ranges))
    reduced = mask_and_reduce([var], Dict{String, Any}(), cfg)
    @test [ClimaAnalysis.short_name(v) for v in reduced] ==
          ["rsut_global_mean", "rsut_annual_pattern", "rsut_seasonal_cycle"]
    @test [ClimaAnalysis.short_name(v) for v in mask_and_reduce(
              [var],
              Dict{String, Any}(),
              cfg;
              names = ["rsut_annual_pattern"],
          )] == ["rsut_annual_pattern"]
    @test isnothing(check_extend(reduced, cfg))

    short_cfg = CalibrationTools.CalibrateConfig(;
        config_file,
        observables = rsut_observables,
        n_iterations = 1,
        sample_date_ranges = [(start_date, start_date)],
        extend = Dates.Month(1),
        spinup = Dates.Month(1),
        output_dir = ".",
        priors = :test_priors,
    )
    @test_throws ErrorException check_extend(reduced, short_cfg)

    # check_extend inspects every range in sample_date_ranges, not only the first.
    second_range_short_cfg = CalibrationTools.CalibrateConfig(;
        config_file,
        observables = rsut_observables,
        n_iterations = 2,
        sample_date_ranges = [
            (start_date, start_date + Dates.Month(11)),
            (start_date, start_date),
        ],
        extend = Dates.Month(1),
        spinup = Dates.Month(1),
        output_dir = ".",
        priors = :test_priors,
    )
    @test_throws ErrorException check_extend(reduced, second_range_short_cfg)

    # A single-date range with a 12-month extend simulates 4 seasons of seasonal_anomaly
    # but selects only DJF, which check_extend rejects.
    unselected_seasons_cfg = CalibrationTools.CalibrateConfig(;
        config_file,
        observables = rsut_observables,
        n_iterations = 1,
        sample_date_ranges = [(start_date, start_date)],
        extend = Dates.Month(12),
        spinup = Dates.Month(1),
        output_dir = ".",
        priors = :test_priors,
    )
    @test_throws ErrorException check_extend(reduced, unselected_seasons_cfg)
    @test_throws ErrorException mask_and_reduce(
        ClimaAnalysis.OutputVar[],
        Dict{String, Any}(),
        cfg,
    )
end

if !isnothing(Base.find_package("ClimaCalibrate")) &&
   !isnothing(Base.find_package("EnsembleKalmanProcesses"))
    include(
        joinpath(
            pkgdir(ClimaCoupler),
            "experiments",
            "calibration",
            "amip",
            "step_size_rule.jl",
        ),
    )
    include(
        joinpath(
            pkgdir(ClimaCoupler),
            "experiments",
            "calibration",
            "amip",
            "generate_observations.jl",
        ),
    )

    @testset "EffectiveEntriesController and model_error_covariance" begin
        start_date = Dates.DateTime(2009, 12, 1)
        months = [start_date + Dates.Month(m) for m in 0:35]
        # 3 years of synthetic data whose spatial pattern crosses zero (latitudinal dipole).
        lats = collect(-75.0:30.0:75.0)
        lons = collect(0.0:60.0:300.0)
        data = [
            Float64(l / 75.0 + 0.1 * sin(t) + 0.05 * cos(lon)) for
            t in 1:length(months), lon in lons, l in lats
        ]
        var =
            Template.TemplateVar() |>
            Template.add_attribs(
                short_name = "rsut",
                units = "W m^-2",
                start_date = string(start_date),
            ) |>
            Template.add_time_dim(
                dim = [
                    Float64(Dates.value(Dates.Second(d - start_date))) for d in months
                ],
            ) |>
            Template.add_lon_dim(dim = lons, units = "degrees_east") |>
            Template.add_lat_dim(dim = lats, units = "degrees_north") |>
            Template.add_data(data = data) |>
            Template.initialize

        pattern = renamed(annual_pattern(var), "rsut_annual_pattern")
        cov_ranges = [
            (start_date + Dates.Year(y), start_date + Dates.Year(y) + Dates.Month(11)) for
            y in 0:2
        ]
        target_ranges = [first(cov_ranges)]

        # model_error_covariance with floor_mean = true raises zero-crossing entries to at
        # least the median squared mean.
        samples = SampleBuilder.build_samples_by_times(
            [pattern],
            snap_date_ranges(cov_ranges, pattern),
        )
        cov_unfloored = ObservationRecipe.covariance(
            model_error_covariance(; rank = 1, beta = 0.5, sigma = 1e-3, floor_mean = false),
            samples,
        )
        cov_floored = ObservationRecipe.covariance(
            model_error_covariance(; rank = 1, beta = 0.5, sigma = 1e-3, floor_mean = true),
            samples,
        )
        @test minimum(cov_floored.diag_cov.diag) > minimum(cov_unfloored.diag_cov.diag)

        # make_observation_vector scales covariance by n_entries / effective_entries.
        obs_untempered = only(
            make_observation_vector(
                pattern,
                target_ranges;
                covariance_date_ranges = cov_ranges,
                beta = 0.5,
                rank = 1,
            ),
        )
        obs_tempered = only(
            make_observation_vector(
                pattern,
                target_ranges;
                covariance_date_ranges = cov_ranges,
                beta = 0.5,
                rank = 1,
                effective_entries = 18.0,
            ),
        )
        n_entries = length(only(EKP.get_samples(obs_untempered)))
        factor = n_entries / 18.0
        @test only(EKP.get_covs(obs_tempered)).diag_cov.diag ≈
              factor .* only(EKP.get_covs(obs_untempered)).diag_cov.diag

        # EffectiveEntriesController uses effective_entries in place of raw entry count M.
        observables = [
            CalibrationTools.Observable(;
                name = "rsut_annual_pattern",
                field = "rsut",
                reduce = annual_pattern,
                model_error_fraction = 0.5,
                effective_entries = 18.0,
            ),
        ]
        @test effective_entries(obs_tempered, observables) == 18.0
        untempered_map = [
            CalibrationTools.Observable(;
                name = "rsut_annual_pattern",
                field = "rsut",
                reduce = annual_pattern,
                model_error_fraction = 0.5,
            ),
        ]
        @test_throws ErrorException effective_entries(obs_untempered, untempered_map)
        scalar_var = renamed(annual_global_mean(var), "rsut_global_mean")
        obs_scalar = only(
            make_observation_vector(
                scalar_var,
                target_ranges;
                noise_sd = 2.0,
            ),
        )
        untempered_scalar = [
            CalibrationTools.Observable(;
                name = "rsut_global_mean",
                field = "rsut",
                reduce = annual_global_mean,
                noise_sd = 2.0,
            ),
        ]
        @test effective_entries(obs_scalar, untempered_scalar) == 1.0
        prior = PD.constrained_gaussian("p", 1.0, 0.5, 0.0, 2.0)
        ekp = EKP.EnsembleKalmanProcess(
            obs_tempered,
            EKP.TransformUnscented(prior; impose_prior = true);
            scheduler = EffectiveEntriesController(18.0),
        )
        y = EKP.get_obs(ekp)
        g = repeat(y, 1, EKP.get_N_ens(ekp)) .+ 0.5
        EKP.update_ensemble!(ekp, g)
        @test length(EKP.get_Δt(ekp)) == 1
        @test only(EKP.get_Δt(ekp)) > 0
    end
end

@testset "Calibration utilities" begin
    # Test adding parameter toml files
    config_dict = Dict()
    parameter_toml_file1 = "parameters1.toml"
    parameter_toml_file2 = "parameters2.toml"
    CalibrationTools.add_parameter_filepath!(config_dict, parameter_toml_file1)

    @test config_dict["coupler_toml"] == [parameter_toml_file1]

    CalibrationTools.add_parameter_filepath!(config_dict, parameter_toml_file2)
    @test config_dict["coupler_toml"] == [parameter_toml_file1, parameter_toml_file2]

    # Test updating tspan
    config_dict = Dict()
    start_date = Dates.DateTime(2010, 12, 13)
    end_date = Dates.DateTime(2012)
    CalibrationTools.update_timespan!(config_dict, start_date, end_date)
    @test config_dict["start_date"] == "20101213"
    @test config_dict["t_end"] == "$((Dates.Second(end_date - start_date)).value)secs"

    config_dict["start_date"] == "20501213"
    config_dict["t_end"] == "0secs"
    CalibrationTools.update_timespan!(config_dict, start_date, end_date)
    @test config_dict["start_date"] == "20101213"
    @test config_dict["t_end"] == "$((Dates.Second(end_date - start_date)).value)secs"

    noon = Dates.DateTime(2019, 12, 31, 12)
    CalibrationTools.update_timespan!(config_dict, noon, noon + Dates.Day(1))
    @test config_dict["start_date"] == "20191231-1200"
end

"""
    check_conventions(var, varname)

Check the variable with the name `varname` that
- the short name is the same as `varname`,
- the longitudes range between -180 and 180,
- the dates are monthly.
"""
function check_conventions(var, varname)
    @test ClimaAnalysis.short_name(var) == varname
    lons = ClimaAnalysis.longitudes(var)
    @test all(lon -> -180.0 <= lon <= 180.0, lons)
    var_dates = ClimaAnalysis.dates(var)
    @test all(Dates.day.(var_dates) .== 1)
    return nothing
end

"""
    artifact_on_disk(name)

Check if the artifact exists on disk given the `name` of the artifact.
"""
function artifact_on_disk(name)
    artifacts_toml = joinpath(pkgdir(ClimaCoupler), "Artifacts.toml")
    hash = Artifacts.artifact_hash(name, artifacts_toml)
    return Artifacts.artifact_exists(hash)
end

@testset "ERA5 data loader" begin
    data_loader = CalibrationTools.ERA5DataLoader()

    varnames = last.(CalibrationTools.ERA5_TO_CLIMA_NAMES)

    @test CalibrationTools.available_vars(data_loader) == Set(varnames)

    irradiance_varnames = Set(["hfls", "hfss", "rsus", "rlus"])
    @test issubset(irradiance_varnames, varnames)

    for varname in varnames
        var = get(data_loader, varname)
        check_conventions(var, varname)
        if ClimaAnalysis.short_name(var) in irradiance_varnames
            @test ClimaAnalysis.units(var) == "W m^-2"
        end
    end

    @test_throws ErrorException get(data_loader, "idk")
end

@testset "CERES data loader" begin
    artifact_on_disk("radiation_obs") || return
    data_loader = CalibrationTools.CERESDataLoader()

    direct_varnames = last.(CalibrationTools.CERES_TO_CLIMA_NAMES)
    all_varnames = Set(vcat(direct_varnames, collect(CalibrationTools.CERES_DERIVED_VARS)))

    @test CalibrationTools.available_vars(data_loader) == all_varnames

    for varname in all_varnames
        var = get(data_loader, varname)
        check_conventions(var, varname)
        @test ClimaAnalysis.units(var) == "W m^-2"
    end

    swcre = get(data_loader, "swcre")
    rsutcs = get(data_loader, "rsutcs")
    rsut = get(data_loader, "rsut")
    @test swcre.data ≈ rsutcs.data .- rsut.data

    lwcre = get(data_loader, "lwcre")
    rlutcs = get(data_loader, "rlutcs")
    rlut = get(data_loader, "rlut")
    @test lwcre.data ≈ rlutcs.data .- rlut.data

    @test_throws ErrorException get(data_loader, "idk")
end

@testset "Calipso data loader" begin
    artifact_on_disk("calipso_cloudsat") || return
    data_loader = CalibrationTools.CalipsoDataLoader()

    varnames = last.(CalibrationTools.CALIPSO_TO_CLIMA_NAMES)
    @test CalibrationTools.available_vars(data_loader) == Set(varnames)

    for varname in varnames
        var = get(data_loader, varname)
        check_conventions(var, varname)
        @test ClimaAnalysis.units(var) == "unitless"
        @test all(0 .<= var.data .<= 1)
        @test !haskey(var.dims, "doop")
    end

    @test_throws ErrorException get(data_loader, "idk")
end

@testset "regrid_to_model_levels" begin

    # Construct a simple synthetic OutputVar with an altitude dimension
    # dims: (lon=2, altitude=5)
    lons = [0.0, 1.0]
    z_i = [0.0, 1000.0, 2000.0, 3000.0, 4000.0]
    # data: cloud fraction = 0.5 everywhere
    data = fill(0.5f0, length(lons), length(z_i))
    obs_var =
        Template.TemplateVar() |>
        Template.add_attribs(short_name = "cl", units = "unitless") |>
        Template.add_lon_dim(dim = lons) |>
        Template.add_z_dim(dim = z_i) |>
        Template.add_data(data = data) |>
        Template.initialize

    # Regrid to 3 model levels
    model_alts = [500.0, 2000.0, 3500.0]
    result = CalibrationTools.regrid_to_model_levels(obs_var, model_alts)

    @test size(result.data, 1) == length(lons)
    @test size(result.data, 2) == 3
    @test ClimaAnalysis.altitudes(result) == model_alts
    # Since data is constant 0.5, weighted average should also be 0.5
    @test all(result.data .≈ 0.5f0)

    # Test with real CALIPSO data
    if artifact_on_disk("calipso_cloudsat")
        cl_var = get(CalibrationTools.CalipsoDataLoader(), "cl")
        z_obs = ClimaAnalysis.altitudes(cl_var)
        # Pick 5 evenly-spaced model levels within the observed altitude range
        model_alts = collect(range(first(z_obs), last(z_obs), length = 5))
        result_calipso = CalibrationTools.regrid_to_model_levels(cl_var, model_alts)
        @test length(ClimaAnalysis.altitudes(result_calipso)) == 5
        @test ClimaAnalysis.altitudes(result_calipso) ≈ model_alts
        @test all(0 .<= result_calipso.data .<= 1)
    end
end

@testset "Other data loaders" begin
    data_loader_and_name_list = []
    gpcp_data_loader = CalibrationTools.GPCPDataLoader()
    gpcp_varnames = last.(CalibrationTools.GPCP_TO_CLIMA_NAMES)
    push!(data_loader_and_name_list, (gpcp_data_loader, gpcp_varnames))

    # This is not automatically downloaded
    if artifact_on_disk("era5_monthly_averages_pressure_levels_1979_2024")
        era5_pressure_level_data_loader = CalibrationTools.ERA5PressureLevelDataLoader()
        era5_pressure_level_names =
            last.(CalibrationTools.ERA5_PRESSURE_LEVEL_TO_CLIMA_NAMES)
        push!(
            data_loader_and_name_list,
            (era5_pressure_level_data_loader, era5_pressure_level_names),
        )
    end

    modis_data_loader = CalibrationTools.ModisDataLoader()
    modis_names = last.(CalibrationTools.MODIS_TO_CLIMA_NAMES)
    push!(data_loader_and_name_list, (modis_data_loader, modis_names))

    for (data_loader, varnames) in data_loader_and_name_list
        @test CalibrationTools.available_vars(data_loader) == Set(varnames)
        for varname in varnames
            var = get(data_loader, varname)
            check_conventions(var, varname)
        end
        @test_throws ErrorException get(data_loader, "idk")
    end
end

@testset "Composite data loader" begin
    era5_data_loader1 = CalibrationTools.ERA5DataLoader()
    era5_data_loader2 = CalibrationTools.ERA5DataLoader()
    gpcp_data_loader = CalibrationTools.GPCPDataLoader()

    composite_data_loader = CalibrationTools.CompositeDataLoader(era5_data_loader1)
    composite_data_loader2 =
        CalibrationTools.CompositeDataLoader(era5_data_loader1, gpcp_data_loader)
    varname_to_loader = Dict(
        v => era5_data_loader1 for v in CalibrationTools.available_vars(era5_data_loader1)
    )
    composite_data_loader3 = CalibrationTools.CompositeDataLoader(
        era5_data_loader1,
        era5_data_loader2;
        varname_to_loader,
    )

    @test CalibrationTools.available_vars(composite_data_loader) ==
          CalibrationTools.available_vars(era5_data_loader1)
    @test CalibrationTools.available_vars(composite_data_loader2) ==
          union(CalibrationTools.available_vars.([era5_data_loader1, gpcp_data_loader])...)
    @test CalibrationTools.available_vars(composite_data_loader3) ==
          CalibrationTools.available_vars(era5_data_loader1)

    rsus_from_e5dl = get(era5_data_loader1, "rsus")
    for cdl in (composite_data_loader, composite_data_loader2, composite_data_loader3)
        @test CalibrationTools.find_source_loader(cdl, "rsus") isa
              CalibrationTools.ERA5DataLoader
        rsus_from_cdl = get(cdl, "rsus")

        @test rsus_from_cdl.data == rsus_from_e5dl.data
        @test ClimaAnalysis.times(rsus_from_cdl) == ClimaAnalysis.times(rsus_from_e5dl)
        @test ClimaAnalysis.longitudes(rsus_from_cdl) ==
              ClimaAnalysis.longitudes(rsus_from_e5dl)
        @test ClimaAnalysis.latitudes(rsus_from_cdl) ==
              ClimaAnalysis.latitudes(rsus_from_e5dl)
    end


    @test_throws r"shared variable names" CalibrationTools.CompositeDataLoader(
        era5_data_loader1,
        era5_data_loader2,
    )

    @test_throws r"is not available in" CalibrationTools.CompositeDataLoader(
        era5_data_loader1,
        era5_data_loader2;
        varname_to_loader = Dict("idk" => era5_data_loader1),
    )

end
