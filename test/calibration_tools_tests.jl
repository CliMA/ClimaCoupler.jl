import Test: @test, @testset, @test_throws, @test_logs
import Artifacts
import Dates
import ClimaCoupler
import ClimaCoupler: CalibrationTools
import ClimaAnalysis

@testset "Calibrate config" begin
    config_file =
        joinpath(pkgdir(ClimaCoupler), "experiments", "ClimaEarth", "test", "amip_test.jl")
    valid_config_args = (;
        config_file,
        short_names = ["pr"],
        minibatch_size = 2,
        n_iterations = 3,
        sample_date_ranges = [["2010-01", "2010-01"], ["2011-01", "2011-01"]],
        extend = Dates.Week(1),
        spinup = Dates.Week(1),
        output_dir = ".",
        rng_seed = 42,
    )
    config = CalibrationTools.CalibrateConfig(; valid_config_args...)

    function make_calibrate_config(; kwargs...)
        CalibrationTools.CalibrateConfig(; valid_config_args..., kwargs...)
    end

    @test_throws ErrorException make_calibrate_config(
        config_file = "this_does_not_exist.yml",
    )

    @test_throws ErrorException make_calibrate_config(short_names = String[])

    @test_throws ErrorException make_calibrate_config(sample_date_ranges = [])

    @test_throws ErrorException make_calibrate_config(extend = -Dates.Week(1))

    @test_throws ErrorException make_calibrate_config(spinup = -Dates.Week(1))

    @test_throws ErrorException make_calibrate_config(
        sample_date_ranges = [["2010-02", "2010-01"], ["2011-02", "2011-01"]],
    )

    @test_throws ErrorException make_calibrate_config(minibatch_size = -1)

    @test_throws ErrorException make_calibrate_config(n_iterations = -1)

    @test_throws ErrorException make_calibrate_config(
        minibatch_size = 3,
        sample_date_ranges = [["2010-01", "2010-01"] for _ in 1:2],
    )

    @test_logs (:warn, r"Number of samples is not divisible by the minibatch size") make_calibrate_config(
        minibatch_size = 3,
        sample_date_ranges = [["2010-01", "2010-01"] for _ in 1:7],
    )
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
    # dims: (altitude=5, lon=2)
    lons = [0.0, 1.0]
    z_i = [0.0, 1000.0, 2000.0, 3000.0, 4000.0]
    # data: cloud fraction = 0.5 everywhere
    data = fill(0.5f0, length(lons), length(z_i))
    dims = Dict("z" => z_i, "longitude" => lons)
    dim_attribs = Dict("z" => Dict{String, Any}(), "longitude" => Dict{String, Any}())
    attribs = Dict{String, Any}("short_name" => "cl", "units" => "unitless")
    obs_var = ClimaAnalysis.OutputVar(attribs, dims, dim_attribs, data)

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
