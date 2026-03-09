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
