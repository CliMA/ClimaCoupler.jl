import Test: @test, @testset, @test_throws, @test_logs
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

@testset "ERA5 data loader" begin
    data_loader = CalibrationTools.ERA5DataLoader()

    varnames = last.(CalibrationTools.ERA5_TO_CLIMA_NAMES)

    @test CalibrationTools.available_vars(data_loader) == Set(varnames)

    irradiance_varnames = Set(["hfls", "hfss", "rsus", "rlus"])
    @test issubset(irradiance_varnames, varnames)

    for varname in varnames
        var = get(data_loader, varname)
        @test ClimaAnalysis.short_name(var) == varname
        lons = ClimaAnalysis.longitudes(var)
        @test all(lon -> -180.0 <= lon <= 180.0, lons)
        var_dates = ClimaAnalysis.dates(var)
        @test all(Dates.day.(var_dates) .== 1)
        if ClimaAnalysis.short_name(var) in irradiance_varnames
            @test ClimaAnalysis.units(var) == "W m^-2"
        end
    end

    @test_throws ErrorException get(data_loader, "idk")
end

@testset "CERES data loader" begin
    data_loader = CalibrationTools.CERESDataLoader()

    direct_varnames = last.(CalibrationTools.CERES_TO_CLIMA_NAMES)
    all_varnames = Set(vcat(direct_varnames, collect(CalibrationTools.CERES_DERIVED_VARS)))

    @test CalibrationTools.available_vars(data_loader) == all_varnames

    for varname in direct_varnames
        var = get(data_loader, varname)
        @test ClimaAnalysis.short_name(var) == varname
        lons = ClimaAnalysis.longitudes(var)
        @test all(lon -> -180.0 <= lon <= 180.0, lons)
        var_dates = ClimaAnalysis.dates(var)
        @test all(Dates.day.(var_dates) .== 1)
        @test ClimaAnalysis.units(var) == "W m^-2"
    end

    swcre = get(data_loader, "swcre")
    @test ClimaAnalysis.short_name(swcre) == "swcre"
    @test ClimaAnalysis.units(swcre) == "W m^-2"
    rsutcs = get(data_loader, "rsutcs")
    rsut = get(data_loader, "rsut")
    @test swcre.data ≈ rsutcs.data .- rsut.data

    lwcre = get(data_loader, "lwcre")
    @test ClimaAnalysis.short_name(lwcre) == "lwcre"
    @test ClimaAnalysis.units(lwcre) == "W m^-2"
    rlutcs = get(data_loader, "rlutcs")
    rlut = get(data_loader, "rlut")
    @test lwcre.data ≈ rlutcs.data .- rlut.data

    @test_throws ErrorException get(data_loader, "idk")
end
