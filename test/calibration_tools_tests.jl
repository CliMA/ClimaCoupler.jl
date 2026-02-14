import Test: @test, @testset
import Dates
import ClimaCoupler: CalibrationTools
import ClimaAnalysis

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
