# Opt-in end-to-end test of the CDS fetch path. Not part of the default test
# suite because it needs CDS credentials and network access, and CDS requests
# can queue for a long time. Run with:
#
#     CLIMACOUPLER_ERA5_NETWORK_TESTS=true julia --project=test test/era5_network_tests.jl

import Test: @test, @testset
import Dates
import ClimaCoupler
import ClimaCoupler.Era5InitialConditions as ERA5IC

run_network_tests =
    get(ENV, "CLIMACOUPLER_ERA5_NETWORK_TESTS", "false") == "true"

if !run_network_tests
    @info "Skipping ERA5 network tests. Set CLIMACOUPLER_ERA5_NETWORK_TESTS=true to run them."
elseif !ERA5IC.cds_credentials_available()
    error("ERA5 network tests need CDS credentials (~/.cdsapirc or CDSAPI_URL/CDSAPI_KEY)")
else
    @testset "ERA5 CDS fetch end to end" begin
        date = Dates.DateTime(2023, 6, 1)
        mktempdir() do dir
            fetched_dir =
                ERA5IC.fetch_era5_initial_conditions(date; dir, wait = 30.0)
            @test fetched_dir == dir
            @test ERA5IC.era5_files_complete(dir, date)
            ERA5IC.validate_era5_dir(dir, date)

            # The second call is a cache hit
            ERA5IC.fetch_era5_initial_conditions(date; dir)
            @test ERA5IC.era5_files_complete(dir, date)
        end
    end
end
