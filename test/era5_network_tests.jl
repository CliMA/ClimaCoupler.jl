# Tests the ERA5 initial condition resolution chain end to end, including the
# CDS download that `InitialConditions.ERA5` performs. The download itself is
# tested in that package; this checks that the coupler reaches a usable
# directory through it.
#
# Not part of the default test suite, because it needs CDS credentials and
# network access, and because CDS requests can queue for a long time. Run with:
#
#     CLIMACOUPLER_ERA5_NETWORK_TESTS=true julia --project=test test/era5_network_tests.jl

import Test: @test, @testset
import Dates
import InitialConditions.ERA5
import ClimaCoupler
import ClimaCoupler.Input

run_network_tests = get(ENV, "CLIMACOUPLER_ERA5_NETWORK_TESTS", "false") == "true"

if !run_network_tests
    @info "Skipping the ERA5 network tests. Set CLIMACOUPLER_ERA5_NETWORK_TESTS=true to run them."
elseif !ERA5.credentials_available()
    error(
        "The ERA5 network tests need CDS credentials in ~/.cdsapirc or in CDSAPI_URL and CDSAPI_KEY",
    )
else
    DATE = Dates.DateTime(2023, 6, 1)

    @testset "resolution chain reaches a usable directory" begin
        # Point the cache somewhere empty so the chain runs from scratch. It
        # resolves to the artifact or to a CDS download, depending on whether
        # the artifact covers this date. Either way the result must hold files
        # that validate for the date.
        mktempdir() do cache_dir
            withenv("INITIAL_CONDITIONS_CACHE_DIR" => cache_dir) do
                dir = Input.resolve_era5_dir(nothing, DATE)
                @test isdir(dir)
                @test ERA5.files_complete(dir, DATE)
                ERA5.validate_dir(dir, DATE)

                # The paths the component models read all exist
                paths = Input.get_era5_filepaths(
                    ClimaCoupler.Interfacer.SubseasonalMode,
                    dir,
                    DATE,
                    "",
                )
                for path in (
                    paths.sst_path,
                    paths.sic_path,
                    paths.land_ic_path,
                    paths.albedo_path,
                    paths.bucket_initial_condition,
                )
                    @test isfile(path)
                end
            end
        end
    end
end
