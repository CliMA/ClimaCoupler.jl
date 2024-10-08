import Test: @test, @testset, @test_throws
import ClimaComms
@static pkgversion(ClimaComms) >= v"0.6" && ClimaComms.@import_required_backends
import ClimaCoupler
import ClimaAnalysis
import ClimaUtilities.ClimaArtifacts: @clima_artifact
import Dates

# Load file to test
include("../../experiments/ClimaEarth/user_io/leaderboard.jl")
# Data
include(joinpath(pkgdir(ClimaCoupler), "artifacts", "artifact_funcs.jl"))

@testset "Leaderboard utils" begin
    @test Leaderboard.isequispaced([1, 2, 3])
    @test !Leaderboard.isequispaced([1, 2, 4])

    input_matrix = reshape(1.0:16, (4, 4))

    @test Leaderboard.resample(input_matrix, (1.0:4, 1.0:4), ([2.8], [3.7]))[1] == 13.6

    @test_throws ErrorException Leaderboard.integration_weights(([1.0], [10.0]))
    @test_throws ErrorException Leaderboard.integration_weights(([10.0], [1.0]))
    @test_throws ErrorException Leaderboard.integration_weights(([10.0, 11.0, 13.0], [10.0]))
    @test_throws ErrorException Leaderboard.integration_weights(([10.0, 20.0], [10.0, 11.0, 13.0]))

    @test Leaderboard.integration_weights(([10.0, 20.0], [20.0, 35.0]))[1] ≈ deg2rad(10.0) * deg2rad(15.0) * cosd(20.0)

    @test Leaderboard.integrate_on_sphere(ones(361, 181), (collect(-180.0:1:180.0), collect(-90.0:1:90.0))) ≈ 1

    @test_throws ErrorException Leaderboard.mse([1], [2, 3], ([1], [2]))
    @test_throws ErrorException Leaderboard.mse([1, 2], [2, 3, 4], ([1], [2]))

    @test_throws ErrorException Leaderboard.bias([1], [2, 3], ([1], [2]))
    @test_throws ErrorException Leaderboard.bias([1, 2], [2, 3, 4], ([1], [2]))

    dates = [
        Dates.DateTime(2015, 1, 13),
        Dates.DateTime(2018, 2, 13),
        Dates.DateTime(1981, 7, 6),
        Dates.DateTime(1993, 11, 19),
        Dates.DateTime(2040, 4, 1),
        Dates.DateTime(2000, 8, 18),
    ]

    expected_dates = (
        [Dates.DateTime(2040, 4, 1)],
        [Dates.DateTime(1981, 7, 6), Dates.DateTime(2000, 8, 18)],
        [Dates.DateTime(1993, 11, 19)],
        [Dates.DateTime(2015, 1, 13), Dates.DateTime(2018, 2, 13)],
    )

    @test Leaderboard.split_by_season(dates) == expected_dates
end

@testset "Leaderboard" begin
    simdir = ClimaAnalysis.SimDir(@__DIR__)

    preprocess_fn = (data) -> data .* Float32(-1 / 86400)

    # The conversion is technically not correct for this data source, but what
    # we care about here is that preprocess_data_fn works
    sim_datasource = Leaderboard.SimDataSource(path = @__DIR__, short_name = "pr", preprocess_data_fn = preprocess_fn)

    pr = get(simdir, "pr")

    @test sim_datasource.lonlat[1] == pr.dims["lon"]
    @test sim_datasource.lonlat[2] == pr.dims["lat"]

    @test Leaderboard.data_at_date(sim_datasource, Dates.DateTime(1979, 1, 2)) == preprocess_fn(pr.data[1, :, :])

    obs_datasource = Leaderboard.ObsDataSource(;
        path = joinpath(@clima_artifact("precipitation_obs"), "gpcp.precip.mon.mean.197901-202305.nc"),
        var_name = "precip",
        preprocess_data_fn = preprocess_fn,
    )

    lat = obs_datasource.ncdataset["lat"][20]
    lon = obs_datasource.ncdataset["lon"][30]

    @test Leaderboard.find_and_resample(obs_datasource, Dates.DateTime(1979, 1, 5), ([lon], [lat]))[1] ==
          preprocess_fn(obs_datasource.ncdataset["precip"][30, 20, 1])

    # A very weak test, to make sure we can call the function.
    # We assume the key functions bias and rmse are tested elsewhere.
    computed_bias =
        Leaderboard.bias(obs_datasource, sim_datasource, [Dates.DateTime(1979, 1, 1), Dates.DateTime(1979, 1, 2)])
    @test haskey(computed_bias.attributes, "rmse")
    @test haskey(computed_bias.attributes, "bias")
end
