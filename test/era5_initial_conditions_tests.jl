import Test: @test, @testset, @test_throws
import Dates
import NCDatasets
import ClimaCoupler
import ClimaCoupler.Era5InitialConditions as ERA5IC

const TEST_DATE = Dates.DateTime(2024, 3, 15)
const NLON = 8
const NLAT = 5
const NLEVELS = 3

# Longitudes as CDS delivers them for area [90, -180, -90, 180]
const TEST_LON = collect(range(-180.0, 135.0, length = NLON))
# Latitudes decreasing, as in ERA5 files
const TEST_LAT = collect(range(90.0, -90.0, length = NLAT))

# A land mask for the fixtures: the first three longitudes are "land"
is_land(i, j) = i <= 3

"""
Write a fake CDS pressure-levels download.
"""
function write_fake_pressure_file(path)
    NCDatasets.NCDataset(path, "c") do ds
        NCDatasets.defDim(ds, "longitude", NLON)
        NCDatasets.defDim(ds, "latitude", NLAT)
        NCDatasets.defDim(ds, "pressure_level", NLEVELS)
        NCDatasets.defDim(ds, "valid_time", 1)
        NCDatasets.defVar(ds, "longitude", TEST_LON, ("longitude",))
        NCDatasets.defVar(ds, "latitude", TEST_LAT, ("latitude",))
        NCDatasets.defVar(
            ds,
            "pressure_level",
            [500.0, 850.0, 1000.0],
            ("pressure_level",),
        )
        NCDatasets.defVar(ds, "valid_time", [TEST_DATE], ("valid_time",))
        dims = ("longitude", "latitude", "pressure_level", "valid_time")
        for (name, value) in (
            ("z", 5.0e4),
            ("t", 280.0),
            ("q", 0.005),
            ("u", 10.0),
            ("v", -5.0),
            ("w", 0.1),
            ("clwc", 1.0e-5),
            ("ciwc", 1.0e-6),
            ("crwc", 1.0e-6),
            ("cswc", 1.0e-6),
        )
            data = fill(value, NLON, NLAT, NLEVELS, 1)
            NCDatasets.defVar(ds, name, data, dims)
        end
    end
    return path
end

"""
Write a fake CDS single-levels surface download. Masked fields (`sst`,
`siconc`, `istl*`) are missing over land.
"""
function write_fake_surface_file(path)
    NCDatasets.NCDataset(path, "c") do ds
        NCDatasets.defDim(ds, "longitude", NLON)
        NCDatasets.defDim(ds, "latitude", NLAT)
        NCDatasets.defDim(ds, "valid_time", 1)
        NCDatasets.defVar(ds, "longitude", TEST_LON, ("longitude",))
        NCDatasets.defVar(ds, "latitude", TEST_LAT, ("latitude",))
        NCDatasets.defVar(ds, "valid_time", [TEST_DATE], ("valid_time",))
        dims = ("longitude", "latitude", "valid_time")
        for (name, value) in (("skt", 285.0), ("sp", 1.0e5), ("z", 100.0))
            NCDatasets.defVar(ds, name, fill(value, NLON, NLAT, 1), dims)
        end
        masked = (
            ("sst", 290.0),
            ("siconc", 0.5),
            ("istl1", 271.0),
            ("istl2", 271.5),
            ("istl3", 272.0),
            ("istl4", 272.5),
        )
        for (name, value) in masked
            data = Array{Union{Missing, Float64}}(undef, NLON, NLAT, 1)
            for j in 1:NLAT, i in 1:NLON
                data[i, j, 1] = is_land(i, j) ? missing : value
            end
            NCDatasets.defVar(ds, name, data, dims; fillvalue = -9.0e33)
        end
    end
    return path
end

"""
Write a fake CDS single-levels land download. Soil fields are missing over
ocean.
"""
function write_fake_land_file(path)
    NCDatasets.NCDataset(path, "c") do ds
        NCDatasets.defDim(ds, "longitude", NLON)
        NCDatasets.defDim(ds, "latitude", NLAT)
        NCDatasets.defDim(ds, "valid_time", 1)
        NCDatasets.defVar(ds, "longitude", TEST_LON, ("longitude",))
        NCDatasets.defVar(ds, "latitude", TEST_LAT, ("latitude",))
        NCDatasets.defVar(ds, "valid_time", [TEST_DATE], ("valid_time",))
        dims = ("longitude", "latitude", "valid_time")
        land_only = Dict(
            "swvl1" => 0.3,
            "swvl2" => 0.32,
            "swvl3" => 0.34,
            "swvl4" => 0.36,
            "stl1" => 272.0,
            "stl2" => 281.0,
            "stl3" => 282.0,
            "stl4" => 283.0,
            "sd" => 0.1,
            "src" => 0.001,
            "tsn" => 265.0,
            "lai_lv" => 1.5,
            "lai_hv" => 2.0,
        )
        for (name, value) in land_only
            data = Array{Union{Missing, Float64}}(undef, NLON, NLAT, 1)
            for j in 1:NLAT, i in 1:NLON
                data[i, j, 1] = is_land(i, j) ? value : missing
            end
            NCDatasets.defVar(ds, name, data, dims; fillvalue = -9.0e33)
        end
        for (name, value) in
            (("skt", 285.0), ("fal", 0.2), ("fsr", 0.05), ("flsr", -4.6))
            NCDatasets.defVar(ds, name, fill(value, NLON, NLAT, 1), dims)
        end
    end
    return path
end

"""
A fake `CDSAPI.retrieve` that writes fixture files instead of calling the
network. Records the requests it receives.
"""
function make_fake_retrieve(recorded)
    return function (dataset, request, path; wait = 1.0)
        push!(recorded, (dataset, request))
        if haskey(request, "pressure_level")
            write_fake_pressure_file(path)
        elseif "sea_surface_temperature" in request["variable"]
            write_fake_surface_file(path)
        else
            write_fake_land_file(path)
        end
        return path
    end
end

@testset "filenames and cache bookkeeping" begin
    @test ERA5IC.raw_filename(TEST_DATE) == "era5_raw_20240315_0000.nc"
    @test ERA5IC.sst_filename(TEST_DATE) == "sst_processed_20240315_0000.nc"
    @test ERA5IC.sic_filename(TEST_DATE) == "sic_processed_20240315_0000.nc"
    @test ERA5IC.land_filename(TEST_DATE) ==
          "era5_land_processed_20240315_0000.nc"
    @test ERA5IC.bucket_filename(TEST_DATE) ==
          "era5_bucket_processed_20240315_0000.nc"
    @test ERA5IC.albedo_filename(TEST_DATE) ==
          "albedo_processed_20240315_0000.nc"
    @test length(ERA5IC.output_filenames(TEST_DATE)) == 6
end

@testset "field helpers" begin
    field = Union{Missing, Float64}[1.0 missing; missing 4.0]
    filled = ERA5IC.nearest_neighbor_fill(field)
    @test !any(isnan, filled)
    @test filled[1, 1] == 1.0
    @test filled[2, 2] == 4.0

    @test ERA5IC.zero_fill(field) == [1.0 0.0; 0.0 4.0]

    lon360, perm = ERA5IC.roll_longitudes([-180.0, -90.0, 0.0, 90.0])
    @test lon360 == [0.0, 90.0, 180.0, 270.0]
    @test perm == [3, 4, 1, 2]

    points = ERA5IC.monthly_time_points(Dates.DateTime(2024, 3, 15))
    @test length(points) == 4
    @test points[2] ==
          Dates.value(Dates.Second(Dates.DateTime(2024, 3, 15) - Dates.DateTime(1970, 1, 1)))
    @test issorted(points)
end

@testset "requests" begin
    request = ERA5IC.pressure_levels_request(TEST_DATE)
    @test request["year"] == "2024"
    @test request["month"] == "03"
    @test request["day"] == "15"
    @test request["data_format"] == "netcdf"
    @test length(request["pressure_level"]) == 37
    @test "temperature" in request["variable"]
    @test "vertical_velocity" in request["variable"]

    surface = ERA5IC.surface_request(TEST_DATE)
    @test "sea_surface_temperature" in surface["variable"]
    @test "sea_ice_cover" in surface["variable"]

    land = ERA5IC.land_request(TEST_DATE)
    @test "volumetric_soil_water_layer_1" in land["variable"]
    @test "forecast_albedo" in land["variable"]
end

@testset "end-to-end fetch with fake retrieval" begin
    mktempdir() do dir
        recorded = []
        fetched_dir = ERA5IC.fetch_era5_initial_conditions(
            TEST_DATE;
            dir,
            retrieve_fn = make_fake_retrieve(recorded),
        )
        @test fetched_dir == dir
        @test length(recorded) == 3
        @test ERA5IC.era5_files_complete(dir, TEST_DATE)
        for name in ERA5IC.output_filenames(TEST_DATE)
            @test isfile(joinpath(dir, name))
        end
        # No leftover temporary directories
        @test !any(startswith("tmp_era5_"), readdir(dir))

        # Validation passes on the finished cache
        ERA5IC.validate_era5_dir(dir, TEST_DATE)

        # A second fetch is a cache hit: no new requests
        ERA5IC.fetch_era5_initial_conditions(
            TEST_DATE;
            dir,
            retrieve_fn = make_fake_retrieve(recorded),
        )
        @test length(recorded) == 3

        # force = true refetches
        ERA5IC.fetch_era5_initial_conditions(
            TEST_DATE;
            dir,
            force = true,
            retrieve_fn = make_fake_retrieve(recorded),
        )
        @test length(recorded) == 6

        # Raw file: dims and merged surface variables
        NCDatasets.NCDataset(joinpath(dir, ERA5IC.raw_filename(TEST_DATE))) do ds
            for dim in ("longitude", "latitude", "pressure_level", "valid_time")
                @test haskey(ds.dim, dim)
            end
            for name in
                ("u", "v", "w", "t", "q", "z", "skt", "sp", "surface_geopotential")
                @test haskey(ds, name)
            end
            @test size(Array(ds["t"])) == (NLON, NLAT, NLEVELS, 1)
            @test all(Array(ds["surface_geopotential"]) .== 100.0f0)
        end

        # SST: Celsius, land filled, rolled longitudes, 4 time points
        NCDatasets.NCDataset(joinpath(dir, ERA5IC.sst_filename(TEST_DATE))) do ds
            sst = Array(ds["SST"])
            @test size(sst) == (NLON, NLAT, 4)
            @test all(sst .≈ 290.0 - 273.15)
            lon = Array(ds["lon"])
            @test issorted(lon)
            @test all(lon .>= 0)
            @test length(Array(ds["time"])) == 4
        end

        # SIC: percent, missing treated as 0 over land
        NCDatasets.NCDataset(joinpath(dir, ERA5IC.sic_filename(TEST_DATE))) do ds
            sic = Array(ds["SEAICE"])
            @test maximum(sic) ≈ 50.0
            @test minimum(sic) == 0.0
            @test haskey(ds, "ISTL1")
        end

        # Land: sorted latitude, ocean zeros, negative sorted z, ice partition
        NCDatasets.NCDataset(joinpath(dir, ERA5IC.land_filename(TEST_DATE))) do ds
            @test issorted(Array(ds["lat"]))
            z = Array(ds["z"])
            @test issorted(z)
            @test all(z .< 0)
            swvl = Array(ds["swvl"])
            stl = Array(ds["stl"])
            si = Array(ds["si"])
            @test size(swvl) == (NLON, NLAT, 4)
            # Ocean points are exactly 0
            @test all(swvl[4:end, :, :] .== 0)
            @test all(stl[4:end, :, :] .== 0)
            # Layer 1 (stl1 = 272 K < 273.15 K) is frozen: si == swvl there.
            # After the z reversal, layer 1 is the last z index.
            @test all(si[1:3, :, 4] .== swvl[1:3, :, 4])
            @test all(si[1:3, :, 1:3] .== 0)
            @test all(Array(ds["swe"])[1:3, :] .≈ 0.1)
            @test haskey(ds, "sie")
            @test haskey(ds, "lai")
        end

        # Bucket: W column sum with 0.5 m cutoff, filled T
        NCDatasets.NCDataset(joinpath(dir, ERA5IC.bucket_filename(TEST_DATE))) do ds
            W = Array(ds["W"])
            expected_W = 0.3 * 0.07 + 0.32 * 0.21 + 0.34 * 0.22
            @test all(W[1:3, :] .≈ expected_W)
            @test all(W[4:end, :] .== 0)
            T = Array(ds["T"])
            # Nearest-neighbor fill gives ocean points land values
            @test minimum(T) > 100
            @test all(Array(ds["Ws"])[1:3, :] .≈ 0.001)
            @test all(Array(ds["S"])[1:3, :] .≈ 0.1)
        end

        # Albedo: fal renamed, in [0, 1]
        NCDatasets.NCDataset(joinpath(dir, ERA5IC.albedo_filename(TEST_DATE))) do ds
            albedo = Array(ds["sw_alb_clr"])
            @test size(albedo, 3) == 4
            @test all(albedo .≈ 0.2)
            @test haskey(ds, "fsr")
            @test haskey(ds, "flsr")
        end
    end
end

@testset "validation failure modes" begin
    mktempdir() do dir
        recorded = []
        ERA5IC.fetch_era5_initial_conditions(
            TEST_DATE;
            dir,
            retrieve_fn = make_fake_retrieve(recorded),
        )

        # Missing file
        missing_dir = joinpath(dir, "incomplete")
        mkpath(missing_dir)
        @test_throws Exception ERA5IC.validate_era5_dir(missing_dir, TEST_DATE)

        # Corrupt a file: SEAICE outside [0, 100]
        sic_path = joinpath(dir, ERA5IC.sic_filename(TEST_DATE))
        NCDatasets.NCDataset(sic_path, "a") do ds
            ds["SEAICE"][1, 1, 1] = 150.0
        end
        @test_throws Exception ERA5IC.validate_era5_dir(dir, TEST_DATE)
    end
end

@testset "failed fetch leaves no partial cache" begin
    mktempdir() do dir
        throwing_retrieve =
            (dataset, request, path; wait = 1.0) -> error("network down")
        @test_throws Exception ERA5IC.fetch_era5_initial_conditions(
            TEST_DATE;
            dir,
            retrieve_fn = throwing_retrieve,
        )
        @test !ERA5IC.era5_files_complete(dir, TEST_DATE)
        for name in ERA5IC.output_filenames(TEST_DATE)
            @test !isfile(joinpath(dir, name))
        end
    end
end

@testset "credentials detection" begin
    withenv("CDSAPI_URL" => "https://example.invalid", "CDSAPI_KEY" => "x") do
        @test ERA5IC.cds_credentials_available()
    end
end
