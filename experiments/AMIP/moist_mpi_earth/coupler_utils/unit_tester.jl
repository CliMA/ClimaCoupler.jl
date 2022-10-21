# unit_tester - contains temporaty unit tests that will be moved to the `/test` folder

import Test: @test
using Pkg
using Dates

# Pkg.activate("../")
# include("../artifacts.jl")

# import coupler utils
include("regridder.jl")
include("masker.jl")
include("../artifacts.jl")
include("../atmos/atmos_init.jl")
include("general_helper.jl")
include("bcfile_reader.jl")

# test setup
TEST_DATA_DIR = "tmp_test_data/"
isdir(TEST_DATA_DIR) ? rm(TEST_DATA_DIR; recursive = true, force = true) : nothing
REGRID_DIR = joinpath(TEST_DATA_DIR, "regrid_tmp/")

FT = Float32

date1 = DateTime(1979, 01, 01, 01, 00, 00)
date = DateTime(1979, 01, 01, 00, 00, 00)
tspan = (1, 90 * 86400) # Jan-Mar
Δt = 1 * 3600

radius = FT(6731e3)
Nq = 4

domain = Domains.SphereDomain(radius)
mesh = Meshes.EquiangularCubedSphere(domain, 4)
topology = Topologies.Topology2D(mesh)
quad = Spaces.Quadratures.GLL{Nq}()
boundary_space_t = Spaces.SpectralElementSpace2D(topology, quad)

# IO monthly
# unit test for update_midmonth_data!

land_mask_t = ClimaCore.Fields.zeros(boundary_space_t)
dummy_data = (; test_data = zeros(axes(land_mask_t)))

_info = bcfile_info_init(
    FT,
    ClimaComms.SingletonCommsContext(),
    sst_data,
    "SST",
    boundary_space_t,
    segment_idx0 = [Int(1309)],
    interpolate_daily = false,
    land_mask = land_mask_t,
)

cs_t = (;
    _info = _info,
    date = [date],
    date1 = [date1],
    SST_all = [],
    updating_dates = [],
    monthly_state_diags = (; fields = deepcopy(dummy_data), ct = [0]),
)

@info "test_update_midmonth_callback()"
function test_update_midmonth_callback()
    # step in time
    walltime = @elapsed for t in ((tspan[1] + Δt):Δt:tspan[end])
        cs_t.date[1] = current_date(t) # if not global, `date`` is not updated. Check that this still runs when distributed.

        # monthly read of boundary condition data
        @calendar_callback :(
            update_midmonth_data!(cs_t.date[1], cs_t._info),
            push!(cs_t.SST_all, deepcopy(cs_t._info.monthly_fields[1])),
            push!(cs_t.updating_dates, deepcopy(cs_t.date[1])),
        ) cs_t.date[1] next_date_in_file(cs_t._info)

    end
    @show walltime
    @test cs_t.SST_all[end] !== cs_t.SST_all[end - 1] # test if the SST field was modified 
    @test Date(cs_t.updating_dates[end]) == Date(1979, 03, 16) # check that the final file date is as expected
end
test_update_midmonth_callback()

# IO daily
@info "test_interpol()"
function test_interpol(boundary_space)
    f1 = zeros(boundary_space)
    f2 = ones(boundary_space)
    out = interpol.(f1, f2, FT(15), FT(30))
    @test parent(out)[1] ≈ FT(0.5)
end
test_interpol(boundary_space_t)

# monthly data averaging and dump
@info "test_monthly_accumulation()"
function test_monthly_accumulation()

    walltime = @elapsed for t in ((tspan[1] + Δt):Δt:tspan[end])
        cs_t.date[1] = current_date(t) # if n

        dummy_data.test_data .= 1

        # accumulate data at each timestep
        accumulate_diags(dummy_data, cs_t.monthly_state_diags)

        # monthly averages
        @calendar_callback :(
            map(x -> x ./= cs_t.monthly_state_diags.ct[1], cs_t.monthly_state_diags.fields),
            save_hdf5(cs.comms_ctx, cs_t.monthly_state_diags.fields, cs_t.date[1], TEST_DATA_DIR),
            map(x -> x .= FT(0), cs_t.monthly_state_diags.fields),
            cs_t.monthly_state_diags.ct .= FT(0),
        ) cs_t.date[1] cs_t.date1[1]
        ## step to next calendar month
        @calendar_callback :(cs_t.date1[1] += Dates.Month(1)) cs_t.date[1] cs_t.date1[1]

    end

    # read Feb data
    hdfreader = InputOutput.HDF5Reader(joinpath(TEST_DATA_DIR, "test_data.monthly_1979-02-01T00:00:00.hdf5"))
    hdf5_data = InputOutput.read_field(hdfreader, "test_data")
    @test parent(hdf5_data)[1] ≈ 1 # no hours in March 

end
test_monthly_accumulation()

@info "test_surfacearea_ones()"
# create dataset of all ones with lat/lon dimensions based on seamask.nc
function gen_data_all_ones!(path, varname)
    # Create dataset of all ones
    ds = NCDataset(path, "c")

    # define the dimensions "lon" and "lat"
    defDim(ds, "lon", 64)
    defDim(ds, "lat", 32)

    # define a global attribute
    ds.attrib["title"] = "this is an NCDataset file containing all 1s"

    # define variables
    lon = defVar(ds, "lon", Float64, ("lon",))
    lat = defVar(ds, "lat", Float64, ("lat",))
    v = defVar(ds, varname, Int64, ("lon", "lat"))

    # populate lon and lat
    lon[:] = [i for i in 0.0:(360 / 64):(360 - (360 / 64))]
    lat[:] = [i for i in (-90 + (180 / 64)):(180 / 32):(90 - (180 / 64))]

    # generate some example data and write it to v
    v[:, :] = [1 for i in 1:ds.dim["lon"], j in 1:ds.dim["lat"]]

    # write attributes
    v.attrib["comments"] = "arbitrary variable with all values 1"

    close(ds)
end
# test that remapping a dataset of all ones conserves surface area
function test_surfacearea_ones()
    varname = "surfacevar"

    mask_data = joinpath(TEST_DATA_DIR, "all_ones_.nc")
    gen_data_all_ones!(mask_data, varname)

    # init land-sea mask
    land_mask = LandSeaMask(FT, ClimaComms.SingletonCommsContext(), mask_data, varname, boundary_space_t)

    @test sum(land_mask) ≈ 4 * π * (radius^2)

end
test_surfacearea_ones()

# clean
rm(TEST_DATA_DIR; recursive = true, force = true)
