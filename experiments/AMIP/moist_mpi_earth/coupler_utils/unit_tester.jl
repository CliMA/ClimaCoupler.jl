# unit_tester - contains temporaty unit tests that will be moved to the `/test` folder

import Test: @test
using Pkg
using Dates

# Pkg.activate("../")
# include("../artifacts.jl")

# import coupler utils
# include("regridder.jl")
# include("masker.jl")
# include("calendar_timer.jl")
# include("general_helper.jl")
# include("bcfile_reader.jl")

# test setup
FT = Float32

date = date0 = DateTime(1979, 01, 01)
tspan = (1, 90 * 86400) # Jan-Mar
Δt = 1 * 3600

domain = Domains.SphereDomain(FT(6731e3))
mesh = Meshes.EquiangularCubedSphere(domain, 4)
topology = Topologies.Topology2D(mesh)
quad = Spaces.Quadratures.GLL{5}()
boundary_space_t = Spaces.SpectralElementSpace2D(topology, quad)

# IO monthly
# unit test for update_midmonth_data!

land_mask_t = ClimaCore.Fields.zeros(boundary_space_t)

_info = bcfile_info_init(
    FT,
    sst_data,
    "SST",
    boundary_space_t,
    segment_idx0 = [Int(1309)],
    interpolate_daily = false,
    land_mask = land_mask_t,
)

cs_t = (; _info = _info, date = [date], SST_all = [], updating_dates = [])

@show "Starting coupling loop"
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
function test_interpol(boundary_space)
    f1 = zeros(boundary_space)
    f2 = ones(boundary_space)
    out = interpol.(f1, f2, FT(15), FT(30))
    @test parent(out)[1] ≈ FT(0.5)
end

test_interpol(boundary_space_t)
