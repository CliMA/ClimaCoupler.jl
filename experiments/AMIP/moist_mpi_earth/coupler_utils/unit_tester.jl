# unit_tester - contains temporaty unit tests that will be moved to the `/test` folder

# unit test for update_midmonth_data! (move during interface revamp)

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

# 1. update_midmonth_data! loop
FT = Float32

date = date0 = DateTime(1979, 01, 01)
tspan = (1, 90 * 86400) # Jan-Mar
Δt = 1 * 86400

domain = Domains.SphereDomain(FT(6731e3))
mesh = Meshes.EquiangularCubedSphere(domain, 4)
topology = Topologies.Topology2D(mesh)
quad = Spaces.Quadratures.GLL{5}()
boundary_space = Spaces.SpectralElementSpace2D(topology, quad)

land_mask = ClimaCore.Fields.zeros(boundary_space)

_info = bcfile_info_init(
    FT,
    sst_data,
    "SST",
    boundary_space,
    segment_idx0 = [Int(1309)],
    interpolate_daily = false,
    land_mask = land_mask,
)

cs = (; _info = _info, date = [date], SST_all = [], updating_dates = [])

@show "Starting coupling loop"
function test_update_midmonth_callback()
    # step in time
    ct = 0
    walltime = @elapsed for t in ((tspan[1] + Δt):Δt:tspan[end])
        cs.date[1] = current_date(t) # if not global, `date`` is not updated. Check that this still runs when distributed.

        # monthly read of boundary condition data
        @calendar_callback :(
            update_midmonth_data!(cs.date[1], cs._info),
            push!(cs.SST_all, deepcopy(cs._info.monthly_fields[1])),
            push!(cs.updating_dates, deepcopy(cs.date[1])),
        ) cs.date[1] next_date_in_file(cs._info)

    end
    @show walltime
    @test cs.SST_all[end] !== cs.SST_all[end - 1] # test if the SST field was modified 
    @test Date(cs.updating_dates[end]) == Date(1979, 03, 16) # check that the final file date is as expected
end

test_update_midmonth_callback()
