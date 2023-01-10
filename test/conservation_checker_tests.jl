#= 
    Unit tests for ClimaCoupler ConservationChecker, with parsed objects mimicking those in the full coupled system
=#

using ClimaCoupler: Utilities, Regridder, TestHelper
using ClimaCoupler.ConservationChecker:
    EnergyConservationCheck, WaterConservationCheck, check_conservation!, plot_global_conservation
using ClimaCore: ClimaCore, Geometry, Meshes, Domains, Topologies, Spaces, Fields, InputOutput
import ClimaCore.InputOutput: read_field
using ClimaComms
using Test
using NCDatasets
using Dates
using Downloads

REGRID_DIR = @isdefined(REGRID_DIR) ? REGRID_DIR : joinpath("", "regrid_tmp/")
FT = Float64

# the helper functions below are temporary. They will be generalized as part of the optimization re-design (where model simulations will contain model specific info)
get_slab_energy(slab_sim, T_sfc) =
    slab_sim.integrator.p.params.ρ .* slab_sim.integrator.p.params.c .* T_sfc .* slab_sim.integrator.p.params.h


function read_field(reader::InputOutput.HDF5Reader, name::AbstractString, get_params::Bool)
    if get_params
        obj = reader.file["fields/$name"]
        pkeys = keys(obj)
        vals = ntuple(x -> obj[pkeys[x]][:][1], length(pkeys))
        (; zip(Symbol.(Tuple(pkeys)), vals)...)
    else
        read_field(reader::InputOutput.HDF5Reader, name::AbstractString)
    end
end

function coupler_sim_from_file(
    hdf5_filename;
    Δt_cpl = 1,
    t = 0,
    tspan = (1, 10),
    dates = (;),
    coupler_fields = (;),
    model_sims = (;),
    mode_specifics = (;),
    parsed_args = (;),
    monthly_3d_diags = (;),
    monthly_2d_diags = (;),
    conservation_checks = (;),
)

    reader = InputOutput.HDF5Reader(hdf5_filename)
    atmos_u = InputOutput.read_field(reader, "atmos_u")
    ocean_u = InputOutput.read_field(reader, "ocean_u")
    land_mask = InputOutput.read_field(reader, "land_mask")
    ocean_params = InputOutput.read_field(reader, "ocean_params", true)
    coupler_fields = InputOutput.read_field(reader, "coupler_fields")
    close(reader)

    as = (; integrator = (; p = (; radiation_model = nothing), u = atmos_u))
    ls = nothing
    os = (; integrator = (; p = (; params = ocean_params), u = ocean_u))
    model_sims = (; atmos_sim = as, land_sim = ls, ocean_sim = os, ice_sim = nothing)

    boundary_space = axes(land_mask)
    FT = eltype(land_mask)

    Utilities.CoupledSimulation{FT}(
        ClimaComms.SingletonCommsContext(),
        tspan,
        dates,
        boundary_space,
        coupler_fields,
        parsed_args,
        conservation_checks,
        t,
        Δt_cpl,
        (; land = land_mask, ocean = FT(1) .- land_mask, ice = land_mask .* FT(0)),
        model_sims,
        mode_specifics,
        monthly_3d_diags,
        monthly_2d_diags,
    )
end

@testset "test check_conservation" begin

    tmp_dir = "cons_tmp"
    mkpath(tmp_dir)
    local_file_0 = joinpath(tmp_dir, "coupler_conservation_0.hdf5")
    local_file_end = joinpath(tmp_dir, "coupler_conservation_end.hdf5")
    Downloads.download("https://caltech.box.com/shared/static/2oaft4v7elhmzlyk571vy8h0cbhzff7h.hdf5", local_file_0)
    Downloads.download("https://caltech.box.com/shared/static/zyh0u00q9ldlr03ue2m8u25hou4p2pdn.hdf5", local_file_end)
    conservation_checks =
        (; energy = EnergyConservationCheck([], [], [], [], [], []), water = WaterConservationCheck([], [], [], []))

    # init
    coupler_sim = coupler_sim_from_file(local_file_0, conservation_checks = conservation_checks)
    rad_source = deepcopy(coupler_sim.fields.F_R_TOA)
    surface_water = deepcopy(coupler_sim.fields.P_net)
    check_conservation!(coupler_sim, get_slab_energy, get_slab_energy)
    # 1 day later
    coupler_sim = coupler_sim_from_file(local_file_end, conservation_checks = conservation_checks)
    check_conservation!(coupler_sim, get_slab_energy, get_slab_energy)
    parent(rad_source) .-= parent(deepcopy(coupler_sim.fields.F_R_TOA))
    parent(surface_water) .+= parent(deepcopy(coupler_sim.fields.P_net))

    cc = conservation_checks.energy
    tot_energy = @. cc.ρe_tot_atmos + cc.ρe_tot_ocean + cc.ρe_tot_land + cc.ρe_tot_seaice + cc.toa_net_source
    err = ((tot_energy[2] - tot_energy[1]) - sum(rad_source)) / (tot_energy[2])
    @test abs(err) < 1e-5

    tot_water = conservation_checks.water.ρq_tot_atmos
    err = ((tot_water[2] - tot_water[1]) - sum(surface_water)) / (tot_water[2])
    @test abs(err) < 1e-1 # TODO: this could be due to limiters. Need more testing to reduce this error.
    rm(tmp_dir, recursive = true)
end

@testset "test plot_global_conservation" begin

    tmp_dir = "cons_tmp"
    mkpath(tmp_dir)
    local_file_0 = joinpath(tmp_dir, "coupler_conservation_0.hdf5")
    Downloads.download("https://caltech.box.com/shared/static/2oaft4v7elhmzlyk571vy8h0cbhzff7h.hdf5", local_file_0)

    conservation_checks =
        (; energy = EnergyConservationCheck([], [], [], [], [], []), water = WaterConservationCheck([], [], [], []))
    coupler_sim = coupler_sim_from_file(local_file_0, conservation_checks = conservation_checks)
    check_conservation!(coupler_sim, get_slab_energy, get_slab_energy)

    # check that plot files are being generated 
    tmpname1 = joinpath(tmp_dir, tempname(".") * ".png")
    tmpname2 = joinpath(tmp_dir, tempname(".") * ".png")
    plot_global_conservation(conservation_checks.energy, coupler_sim, figname1 = tmpname1, figname2 = tmpname2)
    @test (isfile(tmpname1) & isfile(tmpname2))
    rm(tmpname1)
    rm(tmpname2)
    plot_global_conservation(conservation_checks.water, coupler_sim, figname1 = tmpname1, figname2 = tmpname2)
    @test (isfile(tmpname1) & isfile(tmpname2))
    rm(tmp_dir, recursive = true)
end
