#=
    Unit tests for ClimaCoupler ConservationChecker, with parsed objects mimicking those in the full coupled system
=#

using ClimaCoupler: Utilities, Regridder, TestHelper, Interfacer
# using ClimaCoupler.ConservationChecker:
#     EnergyConservationCheck, WaterConservationCheck, check_conservation!, plot_global_conservation
using ClimaCore: ClimaCore, Geometry, Meshes, Domains, Topologies, Spaces, Fields, InputOutput
import ClimaCore.InputOutput: read_field
using ClimaLSM
using ClimaComms
using Test
using NCDatasets
using Dates
using Downloads

import ClimaCoupler.Interfacer: AtmosModelSimulation, SurfaceModelSimulation, SurfaceStub, get_field, name
# import ConservationChecker:

include("../src/ConservationChecker_new.jl")

REGRID_DIR = @isdefined(REGRID_DIR) ? REGRID_DIR : joinpath("", "regrid_tmp/")
FT = Float64

# the helper functions below are temporary. They will be generalized as part of the optimization re-design (where model simulations will contain model specific info)
get_slab_energy(slab_sim, T_sfc) =
    slab_sim.integrator.p.params.ρ .* slab_sim.integrator.p.params.c .* T_sfc .* slab_sim.integrator.p.params.h


# function read_field(reader::InputOutput.HDF5Reader, name::AbstractString, get_params::Bool)
#     if get_params
#         obj = reader.file["fields/$name"]
#         pkeys = keys(obj)
#         vals = ntuple(x -> obj[pkeys[x]][:][1], length(pkeys))
#         (; zip(Symbol.(Tuple(pkeys)), vals)...)
#     else
#         read_field(reader::InputOutput.HDF5Reader, name::AbstractString)
#     end
# end

# function coupler_sim_from_file(
#     hdf5_filename;
#     Δt_cpl = 1,
#     t = 0,
#     tspan = (1, 10),
#     dates = (;),
#     coupler_fields = (;),
#     model_sims = (;),
#     mode_specifics = (;),
#     parsed_args = (;),
#     diagnostics = (),
#     conservation_checks = (;),
#     land_fraction = (),
# )

#     if isempty(model_sims)

#         reader = InputOutput.HDF5Reader(hdf5_filename)
#         atmos_u = InputOutput.read_field(reader, "atmos_u")
#         ocean_u = InputOutput.read_field(reader, "ocean_u")
#         land_fraction = InputOutput.read_field(reader, "land_mask")
#         ocean_params = InputOutput.read_field(reader, "ocean_params", true)
#         orig_coupler_fields = InputOutput.read_field(reader, "coupler_fields")
#         orig_names = (propertynames(orig_coupler_fields)...,)
#         orig_data = map(x -> getproperty(orig_coupler_fields, x), orig_names)

#         new_coupler_fields = (;
#             F_turb_energy = orig_coupler_fields.F_A,
#             F_turb_moisture = orig_coupler_fields.F_E,
#             F_radiative = orig_coupler_fields.F_R,
#             F_radiative_TOA = orig_coupler_fields.F_R_TOA,
#         )
#         new_names = (propertynames(new_coupler_fields)...,)
#         new_data = map(x -> getproperty(new_coupler_fields, x), new_names)

#         all_names = (orig_names..., new_names...)
#         all_data = (orig_data..., new_data...)

#         coupler_fields = ClimaCore.Fields.FieldVector(; zip(all_names, all_data)...)

#         close(reader)

#         as = (; integrator = (; p = (; radiation_model = nothing), u = atmos_u))
#         ls = nothing
#         os = (; integrator = (; p = (; params = ocean_params), u = ocean_u))
#         model_sims = (; atmos_sim = as, land_sim = ls, ocean_sim = os, ice_sim = nothing)
#     end

#     boundary_space = axes(land_fraction)
#     FT = eltype(land_fraction)

#     Utilities.CoupledSimulation{FT}(
#         ClimaComms.SingletonCommsContext(),
#         dates,
#         boundary_space,
#         coupler_fields,
#         parsed_args,
#         conservation_checks,
#         tspan,
#         t,
#         Δt_cpl,
#         (; land = land_fraction, ocean = FT(1) .- land_fraction, ice = land_fraction .* FT(0)),
#         model_sims,
#         mode_specifics,
#         diagnostics,
#     )
# end

struct TestAtmos{I} <: Interfacer.AtmosModelSimulation
    i::I
end
name(s::TestAtmos) = "TestAtmos"
get_field(s::TestAtmos, ::Val{:F_radiative_TOA}) = ones(s.i.space) .* 200
get_field(s::TestAtmos, ::Val{:energy}) = ones(s.i.space) .* 1e6
get_field(s::TestAtmos, ::Val{:water}) = ones(s.i.space) .* 1

struct TestOcean{I} <: Interfacer.SurfaceModelSimulation
    i::I
end
name(s::TestOcean) = "TestOcean"
get_field(s::TestOcean, ::Val{:energy}) = ones(s.i.space) .* 1e6
get_field(s::TestOcean, ::Val{:water}) = ones(s.i.space) .* 0
get_field(s::TestOcean, ::Val{:area_fraction}) = ones(s.i.space) .* 0.25

struct TestLand{I} <: Interfacer.SurfaceModelSimulation
    i::I
end
name(s::TestLand) = "TestLand"
get_field(s::TestLand, ::Val{:energy}) = ones(s.i.space) .* 0
get_field(s::TestLand, ::Val{:water}) = ones(s.i.space) .* 0
get_field(s::TestLand, ::Val{:area_fraction}) = ones(s.i.space) .* 0.25

@testset "test check_conservation for conservation" begin

    space = TestHelper.create_space(FT)
    # tmp_dir = "cons_tmp"
    # mkpath(tmp_dir)
    # local_file_0 = joinpath(tmp_dir, "coupler_conservation_0.hdf5")
    # local_file_end = joinpath(tmp_dir, "coupler_conservation_end.hdf5")
    # Downloads.download("https://caltech.box.com/shared/static/2oaft4v7elhmzlyk571vy8h0cbhzff7h.hdf5", local_file_0)
    # Downloads.download("https://caltech.box.com/shared/static/zyh0u00q9ldlr03ue2m8u25hou4p2pdn.hdf5", local_file_end)

    # set up model simulations
    atmos = TestAtmos( (; space = space))
    land = TestOcean( (; space = space))
    ocean = TestLand( (; space = space))
    ice = SurfaceStub((; area_fraction = Fields.ones(space) .* 0.5 ) )
    model_sims = (; atmos_sim = atmos, land_sim = land, ocean_sim = ocean, ice_sim = ice)

    # conservation checkers
    cc =
        (; energy = EnergyConservationCheck(model_sims), water = WaterConservationCheck(model_sims))

    # coupler fields
    cf = (;
        F_radiative_TOA = Fields.ones(space),
        P_net = Fields.zeros(space),
        P_liq = Fields.zeros(space),
        P_snow = Fields.zeros(space),
        F_turb_moisture = Fields.zeros(space),
    )
    @. cf.F_radiative_TOA = 200
    @. cf.P_liq = - 100

    # init
    cs = Utilities.CoupledSimulation{FT}(
        nothing, # comms_ctx
        nothing, # dates
        space, # boundary_space
        cf, # fields
        nothing, # parsed_args
        cc, # conservation_checks
        (Int(0), Int(1000)), # tspan
        Int(200), # t
        Int(200), # Δt_cpl
        (;), # surface_masks
        model_sims, # model_sims
        (;), # mode
        (), # diagnostics
    );

    F_r = cf.F_radiative_TOA
    P = cf.P_liq
    Δt = cs.Δt_cpl
    tot_energy_an = sum(F_r .* Δt .+ 1e6 .* 1.25)
    tot_water_an = sum(.- P .* Δt .* 0.5 .+ Fields.ones(space))
    tot_energy, tot_water = check_conservation!(cs )

    @test abs((tot_energy[1] .- tot_energy_an) / tot_energy[end]) < 1e-4
    @test abs((tot_water[1] .- tot_water_an) / tot_water[end]) < 1e-4

end

@testset "test plot_global_conservation with dummy zero models" begin
end

# @testset "test plot_global_conservation with dummy zero models" begin

#     center_3d_space = TestHelper.create_space(FT, nz = 2)
#     face_3d_space = Spaces.FaceExtrudedFiniteDifferenceSpace(center_3d_space)
#     face_var = zeros(face_3d_space)
#     face_array = reshape(parent(face_var), size(parent(face_var), 1), :)

#     land_fraction = Fields.level(zeros(center_3d_space), 1)
#     coupler_fields = (; F_radiative_TOA = land_fraction)
#     ls = nothing
#     ice_u = (; T_sfc = Fields.level(zeros(center_3d_space), 1))
#     is = (; integrator = (; p = (; params = (; h = FT(0), c = FT(0), ρ = FT(0))), u = ice_u))

#     # radiative fluxes balancing
#     radiation_model = (;
#         face_lw_flux_dn = face_array .+ FT(1),
#         face_lw_flux_up = face_array,
#         face_sw_flux_dn = face_array .- FT(1),
#         face_sw_flux_up = face_array,
#     )
#     as = (;
#         integrator = (;
#             p = (; radiation_model = radiation_model),
#             u = (; c = (; ρe_tot = zeros(center_3d_space)), f = face_var),
#         )
#     )
#     model_sims = (; atmos_sim = as, land_sim = ls, ocean_sim = nothing, ice_sim = is)

#     conservation_checks = (; energy = EnergyConservationCheck([], [], [], [], [], []))
#     coupler_sim = coupler_sim_from_file(
#         nothing,
#         conservation_checks = conservation_checks,
#         model_sims = model_sims,
#         land_fraction = land_fraction,
#         coupler_fields = coupler_fields,
#     )
#     tot = check_conservation!(coupler_sim, get_slab_energy, get_slab_energy)
#     @test tot.energy[1] ≈ 0.0

#     # radiative fluxes not balancing
#     radiation_model = (;
#         face_lw_flux_dn = face_array .+ FT(2),
#         face_lw_flux_up = face_array,
#         face_sw_flux_dn = face_array .- FT(1),
#         face_sw_flux_up = face_array,
#     )
#     as = (;
#         integrator = (;
#             p = (; radiation_model = radiation_model),
#             u = (; c = (; ρe_tot = zeros(center_3d_space)), f = face_var),
#         )
#     )
#     model_sims = (; atmos_sim = as, land_sim = ls, ocean_sim = nothing, ice_sim = is)

#     conservation_checks = (; energy = EnergyConservationCheck([], [], [], [], [], []))
#     coupler_sim = coupler_sim_from_file(
#         nothing,
#         conservation_checks = conservation_checks,
#         model_sims = model_sims,
#         land_fraction = land_fraction,
#         coupler_fields = coupler_fields,
#     )
#     tot = check_conservation!(coupler_sim, get_slab_energy, get_slab_energy)
#     @test tot.energy[1] !== 0.0
# end
