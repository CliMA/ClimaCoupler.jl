# #=
#     Unit tests for ClimaCoupler Regridder module
# =#
# import CLIMAParameters as CP
# import Thermodynamics.Parameters as TP
# import Thermodynamics as TD

# using ClimaCore: Fields
# using ClimaCoupler: Utilities, TestHelper
# using Test

# include("../experiments/AMIP/modular/components/flux_calculator.jl")

# struct TestAtmos{P, Y, D, I} <: AtmosModelSimulation
#     params::P
#     Y_init::Y
#     domain::D
#     integrator::I
# end
# function TestAtmos()
#     TestAtmos(params, )
# end

# struct TestOcean{M, Y, D, I}
#     model::M
#     Y_init::Y
#     domain::D
#     integrator::I
# end

# struct TestLand{M, Y, D, I}
#     model::M
#     Y_init::Y
#     domain::D
#     integrator::I
# end



# function create_thermo_params(FT)
#     aliases = string.(fieldnames(TD.Parameters.ThermodynamicsParameters))
#     toml_dict = CP.create_toml_dict(FT; dict_type = "alias")
#     pairs = CP.get_parameter_values!(toml_dict, aliases, "Thermodynamics")
#     TD.Parameters.ThermodynamicsParameters{FT}(; pairs...)
# end

# FT = Float64
# #for FT in (Float32, Float64) # bugs in SF Float32
#     # @testset "test calculate_and_send_turbulent_fluxes" begin
#         test_space = TestHelper.create_space(FT)
#         # Construct land fraction of 0s in top half, 1s in bottom half
#         land_fraction = Fields.ones(test_space)
#         dims = size(parent(land_fraction))
#         m = dims[1]
#         n = dims[2]
#         parent(land_fraction)[1:(m ÷ 2), :, :, :] .= FT(0)

#         # Construct ice fraction of 0s on left, 0.5s on right
#         ice_d = Fields.zeros(test_space)
#         parent(ice_d)[:, (n ÷ 2 + 1):n, :, :] .= FT(0.5)

#         coupler_cache_names = (:T_S, :albedo, :F_R_sfc, :F_R_toa, :P_liq, :P_snow, :P_net, :F_lhf, :shf, :ρτxz, :ρτyz,  :evap)
#         coupler_fields = NamedTuple{coupler_cache_names}(ntuple(i -> Fields.zeros(test_space), length(coupler_cache_names)))

#         # Fill in only the necessary parts of the simulation
#         cs = Utilities.CoupledSimulation{FT}(
#             nothing, # comms_ctx
#             nothing, # dates
#             test_space, # boundary_space
#             coupler_fields, # fields
#             nothing, # parsed_args
#             nothing, # conservation_checks
#             (Int(0), Int(1000)), # tspan
#             Int(200), # t
#             Int(200), # Δt_cpl
#             (; land = land_fraction, ice = Fields.zeros(test_space), ocean = Fields.zeros(test_space)), # surface_fractions
#             (; atmos = (; integrator = (; p = (surface_scheme = BulkScheme()),), ), ice_sim = (; integrator = (; p = (; ice_fraction = ice_d)))), # model_sims
#             (;), # mode
#             (), # diagnostics
#         )

#         thermo_params = create_thermo_params(FT)

#         calculate_and_send_turbulent_fluxes!(cs)
#         # Test that sum of fractions is 1 everywhere
#         # @test all(parent(cs.surface_fractions.ice .+ cs.surface_fractions.land .+ cs.surface_fractions.ocean) .== FT(1))
#     # end
# # end


# get_height_int_point(atmos_sim::ClimaAtmosSimulation, colidx)
# get_uv_int_point(atmos_sim::ClimaAtmosSimulation, colidx)
# get_thermo_state_point(atmos_sim::ClimaAtmosSimulation, colidx)
# get_height_sfc_point(atmos_sim::ClimaAtmosSimulation, colidx)
# get_air_density(atmos_sim::ClimaAtmosSimulation, thermo_params, thermo_state) =
# get_surface_params(atmos_sim::ClimaAtmosSimulation)

# get_temperature_point(sim::TestLand, colidx)
# get_humidity_point(sim::TestLand, colidx)
# get_z0m_point(sim::TestLand, colidx)
# get_z0b_point(sim::TestLand, colidx)
# get_beta_point(sim::TestLand, colidx)

# get_temperature_point(sim::TestOcean, colidx)
# get_humidity_point(sim::TestOcean, colidx)
# get_z0m_point(sim::TestOcean, colidx)
# get_z0b_point(sim::TestOcean, colidx)
# get_beta_point(sim::TestOcean, colidx)


# update_turbulent_fluxes_point!(sim::TestLand, fields, colidx)
# update_turbulent_fluxes_point!(sim::TestOcean, fields, colidx)

