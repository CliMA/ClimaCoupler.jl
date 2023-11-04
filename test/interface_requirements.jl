# 1. atmosphere requirements the Interfacer and FluxCalculator, ConservationChecker

#   get_field(atmos_sim::ClimaAtmosSimulation, ::Val{:field})
(:radiative_energy_flux, :liquid_precipitation, :snow_precipitation, :turbulent_energy_flux, :turbulent_moisture_flux,
 :thermo_state_int, :height_int, :height_sfc, :uv_int, :air_density, :air_temperature, :cv_m, :gas_constant_air, :F_radiative_TOA, :energy, :water)

 get_thermo_params(sim::AtmosModelSimulation) # this is required by the FluxCalculator, and should be generalized
 get_surface_params(sim::AtmosModelSimulation) # this is required by the FluxCalculator, and should be extracted using SF.jl directly

#   update_field!(atmos_sim::ClimaAtmosSimulation, ::Val{:field}, field)
(:co2_gm, :surface_temperature, :albedo, :turbulent_fluxes)

step!(sim::AtmosModelSimulation, t)
reinit!(sim::AtmosModelSimulation)

atmos_turbulent_fluxes!(atmos_sim::ClimaAtmosSimulation, csf) # if using the CombinedStateFluxes (i.e. atmos functions to calculate turbulent fluxes)
get_model_state_vector(sim::ClimaAtmosSimulation) # required by Checkpointer

# 2. surface model requirements
#   get_field(atmos_sim::ClimaAtmosSimulation, ::Val{:field})
(:surface_temperature, :surface_humidity, :roughness_momentum, :roughness_buoyancy, :beta, :albedo, :area_fraction, :air_density, :energy, :water)

#   update_field!(atmos_sim::ClimaAtmosSimulation, ::Val{:field}, field)
(:turbulent_energy_flux, :turbulent_moisture_flux, :radiative_energy_flux, :liquid_precipitation, :snow_precipitation, :air_density)

step!(sim::SurfaceModelSimulation, t)
reinit!(sim::SurfaceModelSimulation)

update_turbulent_fluxes_point!(sim::SurfaceModelSimulation, fields::NamedTuple, colidx::Fields.ColumnIndex) # if using PartitionedStateFluxes
get_model_state_vector(sim::SurfaceModelSimulation) # required by Checkpointer
