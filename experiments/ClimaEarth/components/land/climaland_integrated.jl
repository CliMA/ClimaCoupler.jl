import ClimaParams
using ClimaLand
import Dates
import ClimaUtilities.TimeVaryingInputs: LinearInterpolation, PeriodicCalendar
import ClimaCoupler: Checkpointer, FieldExchanger, FluxCalculator, Interfacer, Utilities
import ClimaCore as CC
import SciMLBase
import ClimaTimeSteppers as CTS

include("climaland_helpers.jl")

"""
    ClimaLandSimulation{M, P, D, I, A}

The integrated ClimaLand model simulation object.
"""
struct ClimaLandSimulation{M, D, I, A} <: Interfacer.LandModelSimulation
    model::M
    domain::D
    integrator::I
    area_fraction::A
end

"""
    ClimaLandSimulation(
        ::Type{FT};
        tspan::Tuple{Float64, Float64},
        start_date::Dates.DateTime,
        dt,
        boundary_space,
        n_vertical_elements,
        area_fraction,
    )

Creates a land domain, ClimaLand.LandModel, and integrator, and then creates a
ClimaLandSimulation object with them.
"""
function ClimaLandSimulation(
    ::Type{FT};
    tspan::Tuple{Float64, Float64},
    start_date::Dates.DateTime,
    dt,
    boundary_space,
    area_fraction,
    n_vertical_elements = 10,
    energy_check::Bool = false,
    stepper = CTS.RK4(),
    # land_temperature_anomaly::String = "amip",
) where {FT}
    earth_param_set = ClimaLand.Parameters.LandParameters(FT)

    d_soil = FT(3.5) # soil depth

    domain = make_land_domain(boundary_space, (-d_soil, FT(0.0)), n_vertical_elements)
    spatially_varying_soil_params =
        ClimaLand.default_spatially_varying_soil_parameters(domain.space.subsurface, domain.space.surface, FT)
    # Set up the soil model args
    soil_model_type = ClimaLand.Soil.EnergyHydrology{FT}
    soil_args = create_soil_args(FT, domain, spatially_varying_soil_params)
    # Set up the soil microbes model args
    soilco2_type = ClimaLand.Soil.Biogeochemistry.SoilCO2Model{FT}
    soilco2_args = create_soilco2_args(FT, domain)
    # Set up the canopy model args
    canopy_model_args, canopy_component_args = create_canopy_args(FT, domain, earth_param_set, start_date)
    canopy_component_types = (;
        autotrophic_respiration = ClimaLand.Canopy.AutotrophicRespirationModel{FT},
        radiative_transfer = ClimaLand.Canopy.TwoStreamModel{FT},
        photosynthesis = ClimaLand.Canopy.FarquharModel{FT},
        conductance = ClimaLand.Canopy.MedlynConductanceModel{FT},
        hydraulics = ClimaLand.Canopy.PlantHydraulicsModel{FT},
        energy = ClimaLand.Canopy.BigLeafEnergyModel{FT},
    )
    # setup snow model args
    snow_args = create_snow_args(FT, domain, earth_param_set, dt)
    snow_model_type = ClimaLand.Snow.SnowModel

    # Setup overall land model
    f_over = FT(3.28) # 1/m
    R_sb = FT(1.484e-4 / 1000) # m/s
    runoff_model = ClimaLand.Soil.Runoff.TOPMODELRunoff{FT}(;
        f_over = f_over,
        f_max = spatially_varying_soil_params.f_max,
        R_sb = R_sb,
    )

    Csom = ClimaLand.PrescribedSoilOrganicCarbon{FT}(TimeVaryingInput((t) -> 5))

    # TODO remove this once CL functions are extended for CoupledAtmosphere
    precip = TimeVaryingInput((t) -> -1.0e-7)
    atmos_q = TimeVaryingInput((t) -> 0.002)
    atmos_T = TimeVaryingInput((t) -> 298.0)
    atmos_p = TimeVaryingInput((t) -> 101320)
    atmos_u = TimeVaryingInput((t) -> 3.0)
    LW_IN = TimeVaryingInput((t) -> 5.67e-8 * 298.0^4)
    SW_IN = TimeVaryingInput((t) -> 500.0)
    snow_precip = TimeVaryingInput((t) -> 0.0)
    atmos_h = FT(32)
    atmos = ClimaLand.PrescribedAtmosphere(
        precip,
        snow_precip,
        atmos_T,
        atmos_u,
        atmos_q,
        atmos_p,
        start_date,
        atmos_h,
        earth_param_set;
    )
    radiation = ClimaLand.PrescribedRadiativeFluxes(
        FT,
        SW_IN,
        LW_IN,
        start_date,
        # θs = zenith_angle,
    )
    land_input = (
        atmos = atmos,#ClimaLand.CoupledAtmosphere{FT}(),
        radiation = radiation,#ClimaLand.CoupledRadiativeFluxes{FT}(),
        # atmos = ClimaLand.CoupledAtmosphere{FT, typeof(start_date)}(start_date),
        # radiation = ClimaLand.CoupledRadiativeFluxes{FT, typeof(start_date)}(start_date),
        runoff = runoff_model,
        soil_organic_carbon = Csom,
    )

    land = ClimaLand.LandModel{FT}(;
        soilco2_type = soilco2_type, # TODO this fails because soil CO2 requires PrescribedAtmosphere currently
        soilco2_args = soilco2_args,
        land_args = land_input,
        soil_model_type = soil_model_type,
        soil_args = soil_args,
        canopy_component_types = canopy_component_types,
        canopy_component_args = canopy_component_args,
        canopy_model_args = canopy_model_args,
        snow_args = snow_args,
        snow_model_type = snow_model_type,
    )

    Y, p, cds = ClimaLand.initialize(land)
    p = get_new_cache(p, Y, energy_check)
    (; θ_r, ν) = spatially_varying_soil_params
    @. Y.soil.ϑ_l = θ_r + (ν - θ_r) / 2
    Y.soil.θ_i .= FT(0.0)
    # TODO: I am not sure why we use this T instead of T from atmosphere
    T = FT(276.85)
    ρc_s = ClimaLand.Soil.volumetric_heat_capacity.(Y.soil.ϑ_l, Y.soil.θ_i, soil_args.parameters.ρc_ds, earth_param_set)
    Y.soil.ρe_int .= ClimaLand.Soil.volumetric_internal_energy.(Y.soil.θ_i, ρc_s, T, earth_param_set)
    Y.soilco2.C .= FT(0.000412)

    plant_ν = FT(1.44e-4)
    Y.canopy.hydraulics.ϑ_l.:1 .= plant_ν
    # TODO: This should be read from atmos model?
    # set temp to T for now
    Y.canopy.energy.T .= T
    # evaluate!(Y.canopy.energy.T, atmos.T, t0)

    Y.snow.S .= 0.0
    Y.snow.S_l .= 0.0
    Y.snow.U .= 0.0

    set_initial_cache! = make_set_initial_cache(land)
    exp_tendency! = make_exp_tendency(land)
    imp_tendency! = ClimaLand.make_imp_tendency(land)
    jacobian! = ClimaLand.make_jacobian(land)
    # TODO: The cache does not have all the expected values when using a CoupledAtmosphere
    set_initial_cache!(p, Y, tspan[1])

    # set up jacobian info
    jac_kwargs = (; jac_prototype = ClimaLand.FieldMatrixWithSolver(Y), Wfact = jacobian!)

    prob = SciMLBase.ODEProblem(
        CTS.ClimaODEFunction(
            T_exp! = exp_tendency!,
            T_imp! = SciMLBase.ODEFunction(imp_tendency!; jac_kwargs...),
            dss! = ClimaLand.dss!,
        ),
        Y,
        tspan,
        p,
    )



    # TODO: correct the args here
    integrator = SciMLBase.init(prob, CTS.ExplicitAlgorithm(stepper); saveat = [tspan[1], tspan[2]], dt = dt)
    # TODO: Ignoring diagnostics for now

    # updateat = Array(tspan[1]:(3600 * 3):tspan[2])
    # drivers = ClimaLand.get_drivers(land)
    # # the update func for the coupled radiation and atmoshpere drivers should do nothing
    # updatefunc = ClimaLand.make_update_drivers(drivers)

    sim = ClimaLandSimulation(land, (; domain = domain, soil_depth = d_soil), integrator, area_fraction)
    return sim
end

###############################################################################
# Helper functions for constructing the land model
###############################################################################
"""
    create_canopy_args(
        ::Type{FT},
        domain,
        earth_param_set,
        start_date,
    )

Creates the arguments for the canopy model.
"""
function create_canopy_args(::Type{FT}, domain, earth_param_set, start_date) where {FT}
    surface_space = domain.space.surface
    (; Ω, rooting_depth, is_c3, Vcmax25, g1, G_Function, α_PAR_leaf, τ_PAR_leaf, α_NIR_leaf, τ_NIR_leaf) =
        ClimaLand.clm_canopy_parameters(surface_space)

    # Energy Balance model
    ac_canopy = FT(2.5e3)
    # Plant Hydraulics and general plant parameters
    SAI = FT(0.0) # m2/m2
    f_root_to_shoot = FT(3.5)
    RAI = FT(1.0)
    K_sat_plant = FT(5e-9) # m/s # seems much too small?
    ψ63 = FT(-4 / 0.0098) # / MPa to m, Holtzman's original parameter value is -4 MPa
    Weibull_param = FT(4) # unitless, Holtzman's original c param value
    a = FT(0.05 * 0.0098) # Holtzman's original parameter for the bulk modulus of elasticity
    conductivity_model = ClimaLand.Canopy.PlantHydraulics.Weibull{FT}(K_sat_plant, ψ63, Weibull_param)
    retention_model = ClimaLand.Canopy.PlantHydraulics.LinearRetentionCurve{FT}(a)
    plant_ν = FT(1.44e-4)
    plant_S_s = FT(1e-2 * 0.0098) # m3/m3/MPa to m3/m3/m
    n_stem = 0
    n_leaf = 1
    h_stem = FT(0.0)
    h_leaf = FT(1.0)
    zmax = FT(0.0)
    h_canopy = h_stem + h_leaf
    compartment_midpoints = n_stem > 0 ? [h_stem / 2, h_stem + h_leaf / 2] : [h_leaf / 2]
    compartment_surfaces = n_stem > 0 ? [zmax, h_stem, h_canopy] : [zmax, h_leaf]

    z0_m = FT(0.13) * h_canopy
    z0_b = FT(0.1) * z0_m

    # Individual Component arguments
    # Set up autotrophic respiration
    autotrophic_respiration_args = (; parameters = ClimaLand.Canopy.AutotrophicRespirationParameters(FT))
    # Set up radiative transfer
    radiative_transfer_args = (;
        parameters = ClimaLand.Canopy.TwoStreamParameters(
            FT;
            Ω,
            α_PAR_leaf,
            τ_PAR_leaf,
            α_NIR_leaf,
            τ_NIR_leaf,
            G_Function,
        )
    )
    # Set up conductance
    conductance_args = (; parameters = ClimaLand.Canopy.MedlynConductanceParameters(FT; g1))
    # Set up photosynthesis
    photosynthesis_args = (; parameters = ClimaLand.Canopy.FarquharParameters(FT, is_c3; Vcmax25 = Vcmax25))
    # Set up plant hydraulics
    era5_lai_artifact_path = ClimaLand.Artifacts.era5_lai_forcing_data2008_folder_path()
    era5_lai_ncdata_path = joinpath(era5_lai_artifact_path, "era5_2008_1.0x1.0_lai.nc")
    LAIfunction = ClimaLand.prescribed_lai_era5(
        era5_lai_ncdata_path,
        surface_space,
        start_date;
        time_interpolation_method = LinearInterpolation(PeriodicCalendar()),
    )
    ai_parameterization = ClimaLand.Canopy.PrescribedSiteAreaIndex{FT}(LAIfunction, SAI, RAI)

    plant_hydraulics_ps = ClimaLand.Canopy.PlantHydraulics.PlantHydraulicsParameters(;
        ai_parameterization = ai_parameterization,
        ν = plant_ν,
        S_s = plant_S_s,
        rooting_depth = rooting_depth,
        conductivity_model = conductivity_model,
        retention_model = retention_model,
    )
    plant_hydraulics_args = (
        parameters = plant_hydraulics_ps,
        n_stem = n_stem,
        n_leaf = n_leaf,
        compartment_midpoints = compartment_midpoints,
        compartment_surfaces = compartment_surfaces,
    )

    energy_args = (parameters = ClimaLand.Canopy.BigLeafEnergyParameters{FT}(ac_canopy),)

    canopy_component_args = (;
        autotrophic_respiration = autotrophic_respiration_args,
        radiative_transfer = radiative_transfer_args,
        photosynthesis = photosynthesis_args,
        conductance = conductance_args,
        hydraulics = plant_hydraulics_args,
        energy = energy_args,
    )

    shared_params = ClimaLand.Canopy.SharedCanopyParameters{FT, typeof(earth_param_set)}(z0_m, z0_b, earth_param_set)

    canopy_model_args = (; parameters = shared_params, domain = ClimaLand.obtain_surface_domain(domain))
    return canopy_model_args, canopy_component_args
end

"""
    create_soilco2_args(::Type{FT}, domain)

Creates the arguments for the soil CO2 model.
"""
function create_soilco2_args(::Type{FT}, domain) where {FT}
    soilco2_ps = ClimaLand.Soil.Biogeochemistry.SoilCO2ModelParameters(FT)
    soilco2_top_bc = ClimaLand.Soil.Biogeochemistry.AtmosCO2StateBC()
    soilco2_bot_bc = ClimaLand.Soil.Biogeochemistry.SoilCO2FluxBC((p, t) -> 0.0) # no flux
    soilco2_sources = (ClimaLand.Soil.Biogeochemistry.MicrobeProduction{FT}(),)
    soilco2_boundary_conditions = (; top = soilco2_top_bc, bottom = soilco2_bot_bc)
    return (;
        boundary_conditions = soilco2_boundary_conditions,
        sources = soilco2_sources,
        domain = domain,
        parameters = soilco2_ps,
    )
end

"""
    create_soil_args(::Type{FT}, domain, spatially_varying_soil_params)

Creates the arguments for the soil model.
"""
function create_soil_args(::Type{FT}, domain, spatially_varying_soil_params) where {FT}
    (;
        ν,
        ν_ss_om,
        ν_ss_quartz,
        ν_ss_gravel,
        hydrology_cm,
        K_sat,
        S_s,
        θ_r,
        PAR_albedo_dry,
        NIR_albedo_dry,
        PAR_albedo_wet,
        NIR_albedo_wet,
        f_max,
    ) = spatially_varying_soil_params

    soil_params = ClimaLand.Soil.EnergyHydrologyParameters(
        FT;
        ν,
        ν_ss_om,
        ν_ss_quartz,
        ν_ss_gravel,
        hydrology_cm,
        K_sat,
        S_s,
        θ_r,
        PAR_albedo_dry = PAR_albedo_dry,
        NIR_albedo_dry = NIR_albedo_dry,
        PAR_albedo_wet = PAR_albedo_wet,
        NIR_albedo_wet = NIR_albedo_wet,
    )

    return (domain = domain, parameters = soil_params)
end

"""
    create_snow_args(::Type{FT}, domain, earth_param_set, dt)

Creates the arguments for the snow model.
"""
function create_snow_args(::Type{FT}, domain, earth_param_set, dt) where {FT}
    snow_parameters = ClimaLand.Snow.SnowParameters{FT}(dt; earth_param_set = earth_param_set)
    return (; parameters = snow_parameters, domain = ClimaLand.obtain_surface_domain(domain))
end

###############################################################################
### Functions required by ClimaCoupler.jl for a SurfaceModelSimulation
###############################################################################

Interfacer.step!(sim::ClimaLandSimulation, t) = Interfacer.step!(sim.integrator, t - sim.integrator.t, true)
Interfacer.reinit!(sim::ClimaLandSimulation, t) = Interfacer.reinit!(sim.integrator, t)

function Checkpointer.get_model_prog_state(sim::ClimaLandSimulation)
    error("get_model_prog_state not implemented")
end

Interfacer.get_field(sim::ClimaLandSimulation, ::Val{:area_fraction}) = sim.area_fraction
Interfacer.get_field(sim::ClimaLandSimulation, ::Val{:beta}) =
    CL.surface_evaporative_scaling(sim.model, sim.integrator.u, sim.integrator.p)

function Interfacer.get_field(sim::ClimaLandSimulation, ::Val{:energy})
    # TODO this will be implemented in ClimaLand - see https://github.com/CliMA/ClimaLand.jl/issues/1038
    CL.total_energy(sim.integrator.u, sim.integrator.p)
    # soil_energy = sim.integrator.u.soil.ρe_int
    # canopy_energy = sim.integrator.u.canopy.energy.ρe_int
    # snow_energy = sim.integrator.u.snow.U
end
function Interfacer.get_field(sim::ClimaLandSimulation, ::Val{:roughness_buoyancy})
    z_0b_canopy = sim.model.canopy.parameters.z_0b
    z_0b_snow = sim.model.snow.parameters.z_0b
    z_0b_soil = sim.model.soil.parameters.z_0b
    return NamedTuple{(:canopy, :snow, :soil)}((z_0b_canopy, z_0b_snow, z_0b_soil))
end
function Interfacer.get_field(sim::ClimaLandSimulation, ::Val{:roughness_momentum})
    z_0m_canopy = sim.model.canopy.parameters.z_0m
    z_0m_snow = sim.model.snow.parameters.z_0m
    z_0m_soil = sim.model.soil.parameters.z_0m
    return NamedTuple{(:canopy, :snow, :soil)}((z_0m_canopy, z_0m_snow, z_0m_soil))
end
Interfacer.get_field(sim::ClimaLandSimulation, ::Val{:surface_direct_albedo}) =
    CL.surface_albedo(sim.model, sim.integrator.u, sim.integrator.p)

Interfacer.get_field(sim::ClimaLandSimulation, ::Val{:surface_diffuse_albedo}) =
    CL.surface_albedo(sim.model, sim.integrator.u, sim.integrator.p)

function Interfacer.get_field(sim::ClimaLandSimulation, ::Val{:surface_temperature})
    T_eff = sim.integrator.p.T_sfc
    T_canopy = ClimaLand.surface_temperature(sim.model.canopy, sim.integrator.u, sim.integrator.p, sim.integrator.t)
    T_snow = ClimaLand.surface_temperature(sim.model.snow, sim.integrator.u, sim.integrator.p, sim.integrator.t)
    T_soil = ClimaLand.surface_temperature(sim.model.soil, sim.integrator.u, sim.integrator.p, sim.integrator.t)
    return NamedTuple{(:effective, :canopy, :snow, :soil)}((T_eff, T_canopy, T_snow, T_soil))
end

function Interfacer.get_field(sim::ClimaLandSimulation, ::Val{:water})
    # TODO this will be implemented in ClimaLand - see https://github.com/CliMA/ClimaLand.jl/issues/1038
    return CL.total_water(sim.integrator.u, sim.integrator.p)
    # return sim.integrator.u.snow.S +
    #        sim.integrator.p.soil.θ_l +
    #        sim.integrator.u.soil.θ_i +
    #        sim.integrator.u.canopy.hydraulics.ϑ_l
end
Interfacer.get_field(sim::ClimaLandSimulation, ::Val{:emissivity}) = sim.integrator.p.ϵ_sfc
function Interfacer.get_field(sim::ClimaLandSimulation, ::Val{:surface_height})
    h_canopy = CL.surface_height(sim.model.canopy, sim.integrator.u, sim.integrator.p)
    h_snow = CL.surface_height(sim.model.snow, sim.integrator.u, sim.integrator.p)
    h_soil = CL.surface_height(sim.model.soil, sim.integrator.u, sim.integrator.p)
    return NamedTuple{(:canopy, :snow, :soil)}((h_canopy, h_snow, h_soil))
end


function Interfacer.update_field!(::ClimaLandSimulation, ::Val{:air_density}, field)
    parent(sim.integrator.p.drivers.ρ_eff) .= parent(field)
    parent(sim.integrator.p.drivers.ρ_soil) .= parent(field)
    parent(sim.integrator.p.drivers.ρ_snow) .= parent(field)
    parent(sim.integrator.p.drivers.ρ_canopy) .= parent(field)
end

function Interfacer.update_field!(sim::ClimaLandSimulation, ::Val{:air_temperature}, field)
    # should the bottom level be taken in the update_sim! function instead?
    parent(sim.integrator.p.drivers.T) .= parent(CC.Spaces.level(field), 1)
end

function Interfacer.update_field!(sim::ClimaLandSimulation, ::Val{:area_fraction}, field)
    parent(sim.area_fraction) .= parent(field)
end


function Interfacer.update_field!(sim::ClimaLandSimulation, ::Val{:liquid_precipitation}, field)
    ρ_liq = (LP.ρ_cloud_liq(sim.model.soil.parameters.earth_param_set))
    parent(sim.integrator.p.drivers.P_liq) .= parent(field ./ ρ_liq)
end

function Interfacer.update_field!(sim::ClimaLandSimulation, ::Val{:radiative_energy_flux_sfc}, field)
    p = sim.integrator.p
    snow_frac, canopy_frac, soil_frac = ClimaLand.get_component_fractions(sim)
    parent(p.drivers.LW_d) .= parent(field) ./ 2
    parent(p.drivers.SW_d) .= parent(field) ./ 2
end

function Interfacer.update_field!(sim::ClimaLandSimulation, ::Val{:snow_precipitation}, field)
    ρ_liq = (LP.ρ_cloud_liq(sim.model.parameters.earth_param_set))
    parent(sim.integrator.p.drivers.P_snow) .= parent(field ./ ρ_liq)
end

function Interfacer.update_field!(::ClimaLandSimulation, ::Val{:surface_humidity}, field)
    parent(sim.integrator.p.drivers.q_eff) .= parent(field)
    parent(sim.integrator.p.drivers.q_soil) .= parent(field)
    parent(sim.integrator.p.drivers.q_snow) .= parent(field)
    parent(sim.integrator.p.drivers.q_canopy) .= parent(field)
end

function Interfacer.update_field!(sim::ClimaLandSimulation, ::Val{:surface_pressure}, field)
    parent(sim.integrator.p.drivers.P) .= parent(field)
end

function Interfacer.update_field!(sim::ClimaLandSimulation, ::Val{:thermo_state_int}, field)
    parent(sim.integrator.p.drivers.thermal_state) .= parent(field)
end

function Interfacer.update_field!(sim::ClimaLandSimulation, ::Val{:turbulent_energy_flux}, field)
    p = sim.integrator.p
    snow_frac, canopy_frac, soil_frac = ClimaLand.get_component_fractions(sim)
    parent(p.soil.shf) .= parent(field) .* soil_frac
    parent(p.snow.shf) .= parent(field) .* snow_frac
    parent(p.canopy.energy.shf) .= parent(field) .* canopy_frac
end

function Interfacer.update_field!(sim::ClimaLandSimulation, ::Val{:turbulent_moisture_flux}, field)
    p = sim.integrator.p
    snow_frac, canopy_frac, soil_frac = ClimaLand.get_component_fractions(sim)
    parent(p.soil.vapor_flux_liq) .= parent(field) .* soil_frac
    parent(p.snow.vapor_flux) .= parent(field) .* snow_frac
end

function Interfacer.update_field!(sim::ClimaLandSimulation, ::Val{:zenith_angle}, field)
    parent(sim.integrator.p.drivers.θs) .= parent(field)
end



Interfacer.name(::ClimaLandSimulation) = "ClimaLandSimulation"
