import ClimaParams
import ClimaLand as CL
import ClimaLand.Parameters as LP
import Dates
import ClimaUtilities.TimeVaryingInputs: LinearInterpolation, PeriodicCalendar, TimeVaryingInput
import ClimaCoupler: Checkpointer, FieldExchanger, FluxCalculator, Interfacer, Utilities
import ClimaCore as CC
import SciMLBase
import ClimaTimeSteppers as CTS
import ClimaDiagnostics as CD
import ClimaUtilities.TimeManager: ITime

include("climaland_helpers.jl")

"""
    ClimaLandSimulation{M, I, A}

The integrated ClimaLand model simulation object.

It contains the following objects:
- `model::M`: The `ClimaLand.LandModel`.
- `integrator::I`: The integrator used in timestepping this model.
- `area_fraction::A`: A ClimaCore Field representing the surface area fraction of this component model.
- `output_writer::OW`: The diagnostic output writer.
"""
struct ClimaLandSimulation{M <: ClimaLand.LandModel, I <: SciMLBase.AbstractODEIntegrator, A <: CC.Fields.Field, OW} <:
       Interfacer.LandModelSimulation
    model::M
    integrator::I
    area_fraction::A
    output_writer::OW
end

"""
    ClimaLandSimulation(
        ::Type{FT},
        dt::TT,
        tspan::Tuple{TT, TT},
        start_date::Dates.DateTime,
        output_dir::String,
        boundary_space,
        area_fraction;
        saveat::Vector{TT} = [tspan[1], tspan[2]],
        land_temperature_anomaly::String = "amip",
        use_land_diagnostics::Bool = true,
        stepper = CTS.ARS111(),
    )

Creates a ClimaLandSimulation object containing a land domain,
a ClimaLand.LandModel, and an integrator.

This type of model contains a canopy model, soil model, snow model, and
soil CO2 model. Specific details about the complexity of the model
can be found in the ClimaLand.jl documentation.
"""
function ClimaLandSimulation(
    ::Type{FT},
    dt::TT,
    tspan::Tuple{TT, TT},
    start_date::Dates.DateTime,
    output_dir::String,
    boundary_space,
    area_fraction;
    saveat::Vector{TT} = [tspan[1], tspan[2]],
    surface_elevation = CC.Fields.zeros(boundary_space),
    land_temperature_anomaly::String = "amip",
    use_land_diagnostics::Bool = true,
    stepper = CTS.ARS111(),
) where {FT, TT <: Union{Float64, ITime}}
    # Set up domain
    depth = FT(50) # soil depth
    n_vertical_elements = 10
    domain = make_land_domain(boundary_space, (-depth, FT(0.0)), n_vertical_elements)
    surface_space = domain.space.surface
    subsurface_space = domain.space.subsurface

    # Set up spatially-varying parameters
    earth_param_set = CL.Parameters.LandParameters(FT)
    spatially_varying_soil_params = CL.default_spatially_varying_soil_parameters(subsurface_space, surface_space, FT)
    # Set up the soil model args
    soil_model_type = CL.Soil.EnergyHydrology{FT}
    soil_args = create_soil_args(FT, domain, spatially_varying_soil_params)
    # Set up the soil microbes model args
    soilco2_type = CL.Soil.Biogeochemistry.SoilCO2Model{FT}
    soilco2_args = create_soilco2_args(FT, domain)
    # Set up the canopy model args
    canopy_model_args, canopy_component_args = create_canopy_args(FT, domain, earth_param_set, start_date)
    canopy_component_types = (;
        autotrophic_respiration = CL.Canopy.AutotrophicRespirationModel{FT},
        radiative_transfer = CL.Canopy.TwoStreamModel{FT},
        photosynthesis = CL.Canopy.FarquharModel{FT},
        conductance = CL.Canopy.MedlynConductanceModel{FT},
        hydraulics = CL.Canopy.PlantHydraulicsModel{FT},
        energy = CL.Canopy.BigLeafEnergyModel{FT},
    )
    # setup snow model args
    snow_args = create_snow_args(FT, domain, earth_param_set, dt)
    snow_model_type = CL.Snow.SnowModel

    # Setup overall land model
    f_over = FT(3.28) # 1/m
    R_sb = FT(1.484e-4 / 1000) # m/s
    runoff_model =
        CL.Soil.Runoff.TOPMODELRunoff{FT}(; f_over = f_over, f_max = spatially_varying_soil_params.f_max, R_sb = R_sb)

    Csom = CL.PrescribedSoilOrganicCarbon{FT}(TimeVaryingInput((t) -> 5))

    land_input = (
        atmos = CL.CoupledAtmosphere{FT}(boundary_space),
        radiation = CL.CoupledRadiativeFluxes{FT}(),
        runoff = runoff_model,
        soil_organic_carbon = Csom,
    )

    model = CL.LandModel{FT}(;
        soilco2_type = soilco2_type,
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

    Y, p, coords = CL.initialize(model)

    # Set initial conditions
    (; θ_r, ν) = spatially_varying_soil_params

    # Apply temperature anomaly function to initial temperature
    T_functions = Dict("aquaplanet" => temp_anomaly_aquaplanet, "amip" => temp_anomaly_amip)
    haskey(T_functions, land_temperature_anomaly) ||
        error("land temp anomaly function $land_temperature_anomaly not supported")
    temp_anomaly = T_functions[land_temperature_anomaly]
    T_sfc0 = FT(276.85) .+ temp_anomaly.(coords.subsurface)
    lapse_rate = FT(6.5e-3)
    # Adjust initial temperature to account for orography of the surface
    # `surface_elevation` is a ClimaCore.Fields.Field(`half` level)
    orog_adjusted_T_data = CC.Fields.field_values(T_sfc0) .- lapse_rate .* CC.Fields.field_values(surface_elevation)
    orog_adjusted_T = CC.Fields.Field(orog_adjusted_T_data, subsurface_space)
    orog_adjusted_T_surface = CC.Fields.Field(CC.Fields.level(orog_adjusted_T_data, 1), surface_space)

    # Set initial conditions for the state
    @. Y.soil.ϑ_l = θ_r + (ν - θ_r) / 2
    Y.soil.θ_i .= FT(0.0)
    ρc_s = CL.Soil.volumetric_heat_capacity.(Y.soil.ϑ_l, Y.soil.θ_i, soil_args.parameters.ρc_ds, earth_param_set)
    Y.soil.ρe_int .= CL.Soil.volumetric_internal_energy.(Y.soil.θ_i, ρc_s, orog_adjusted_T, earth_param_set)
    Y.soilco2.C .= FT(0.000412) # set to atmospheric co2, mol co2 per mol air
    Y.canopy.hydraulics.ϑ_l.:1 .= canopy_component_args.hydraulics.parameters.ν
    @. Y.canopy.energy.T = orog_adjusted_T_surface

    Y.snow.S .= FT(0)
    Y.snow.S_l .= FT(0)
    Y.snow.U .= FT(0)

    set_initial_cache! = CL.make_set_initial_cache(model)
    exp_tendency! = CL.make_exp_tendency(model)
    imp_tendency! = CL.make_imp_tendency(model)
    jacobian! = CL.make_jacobian(model)
    set_initial_cache!(p, Y, tspan[1])

    # set up jacobian info
    jac_kwargs = (; jac_prototype = CL.FieldMatrixWithSolver(Y), Wfact = jacobian!)

    prob = SciMLBase.ODEProblem(
        CTS.ClimaODEFunction(T_exp! = exp_tendency!, T_imp! = SciMLBase.ODEFunction(imp_tendency!; jac_kwargs...)),
        Y,
        tspan,
        p,
    )

    # Set up diagnostics
    if use_land_diagnostics
        output_writer = CD.Writers.NetCDFWriter(subsurface_space, output_dir; start_date)
        scheduled_diagnostics = CL.default_diagnostics(
            model,
            start_date,
            output_writer = output_writer,
            output_vars = :short,
            average_period = :monthly,
        )
        diagnostic_handler = CD.DiagnosticsHandler(scheduled_diagnostics, Y, p, tspan[1]; dt = dt)
        diag_cb = CD.DiagnosticsCallback(diagnostic_handler)
    else
        output_writer = nothing
        diag_cb = nothing
    end

    # Set up time stepper and integrator
    ode_algo =
        CTS.IMEXAlgorithm(stepper, CTS.NewtonsMethod(max_iters = 3, update_j = CTS.UpdateEvery(CTS.NewNewtonIteration)))
    integrator = SciMLBase.init(
        prob,
        ode_algo;
        dt = dt,
        saveat = saveat,
        adaptive = false,
        callback = SciMLBase.CallbackSet(diag_cb),
    )

    return ClimaLandSimulation(model, integrator, area_fraction, output_writer)
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
        CL.clm_canopy_parameters(surface_space)

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
    conductivity_model = CL.Canopy.PlantHydraulics.Weibull{FT}(K_sat_plant, ψ63, Weibull_param)
    retention_model = CL.Canopy.PlantHydraulics.LinearRetentionCurve{FT}(a)
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
    autotrophic_respiration_args = (; parameters = CL.Canopy.AutotrophicRespirationParameters(FT))
    # Set up radiative transfer
    radiative_transfer_args = (;
        parameters = CL.Canopy.TwoStreamParameters(FT; Ω, α_PAR_leaf, τ_PAR_leaf, α_NIR_leaf, τ_NIR_leaf, G_Function)
    )
    # Set up conductance
    conductance_args = (; parameters = CL.Canopy.MedlynConductanceParameters(FT; g1))
    # Set up photosynthesis
    photosynthesis_args = (; parameters = CL.Canopy.FarquharParameters(FT, is_c3; Vcmax25 = Vcmax25))
    # Set up plant hydraulics
    modis_lai_artifact_path = CL.Artifacts.modis_lai_forcing_data_path()
    modis_lai_ncdata_path = joinpath(modis_lai_artifact_path, "Yuan_et_al_2008_1x1.nc")
    LAIfunction = CL.prescribed_lai_modis(
        modis_lai_ncdata_path,
        surface_space,
        start_date;
        time_interpolation_method = LinearInterpolation(PeriodicCalendar()),
    )
    ai_parameterization = CL.Canopy.PrescribedSiteAreaIndex{FT}(LAIfunction, SAI, RAI)

    plant_hydraulics_ps = CL.Canopy.PlantHydraulics.PlantHydraulicsParameters(;
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

    energy_args = (parameters = CL.Canopy.BigLeafEnergyParameters{FT}(ac_canopy),)

    canopy_component_args = (;
        autotrophic_respiration = autotrophic_respiration_args,
        radiative_transfer = radiative_transfer_args,
        photosynthesis = photosynthesis_args,
        conductance = conductance_args,
        hydraulics = plant_hydraulics_args,
        energy = energy_args,
    )

    shared_params = CL.Canopy.SharedCanopyParameters{FT, typeof(earth_param_set)}(z0_m, z0_b, earth_param_set)

    canopy_model_args = (; parameters = shared_params, domain = CL.obtain_surface_domain(domain))
    return canopy_model_args, canopy_component_args
end

"""
    create_soilco2_args(::Type{FT}, domain)

Creates the arguments for the soil CO2 model.
"""
function create_soilco2_args(::Type{FT}, domain) where {FT}
    soilco2_ps = CL.Soil.Biogeochemistry.SoilCO2ModelParameters(FT)
    soilco2_top_bc = CL.Soil.Biogeochemistry.AtmosCO2StateBC()
    soilco2_bot_bc = CL.Soil.Biogeochemistry.SoilCO2FluxBC((p, t) -> 0.0) # no flux
    soilco2_sources = (CL.Soil.Biogeochemistry.MicrobeProduction{FT}(),)
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

    soil_params = CL.Soil.EnergyHydrologyParameters(
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
    snow_parameters = CL.Snow.SnowParameters{FT}(dt; earth_param_set = earth_param_set)
    return (; parameters = snow_parameters, domain = CL.obtain_surface_domain(domain))
end

###############################################################################
### Functions required by ClimaCoupler.jl for a SurfaceModelSimulation
###############################################################################

Interfacer.get_field(sim::ClimaLandSimulation, ::Val{:area_fraction}) = sim.area_fraction
Interfacer.get_field(sim::ClimaLandSimulation, ::Val{:beta}) =
    CL.surface_evaporative_scaling(sim.model, sim.integrator.u, sim.integrator.p)
Interfacer.get_field(sim::ClimaLandSimulation, ::Val{:emissivity}) = sim.integrator.p.ϵ_sfc
Interfacer.get_field(sim::ClimaLandSimulation, ::Val{:energy}) = CL.total_energy(sim.integrator.u, sim.integrator.p)
Interfacer.get_field(sim::ClimaLandSimulation, ::Val{:surface_direct_albedo}) =
    CL.surface_albedo(sim.model, sim.integrator.u, sim.integrator.p)
Interfacer.get_field(sim::ClimaLandSimulation, ::Val{:surface_diffuse_albedo}) =
    CL.surface_albedo(sim.model, sim.integrator.u, sim.integrator.p)
Interfacer.get_field(sim::ClimaLandSimulation, ::Val{:water}) = CL.total_water(sim.integrator.u, sim.integrator.p)
Interfacer.get_field(sim::ClimaLandSimulation, ::Val{:surface_temperature}) = sim.integrator.p.T_sfc

# Update fields stored in land model sub-components
function Interfacer.update_field!(sim::ClimaLandSimulation, ::Val{:area_fraction}, field)
    parent(sim.area_fraction) .= parent(field)
end
function Interfacer.update_field!(sim::ClimaLandSimulation, ::Val{:diffuse_fraction}, field)
    parent(sim.integrator.p.drivers.frac_diff) .= parent(field)
end

# Update fields stored in land drivers
function Interfacer.update_field!(sim::ClimaLandSimulation, ::Val{:air_temperature}, field)
    parent(sim.integrator.p.drivers.T) .= parent(field)
end
function Interfacer.update_field!(sim::ClimaLandSimulation, ::Val{:air_pressure}, field)
    parent(sim.integrator.p.drivers.P) .= parent(field)
end
function Interfacer.update_field!(sim::ClimaLandSimulation, ::Val{:air_humidity}, field)
    parent(sim.integrator.p.drivers.q) .= parent(field)
end
function Interfacer.update_field!(sim::ClimaLandSimulation, ::Val{:c_co2}, field)
    sim.integrator.p.drivers.c_co2 .= field
end
function Interfacer.update_field!(sim::ClimaLandSimulation, ::Val{:liquid_precipitation}, field)
    # Arbitrarily take parameters from the soil (they are the same for all land sub-components)
    ρ_liq = (LP.ρ_cloud_liq(sim.model.soil.parameters.earth_param_set))
    parent(sim.integrator.p.drivers.P_liq) .= parent(field ./ ρ_liq)
end
function Interfacer.update_field!(sim::ClimaLandSimulation, ::Val{:snow_precipitation}, field)
    # Arbitrarily take parameters from the soil (they are the same for all land sub-components)
    ρ_liq = (LP.ρ_cloud_liq(sim.model.soil.parameters.earth_param_set))
    parent(sim.integrator.p.drivers.P_snow) .= parent(field ./ ρ_liq)
end
function Interfacer.update_field!(sim::ClimaLandSimulation, ::Val{:lw_d}, field)
    parent(sim.integrator.p.drivers.LW_d) .= parent(field)
end
function Interfacer.update_field!(sim::ClimaLandSimulation, ::Val{:sw_d}, field)
    parent(sim.integrator.p.drivers.SW_d) .= parent(field)
end
function Interfacer.update_field!(sim::ClimaLandSimulation, ::Val{:cos_zenith_angle}, field)
    parent(sim.integrator.p.drivers.cosθs) .= parent(field)
end

Interfacer.step!(sim::ClimaLandSimulation, t) = Interfacer.step!(sim.integrator, t - sim.integrator.t, true)
Interfacer.reinit!(sim::ClimaLandSimulation, t) = Interfacer.reinit!(sim.integrator, t)

function FieldExchanger.update_sim!(sim::ClimaLandSimulation, csf, area_fraction)
    # update fields for radiative transfer
    Interfacer.update_field!(sim, Val(:diffuse_fraction), csf.diffuse_fraction)
    Interfacer.update_field!(sim, Val(:sw_d), csf.SW_d)
    Interfacer.update_field!(sim, Val(:lw_d), csf.LW_d)
    Interfacer.update_field!(sim, Val(:cos_zenith_angle), csf.cos_zenith_angle)

    # update fields for canopy conductance and photosynthesis
    Interfacer.update_field!(sim, Val(:c_co2), csf.c_co2)
    Interfacer.update_field!(sim, Val(:air_temperature), csf.T_air)
    Interfacer.update_field!(sim, Val(:air_pressure), csf.P_air)
    Interfacer.update_field!(sim, Val(:air_humidity), csf.q_air)

    # precipitation
    Interfacer.update_field!(sim, Val(:liquid_precipitation), csf.P_liq)
    Interfacer.update_field!(sim, Val(:snow_precipitation), csf.P_snow)
end

function FieldExchanger.import_atmos_fields!(csf, sim::ClimaLandSimulation, atmos_sim)
    FieldExchanger.dummmy_remap!(csf.diffuse_fraction, Interfacer.get_field(atmos_sim, Val(:diffuse_fraction)))
    FieldExchanger.dummmy_remap!(csf.SW_d, Interfacer.get_field(atmos_sim, Val(:SW_d)))
    FieldExchanger.dummmy_remap!(csf.LW_d, Interfacer.get_field(atmos_sim, Val(:LW_d)))
    FieldExchanger.dummmy_remap!(csf.cos_zenith_angle, Interfacer.get_field(atmos_sim, Val(:cos_zenith)))
    FieldExchanger.dummmy_remap!(csf.P_air, Interfacer.get_field(atmos_sim, Val(:air_pressure)))
    FieldExchanger.dummmy_remap!(csf.T_air, Interfacer.get_field(atmos_sim, Val(:air_temperature)))
    FieldExchanger.dummmy_remap!(csf.q_air, Interfacer.get_field(atmos_sim, Val(:specific_humidity)))
    FieldExchanger.dummmy_remap!(csf.P_liq, Interfacer.get_field(atmos_sim, Val(:liquid_precipitation)))
    FieldExchanger.dummmy_remap!(csf.P_snow, Interfacer.get_field(atmos_sim, Val(:snow_precipitation)))
    # CO2 is a scalar for now so it doesn't need remapping
    csf.c_co2 .= Interfacer.get_field(atmos_sim, Val(:co2))
end

"""
Extend Interfacer.add_coupler_fields! to add the fields required for ClimaLandSimulation.

The fields added are:
- `:SW_d` (for radiative transfer)
- `:LW_d` (for radiative transfer)
- `:cos_zenith_angle` (for radiative transfer)
- `:diffuse_fraction` (for radiative transfer)
- `:c_co2` (for photosynthesis, biogeochemistry)
- `:P_air` (for canopy conductance)
- `:T_air` (for canopy conductance)
- `:q_air` (for canopy conductance)
- `P_liq` (for moisture fluxes)
- `P_snow` (for moisture fluxes)
"""
function Interfacer.add_coupler_fields!(coupler_field_names, ::ClimaLandSimulation)
    land_coupler_fields =
        [:SW_d, :LW_d, :cos_zenith_angle, :diffuse_fraction, :c_co2, :P_air, :T_air, :q_air, :P_liq, :P_snow]
    push!(coupler_field_names, land_coupler_fields...)
end

function Checkpointer.get_model_prog_state(sim::ClimaLandSimulation)
    error("get_model_prog_state not implemented")
end

Interfacer.name(::ClimaLandSimulation) = "ClimaLandSimulation"

## Extend functions for land-specific flux calculation
"""
    compute_surface_fluxes!(csf, sim::ClimaLandSimulation, atmos_sim, boundary_space, thermo_params, surface_scheme)

This function computes surface fluxes between the integrated land model
simulation and the atmosphere.

Update the input coupler surface fields `csf` in-place with the computed fluxes
for this model. These are then summed using area-weighting across all surface
models to get the total fluxes. Fluxes where the area fraction is zero are set to zero.

Because the integrated land model is composed of multiple sub-components, the
fluxes are computed for each sub-component and then combined to get the total for this model.
The land model cache is updated with the computed fluxes for each sub-component.

# Arguments
- `csf`: [CC.Fields.Field] containing a NamedTuple of turbulent flux fields: `F_turb_ρτxz`, `F_turb_ρτyz`, `F_turb_energy`, `F_turb_moisture`.
- `sim`: [ClimaLandSimulation] the integrated land simulation to compute fluxes for.
- `atmos_sim`: [Interfacer.AtmosModelSimulation] the atmosphere simulation to compute fluxes with.
- unused arguments: `boundary_space`, `thermo_params`, `surface_scheme`
"""
function FluxCalculator.compute_surface_fluxes!(
    csf,
    sim::ClimaLandSimulation,
    atmos_sim::Interfacer.AtmosModelSimulation,
    _...,
)
    coupled_atmos = sim.model.soil.boundary_conditions.top.atmos

    # `_int` refers to atmos state of center level 1
    z_int = Interfacer.get_field(atmos_sim, Val(:height_int))
    uₕ_int = Interfacer.get_field(atmos_sim, Val(:uv_int))
    thermo_state_int = Interfacer.get_field(atmos_sim, Val(:thermo_state_int))

    # We use `field_values` here because of the mismatch between the lowest atmos level
    #  and the boundary space, even though the horizontal grids are the same.
    @assert axes(coupled_atmos.u).grid == axes(uₕ_int).grid.full_grid.horizontal_grid
    @assert axes(coupled_atmos.thermal_state).grid == axes(thermo_state_int).grid.full_grid.horizontal_grid
    @assert axes(coupled_atmos.h).grid == axes(z_int).grid.full_grid.horizontal_grid
    Fields.field_values(coupled_atmos.u) .= Fields.field_values(uₕ_int)
    Fields.field_values(coupled_atmos.thermal_state) .= Fields.field_values(thermo_state_int)
    Fields.field_values(coupled_atmos.h) .= Fields.field_values(z_int)

    # set the same atmosphere state for all sub-components
    @assert sim.model.soil.boundary_conditions.top.atmos ===
            sim.model.canopy.boundary_conditions.atmos ===
            sim.model.snow.boundary_conditions.atmos ===
            coupled_atmos

    # get area fraction of the land model (min = 0, max = 1)
    area_fraction = Interfacer.get_field(sim, Val(:area_fraction))
    Y, p, t, model = sim.integrator.u, sim.integrator.p, sim.integrator.t, sim.model

    # compute the fluxes for each sub-component and update the land model cache
    soil_dest = p.soil.turbulent_fluxes
    ClimaLand.coupler_compute_turbulent_fluxes!(soil_dest, coupled_atmos, model.soil, Y, p, t)

    snow_dest = p.snow.turbulent_fluxes
    ClimaLand.coupler_compute_turbulent_fluxes!(snow_dest, coupled_atmos, model.snow, Y, p, t)

    canopy_dest = p.canopy.turbulent_fluxes
    ClimaLand.coupler_compute_turbulent_fluxes!(canopy_dest, coupled_atmos, model.canopy, Y, p, t)

    # Combine turbulent energy fluxes from each component of the land model
    # Use temporary variables to avoid allocating
    @. csf.temp1 =
        canopy_dest.lhf + soil_dest.lhf * (1 - p.snow.snow_cover_fraction) + p.snow.snow_cover_fraction * snow_dest.lhf
    @. csf.temp2 =
        canopy_dest.shf + soil_dest.shf * (1 - p.snow.snow_cover_fraction) + p.snow.snow_cover_fraction * snow_dest.shf

    # Zero out the fluxes where the area fraction is zero
    @. csf.temp1 = ifelse(area_fraction == 0, zero(csf.temp1), csf.temp1)
    @. csf.temp2 = ifelse(area_fraction == 0, zero(csf.temp2), csf.temp2)

    # Update the coupler field in-place
    @. csf.F_turb_energy += (csf.temp1 .+ csf.temp2) * area_fraction

    # Combine turbulent moisture fluxes from each component of the land model
    @. csf.temp1 =
        canopy_dest.transpiration +
        (soil_dest.vapor_flux_liq + soil_dest.vapor_flux_ice) * (1 - p.snow.snow_cover_fraction) +
        p.snow.snow_cover_fraction * snow_dest.vapor_flux
    @. csf.temp1 = ifelse(area_fraction == 0, zero(csf.temp1), csf.temp1)
    @. csf.F_turb_moisture += csf.temp1 * area_fraction

    # Combine turbulent momentum fluxes from each component of the land model
    @. csf.temp1 = soil_dest.ρτxz * (1 - p.snow.snow_cover_fraction) + p.snow.snow_cover_fraction * snow_dest.ρτxz
    @. csf.temp1 = ifelse(area_fraction == 0, zero(csf.temp1), csf.temp1)
    @. csf.F_turb_ρτxz += csf.temp1 * area_fraction

    @. csf.temp1 = soil_dest.ρτyz * (1 - p.snow.snow_cover_fraction) + p.snow.snow_cover_fraction * snow_dest.ρτyz
    @. csf.temp1 = ifelse(area_fraction == 0, zero(csf.temp1), csf.temp1)
    @. csf.F_turb_ρτyz += csf.temp1 * area_fraction
    return nothing
end
