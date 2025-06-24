import ClimaParams as CP
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
import SurfaceFluxes as SF
import SurfaceFluxes.Parameters as SFP
import Thermodynamics as TD

include("climaland_helpers.jl")

"""
    ClimaLandSimulation{M, I, A}

The integrated ClimaLand model simulation object.

It contains the following objects:
- `model::M`: The `ClimaLand.LandModel`.
- `integrator::I`: The integrator used in timestepping this model.
- `area_fraction::A`: A ClimaCore Field on the boundary space representing the surface area fraction of this component model.
- `output_writer::OW`: The diagnostic output writer.
"""
struct ClimaLandSimulation{M <: CL.LandModel, I <: SciMLBase.AbstractODEIntegrator, A <: CC.Fields.Field, OW} <:
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
        area_fraction;
        saveat::Vector{TT} = [tspan[1], tspan[2]],
        land_temperature_anomaly::String = "amip",
        use_land_diagnostics::Bool = true,
        stepper = CTS.ARS111(),
        parameter_files = [],
    )

Creates a ClimaLandSimulation object containing a land domain,
a ClimaLand.LandModel, and an integrator.

This type of model contains a canopy model, soil model, snow model, and
soil CO2 model. Specific details about the complexity of the model
can be found in the ClimaLand.jl documentation.
"""
function ClimaLandSimulation(
    ::Type{FT};
    dt::TT,
    tspan::Tuple{TT, TT},
    start_date::Dates.DateTime,
    output_dir::String,
    area_fraction,
    nelements::Tuple{Int, Int} = (101, 10),
    depth::FT = FT(50),
    dz_tuple::Tuple{FT, FT} = FT.((10.0, 0.05)),
    shared_surface_space = nothing,
    saveat::Vector{TT} = [tspan[1], tspan[2]],
    surface_elevation = nothing,
    land_temperature_anomaly::String = "amip",
    use_land_diagnostics::Bool = true,
    stepper = CTS.ARS111(),
    parameter_files = [],
) where {FT, TT <: Union{Float64, ITime}}
    # Note that this does not take into account topography of the surface, which is OK for this land model.
    # But it must be taken into account when computing surface fluxes, for Δz.
    if isnothing(shared_surface_space)
        domain = make_land_domain(depth; nelements, dz_tuple)
    else
        domain = make_land_domain(shared_surface_space, depth)
    end
    surface_space = domain.space.surface
    subsurface_space = domain.space.subsurface

    # If provided, interpolate surface elevation field to surface space; otherwise use zero elevation
    if isnothing(surface_elevation)
        surface_elevation = CC.Fields.zeros(surface_space)
    else
        surface_elevation = Interfacer.remap(surface_elevation, surface_space)
    end

    # Set up spatially-varying parameters
    earth_param_set = CL.Parameters.LandParameters(
        CP.create_toml_dict(FT; override_file = CP.merge_toml_files(parameter_files; override = true)),
    )
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
        atmos = CL.CoupledAtmosphere{FT}(surface_space),
        radiation = CL.CoupledRadiativeFluxes{FT}(
            start_date;
            insol_params = LP.insolation_parameters(earth_param_set),
            latitude = ClimaCore.Fields.coordinate_field(domain.space.surface).lat,
            longitude = ClimaCore.Fields.coordinate_field(domain.space.surface).long,
        ),
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

    # Initialize the surface temperature so the atmosphere can compute radiation.
    # This cache variable is normally computed using `p.drivers.LW_d`, which is
    #  set after the radiation calculation and exchange, but the radiation
    #  code itself requires a surface temperature as input.
    @. p.T_sfc = orog_adjusted_T_surface

    updateat = [promote(tspan[1]:10dt:(tspan[2] + dt)...)...] # add an extra time at end in case sim steps over end
    updatefunc = CL.make_update_drivers(CL.get_drivers(model))
    driver_cb = CL.DriverUpdateCallback(updateat, updatefunc)

    exp_tendency! = CL.make_exp_tendency(model)
    imp_tendency! = CL.make_imp_tendency(model)
    jacobian! = CL.make_jacobian(model)

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
            output_vars = :long,
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
        callback = SciMLBase.CallbackSet(driver_cb, diag_cb),
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

# Update fields stored in land drivers
function Interfacer.update_field!(sim::ClimaLandSimulation, ::Val{:diffuse_fraction}, field)
    Interfacer.remap!(sim.integrator.p.drivers.frac_diff, field)
end
function Interfacer.update_field!(sim::ClimaLandSimulation, ::Val{:air_temperature}, field)
    Interfacer.remap!(sim.integrator.p.drivers.T, field)
end
function Interfacer.update_field!(sim::ClimaLandSimulation, ::Val{:air_pressure}, field)
    Interfacer.remap!(sim.integrator.p.drivers.P, field)
end
function Interfacer.update_field!(sim::ClimaLandSimulation, ::Val{:air_humidity}, field)
    Interfacer.remap!(sim.integrator.p.drivers.q, field)
end
function Interfacer.update_field!(sim::ClimaLandSimulation, ::Val{:c_co2}, field)
    Interfacer.remap!(sim.integrator.p.drivers.c_co2, field)
end
function Interfacer.update_field!(sim::ClimaLandSimulation, ::Val{:liquid_precipitation}, field)
    # Arbitrarily take parameters from the soil (they are the same for all land sub-components)
    ρ_liq = (LP.ρ_cloud_liq(sim.model.soil.parameters.earth_param_set))
    Interfacer.remap!(sim.integrator.p.drivers.P_liq, field ./ ρ_liq)
end
function Interfacer.update_field!(sim::ClimaLandSimulation, ::Val{:snow_precipitation}, field)
    # Arbitrarily take parameters from the soil (they are the same for all land sub-components)
    ρ_liq = (LP.ρ_cloud_liq(sim.model.soil.parameters.earth_param_set))
    Interfacer.remap!(sim.integrator.p.drivers.P_snow, field ./ ρ_liq)
end
function Interfacer.update_field!(sim::ClimaLandSimulation, ::Val{:lw_d}, field)
    Interfacer.remap!(sim.integrator.p.drivers.LW_d, field)
end
function Interfacer.update_field!(sim::ClimaLandSimulation, ::Val{:sw_d}, field)
    Interfacer.remap!(sim.integrator.p.drivers.SW_d, field)
end

Interfacer.step!(sim::ClimaLandSimulation, t) = Interfacer.step!(sim.integrator, t - sim.integrator.t, true)
Interfacer.close_output_writers(sim::ClimaLandSimulation) = isnothing(sim.output_writer) || close(sim.output_writer)

function FieldExchanger.update_sim!(sim::ClimaLandSimulation, csf, area_fraction)
    # update fields for radiative transfer
    Interfacer.update_field!(sim, Val(:diffuse_fraction), csf.diffuse_fraction)
    Interfacer.update_field!(sim, Val(:sw_d), csf.SW_d)
    Interfacer.update_field!(sim, Val(:lw_d), csf.LW_d)

    # update fields for canopy conductance and photosynthesis
    Interfacer.update_field!(sim, Val(:c_co2), csf.c_co2)
    Interfacer.update_field!(sim, Val(:air_temperature), csf.T_atmos)
    Interfacer.update_field!(sim, Val(:air_pressure), csf.P_atmos)
    Interfacer.update_field!(sim, Val(:air_humidity), csf.q_atmos)

    # precipitation
    Interfacer.update_field!(sim, Val(:liquid_precipitation), csf.P_liq)
    Interfacer.update_field!(sim, Val(:snow_precipitation), csf.P_snow)
end

function FieldExchanger.import_atmos_fields!(csf, sim::ClimaLandSimulation, atmos_sim)
    Interfacer.get_field!(csf.diffuse_fraction, atmos_sim, Val(:diffuse_fraction))
    Interfacer.get_field!(csf.SW_d, atmos_sim, Val(:SW_d))
    Interfacer.get_field!(csf.LW_d, atmos_sim, Val(:LW_d))
    Interfacer.get_field!(csf.P_atmos, atmos_sim, Val(:air_pressure))
    Interfacer.get_field!(csf.T_atmos, atmos_sim, Val(:air_temperature))
    Interfacer.get_field!(csf.q_atmos, atmos_sim, Val(:specific_humidity))
    Interfacer.get_field!(csf.ρ_atmos, atmos_sim, Val(:air_density))
    Interfacer.get_field!(csf.P_liq, atmos_sim, Val(:liquid_precipitation))
    Interfacer.get_field!(csf.P_snow, atmos_sim, Val(:snow_precipitation))
    # CO2 is a scalar for now so it doesn't need remapping
    csf.c_co2 .= Interfacer.get_field(atmos_sim, Val(:co2))
    return nothing
end

"""
Extend Interfacer.add_coupler_fields! to add the fields required for ClimaLandSimulation.

The fields added are:
- `:SW_d` (for radiative transfer)
- `:LW_d` (for radiative transfer)
- `:diffuse_fraction` (for radiative transfer)
- `:c_co2` (for photosynthesis, biogeochemistry)
- `:P_atmos` (for canopy conductance)
- `:T_atmos` (for canopy conductance)
- `:q_atmos` (for canopy conductance)
- `P_liq` (for moisture fluxes)
- `P_snow` (for moisture fluxes)
"""
function Interfacer.add_coupler_fields!(coupler_field_names, ::ClimaLandSimulation)
    land_coupler_fields = [:SW_d, :LW_d, :diffuse_fraction, :c_co2, :P_atmos, :T_atmos, :q_atmos, :P_liq, :P_snow]
    push!(coupler_field_names, land_coupler_fields...)
end

function Checkpointer.get_model_prog_state(sim::ClimaLandSimulation)
    return sim.integrator.u
end

function Checkpointer.get_model_cache(sim::ClimaLandSimulation)
    return sim.integrator.p
end

function Checkpointer.restore_cache!(sim::ClimaLandSimulation, new_cache)
    old_cache = Checkpointer.get_model_cache(sim)
    comms_ctx = ClimaComms.context(sim.model.soil)
    restore!(
        old_cache,
        new_cache,
        comms_ctx,
        ignore = Set([:dss_buffer_2d, :dss_buffer_3d, :scratch1, :scratch2, :scratch3, :sfc_scratch, :subsfc_scratch]),
    )
end

## Extend functions for land-specific flux calculation
"""
    compute_surface_fluxes!(csf, sim::ClimaLandSimulation, atmos_sim, thermo_params)

This function computes surface fluxes between the integrated land model
simulation and the atmosphere.

Update the input coupler surface fields `csf` in-place with the computed fluxes
for this model. These are then summed using area-weighting across all surface
models to get the total fluxes. Fluxes where the area fraction is zero are set to zero.

Because the integrated land model is composed of multiple sub-components, the
fluxes are computed for each sub-component and then combined to get the total for this model.
The land model cache is updated with the computed fluxes for each sub-component.

Currently, this calculation is done on the land surface space, and the computed fluxes
are remapped onto the coupler boundary space as the coupler fields are updated. Ideally,
we would compute fluxes on the coupler boundary space directly (as we do for other components),
but this is not done currently because the `coupler_compute_turbulent_fluxes!` functions
internally use land variables that are defined on the land surface space.

# Arguments
- `csf`: [CC.Fields.Field] containing a NamedTuple of turbulent flux fields: `F_turb_ρτxz`, `F_turb_ρτyz`, `F_lh`, `F_sh`, `F_turb_moisture`.
- `sim`: [ClimaLandSimulation] the integrated land simulation to compute fluxes for.
- `atmos_sim`: [Interfacer.AtmosModelSimulation] the atmosphere simulation to compute fluxes with.
- `thermo_params`: [ClimaParams.ThermodynamicParameters] the thermodynamic parameters for the simulation.
"""
function FluxCalculator.compute_surface_fluxes!(
    csf,
    sim::ClimaLandSimulation,
    atmos_sim::Interfacer.AtmosModelSimulation,
    thermo_params,
)
    boundary_space = axes(csf)
    FT = CC.Spaces.undertype(boundary_space)
    Y, p, t, model = sim.integrator.u, sim.integrator.p, sim.integrator.t, sim.model

    # We should change this to be on the boundary_space
    land_space = axes(p.soil.turbulent_fluxes)
    coupled_atmos = sim.model.soil.boundary_conditions.top.atmos

    # Update the land simulation's coupled atmosphere state
    Interfacer.get_field!(coupled_atmos.h, atmos_sim, Val(:height_int))

    # Use scratch space for remapped wind vector components to avoid allocations
    Interfacer.get_field!(p.scratch1, atmos_sim, Val(:u_int)) # u_atmos
    Interfacer.get_field!(p.scratch2, atmos_sim, Val(:v_int)) # v_atmos
    @. coupled_atmos.u = StaticArrays.SVector(p.scratch1, p.scratch2)

    # Use scratch space for remapped atmospheric fields to avoid allocations
    Interfacer.remap!(p.scratch1, csf.ρ_atmos)
    Interfacer.remap!(p.scratch2, csf.T_atmos)
    Interfacer.remap!(p.scratch3, csf.q_atmos)
    @. coupled_atmos.thermal_state = TD.PhaseEquil_ρTq(thermo_params, p.scratch1, p.scratch2, p.scratch3)

    # set the same atmosphere state for all sub-components
    @assert sim.model.soil.boundary_conditions.top.atmos ===
            sim.model.canopy.boundary_conditions.atmos ===
            sim.model.snow.boundary_conditions.atmos ===
            coupled_atmos

    # compute the fluxes for each sub-component and update the land model cache
    soil_dest = p.soil.turbulent_fluxes
    CL.coupler_compute_turbulent_fluxes!(soil_dest, coupled_atmos, model.soil, Y, p, t)

    snow_dest = p.snow.turbulent_fluxes
    CL.coupler_compute_turbulent_fluxes!(snow_dest, coupled_atmos, model.snow, Y, p, t)

    canopy_dest = p.canopy.turbulent_fluxes
    CL.coupler_compute_turbulent_fluxes!(canopy_dest, coupled_atmos, model.canopy, Y, p, t)

    # Get area fraction of the land model (min = 0, max = 1)
    area_fraction = Interfacer.get_field(sim, Val(:area_fraction))

    # Combine turbulent energy fluxes from each component of the land model
    # Use temporary variables to avoid allocating
    Interfacer.remap!(
        csf.temp1,
        canopy_dest.lhf .+ soil_dest.lhf .* (1 .- p.snow.snow_cover_fraction) .+
        p.snow.snow_cover_fraction .* snow_dest.lhf,
    )
    Interfacer.remap!(
        csf.temp2,
        canopy_dest.shf .+ soil_dest.shf .* (1 .- p.snow.snow_cover_fraction) .+
        p.snow.snow_cover_fraction .* snow_dest.shf,
    )

    # Zero out the fluxes where the area fraction is zero
    @. csf.temp1 = ifelse(area_fraction == 0, zero(csf.temp1), csf.temp1)
    @. csf.temp2 = ifelse(area_fraction == 0, zero(csf.temp2), csf.temp2)

    # Update the coupler field in-place
    @. csf.F_lh += csf.temp1 * area_fraction
    @. csf.F_sh += csf.temp2 * area_fraction

    # Combine turbulent moisture fluxes from each component of the land model
    Interfacer.remap!(
        csf.temp1,
        canopy_dest.transpiration .+
        (soil_dest.vapor_flux_liq .+ soil_dest.vapor_flux_ice) .* (1 .- p.snow.snow_cover_fraction) .+
        p.snow.snow_cover_fraction .* snow_dest.vapor_flux,
    )
    @. csf.temp1 = ifelse(area_fraction == 0, zero(csf.temp1), csf.temp1)
    @. csf.F_turb_moisture += csf.temp1 * area_fraction

    # Combine turbulent momentum fluxes from each component of the land model
    # Note that we exclude the canopy component here for now, since we can have nonzero momentum fluxes
    #  where there is zero LAI. This should be fixed in ClimaLand.
    Interfacer.remap!(
        csf.temp1,
        soil_dest.ρτxz .* (1 .- p.snow.snow_cover_fraction) .+ p.snow.snow_cover_fraction .* snow_dest.ρτxz,
    )
    @. csf.temp1 = ifelse(area_fraction == 0, zero(csf.temp1), csf.temp1)
    @. csf.F_turb_ρτxz += csf.temp1 * area_fraction

    Interfacer.remap!(
        csf.temp1,
        soil_dest.ρτyz .* (1 .- p.snow.snow_cover_fraction) .+ p.snow.snow_cover_fraction .* snow_dest.ρτyz,
    )
    @. csf.temp1 = ifelse(area_fraction == 0, zero(csf.temp1), csf.temp1)
    @. csf.F_turb_ρτyz += csf.temp1 * area_fraction

    # Combine the buoyancy flux from each component of the land model
    # Note that we exclude the canopy component here for now, since ClimaLand doesn't
    #  include its extra resistance term in the buoyancy flux calculation.
    Interfacer.remap!(
        csf.temp1,
        soil_dest.buoy_flux .* (1 .- p.snow.snow_cover_fraction) .+ p.snow.snow_cover_fraction .* snow_dest.buoy_flux,
    )
    @. csf.temp1 = ifelse(area_fraction == 0, zero(csf.temp1), csf.temp1)
    @. csf.buoyancy_flux += csf.temp1 * area_fraction

    # Compute ustar from the momentum fluxes and surface air density
    #  ustar = sqrt(ρτ / ρ)
    @. csf.temp1 = sqrt(sqrt(csf.F_turb_ρτxz^2 + csf.F_turb_ρτyz^2) / csf.ρ_atmos)
    @. csf.temp1 = ifelse(area_fraction == 0, zero(csf.temp1), csf.temp1)
    # If ustar is zero, set it to eps to avoid division by zero in the atmosphere
    @. csf.ustar += max(csf.temp1 * area_fraction, eps(FT))

    # Compute the Monin-Obukhov length from ustar and the buoyancy flux
    #  L_MO = -u^3 / (k * buoyancy_flux)
    # Prevent dividing by zero in the case of zero buoyancy flux
    function non_zero(v::FT) where {FT}
        sign_of_v = v == 0 ? 1 : sign(v)
        return abs(v) < eps(FT) ? eps(FT) * sign_of_v : v
    end
    surface_params = LP.surface_fluxes_parameters(sim.model.soil.parameters.earth_param_set)
    @. csf.temp1 = -csf.ustar^3 / SFP.von_karman_const(surface_params) / non_zero(csf.buoyancy_flux)
    @. csf.temp1 = ifelse(area_fraction == 0, zero(csf.temp1), csf.temp1)
    # When L_MO is infinite, avoid multiplication by zero to prevent NaN
    @. csf.L_MO += ifelse(isinf(csf.temp1), csf.temp1, csf.temp1 * area_fraction)

    return nothing
end

"""
    Interfacer.set_cache!(sim::ClimaLandSimulation)

Set cache variables that cannot be initialized before the initial exchange.
This must be called after radiation, so that `p.drivers`
is filled with the initial radiation fluxes, and these can be propagated
to the rest of the cache (e.g. in canopy radative transfer).

This function does not set all the cache variables, because many are computed
as part of the tendendencies.
"""
function Interfacer.set_cache!(sim::ClimaLandSimulation)
    land_set_initial_cache! = CL.make_set_initial_cache(sim.model)
    land_set_initial_cache!(sim.integrator.p, sim.integrator.u, sim.integrator.t)
    return nothing
end
