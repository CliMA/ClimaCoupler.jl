import ClimaParams as CP
import ClimaLand as CL
import ClimaLand.Parameters as LP
import Dates
import ClimaUtilities.TimeVaryingInputs:
    LinearInterpolation, PeriodicCalendar, TimeVaryingInput
import ClimaUtilities.SpaceVaryingInputs: SpaceVaryingInput
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
struct ClimaLandSimulation{
    M <: CL.LandModel,
    I <: SciMLBase.AbstractODEIntegrator,
    A <: CC.Fields.Field,
    OW,
} <: Interfacer.LandModelSimulation
    model::M
    integrator::I
    area_fraction::A
    output_writer::OW
end

"""
    ClimaLandSimulation(
        ::Type{FT};
        dt::TT,
        tspan::Tuple{TT, TT},
        start_date::Dates.DateTime,
        output_dir::String,
        area_fraction,
        nelements::Tuple{Int, Int} = (101, 15),
        depth::FT = FT(50),
        dz_tuple::Tuple{FT, FT} = FT.((10.0, 0.05)),
        shared_surface_space = nothing,
        land_spun_up_ic::Bool = true,
        saveat::Vector{TT} = [tspan[1], tspan[2]],
        surface_elevation = nothing,
        atmos_h,
        land_temperature_anomaly::String = "amip",
        use_land_diagnostics::Bool = true,
        parameter_files = [],
        land_ic_path::Union{Nothing,String} = nothing,
    ) where {FT, TT <: Union{Float64, ITime}}

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
    nelements::Tuple{Int, Int} = (101, 15),
    depth::FT = FT(50),
    dz_tuple::Tuple{FT, FT} = FT.((10.0, 0.05)),
    shared_surface_space = nothing,
    land_spun_up_ic::Bool = true,
    saveat::Vector{TT} = [tspan[1], tspan[2]],
    surface_elevation = nothing,
    atmos_h,
    land_temperature_anomaly::String = "amip",
    use_land_diagnostics::Bool = true,
    parameter_files = [],
    land_ic_path::Union{Nothing, String} = nothing,
) where {FT, TT <: Union{Float64, ITime}}
    # Note that this does not take into account topography of the surface, which is OK for this land model.
    # But it must be taken into account when computing surface fluxes, for Δz.


    if isnothing(shared_surface_space)
        domain = make_land_domain(depth; nelements, dz_tuple)
    else
        domain = make_land_domain(
            shared_surface_space,
            depth;
            nelements_vert = nelements[2],
            dz_tuple,
        )
    end
    surface_space = domain.space.surface
    subsurface_space = domain.space.subsurface

    # If provided, interpolate surface elevation field to surface space; otherwise use zero elevation
    if isnothing(surface_elevation)
        surface_elevation = CC.Fields.zeros(surface_space)
    else
        surface_elevation = Interfacer.remap(surface_elevation, surface_space)
    end
    # Interpolate atmosphere height field to surface space of land model,
    #  since that's where we compute fluxes for this land model
    atmos_h = Interfacer.remap(atmos_h, surface_space)

    # Set up spatially-varying parameters
    default_parameter_file = joinpath(pkgdir(CL), "toml", "default_parameters.toml")
    toml_dict = CP.create_toml_dict(
        FT;
        override_file = CP.merge_toml_files(
            [default_parameter_file, parameter_files...];
            override = true,
        ),
    )
    earth_param_set = CL.Parameters.LandParameters(toml_dict)

    # Set up atmosphere and radiation forcing
    forcing = (;
        atmos = CL.CoupledAtmosphere{FT}(surface_space, atmos_h),
        radiation = CL.CoupledRadiativeFluxes{FT}(
            start_date;
            insol_params = LP.insolation_parameters(earth_param_set),
            latitude = ClimaCore.Fields.coordinate_field(domain.space.surface).lat,
            longitude = ClimaCore.Fields.coordinate_field(domain.space.surface).long,
        ),
    )

    # Set up leaf area index (LAI)
    stop_date = start_date + Dates.Second(float(tspan[2] - tspan[1]))
    LAI = CL.prescribed_lai_modis(
        surface_space,
        start_date,
        stop_date;
        time_interpolation_method = LinearInterpolation(),
    )

    model = CL.LandModel{FT}(forcing, LAI, toml_dict, domain, dt)

    Y, p, coords = CL.initialize(model)

    # Set initial conditions

    # Apply temperature anomaly function to initial temperature only if specified
    T_base = FT(276.85)
    if land_temperature_anomaly != "nothing"
        T_functions =
            Dict("aquaplanet" => temp_anomaly_aquaplanet, "amip" => temp_anomaly_amip)
        haskey(T_functions, land_temperature_anomaly) ||
            error("land temp anomaly function $land_temperature_anomaly not supported")
        temp_anomaly = T_functions[land_temperature_anomaly]
        T_sfc0 = T_base .+ temp_anomaly.(coords.subsurface)
    else
        # constant field on subsurface space
        T_sfc0 = CC.Fields.Field(T_base .* CC.Fields.ones(subsurface_space))
    end
    lapse_rate = FT(6.5e-3)
    # Adjust initial temperature to account for orography of the surface
    # `surface_elevation` is a ClimaCore.Fields.Field(`half` level)
    orog_adjusted_T_data =
        CC.Fields.field_values(T_sfc0) .-
        lapse_rate .* CC.Fields.field_values(surface_elevation)
    orog_adjusted_T = CC.Fields.Field(orog_adjusted_T_data, subsurface_space)
    orog_adjusted_T_surface =
        CC.Fields.Field(CC.Fields.level(orog_adjusted_T_data, 1), surface_space)

    # Set initial conditions that aren't read in from file
    Y.soilco2.C .= FT(0.000412) # set to atmospheric co2, mol co2 per mol air
    Y.canopy.hydraulics.ϑ_l.:1 .= model.canopy.hydraulics.parameters.ν
    @. Y.canopy.energy.T = orog_adjusted_T_surface

    # Read in initial conditions for snow and soil from file, if requested
    (; θ_r, ν, ρc_ds) = model.soil.parameters


    if !land_spun_up_ic && !isnothing(land_ic_path)
        # Subseasonal setup: land_spun_up_ic false, but a land_ic_path is provided
        ic_path = land_ic_path
        @info "ClimaLand: using land IC file" ic_path

        regridder_type = :InterpolationsRegridder
        extrapolation_bc =
            (Interpolations.Periodic(), Interpolations.Flat(), Interpolations.Flat())
        interpolation_method = Interpolations.Linear()

        # Set snow T first to use in computing snow internal energy from IC file
        p.snow.T .= SpaceVaryingInput(
            ic_path,
            "tsn",
            surface_space;
            regridder_type,
            regridder_kwargs = (; extrapolation_bc, interpolation_method),
        )
        # Set canopy temperature to skin temperature
        Y.canopy.energy.T .= SpaceVaryingInput(
            ic_path,
            "skt",
            surface_space;
            regridder_type,
            regridder_kwargs = (; extrapolation_bc, interpolation_method),
        )

        # Initialize the surface temperature so the atmosphere can compute radiation.
        # Set surface temperature to skin temperature
        p.T_sfc .= SpaceVaryingInput(
            ic_path,
            "skt",
            surface_space;
            regridder_type,
            regridder_kwargs = (; extrapolation_bc, interpolation_method),
        )

        CL.Simulations.set_snow_initial_conditions!(
            Y,
            p,
            surface_space,
            ic_path,
            model.snow.parameters,
        )

        T_bounds = extrema(Y.canopy.energy.T)
        CL.Simulations.set_soil_initial_conditions!(
            Y,
            ν,
            θ_r,
            subsurface_space,
            ic_path,
            model.soil,
            T_bounds,
        )
    elseif land_spun_up_ic
        # Use artifact spun-up initial conditions
        ic_path = CL.Artifacts.soil_ic_2008_50m_path()
        @info "ClimaLand: using land IC file" ic_path

        # Set snow T to orography-adjusted surface temperature before computing internal energy
        p.snow.T .= orog_adjusted_T_surface

        CL.Simulations.set_snow_initial_conditions!(
            Y,
            p,
            surface_space,
            ic_path,
            model.snow.parameters,
        )

        T_bounds = extrema(Y.canopy.energy.T)
        CL.Simulations.set_soil_initial_conditions!(
            Y,
            ν,
            θ_r,
            subsurface_space,
            ic_path,
            model.soil,
            T_bounds,
        )

        # Initialize the surface temperature so the atmosphere can compute radiation.
        @. p.T_sfc = orog_adjusted_T_surface
    else
        # Set initial conditions for the state
        @. Y.soil.ϑ_l = θ_r + (ν - θ_r) / 2
        Y.soil.θ_i .= FT(0.0)
        ρc_s =
            CL.Soil.volumetric_heat_capacity.(
                Y.soil.ϑ_l,
                Y.soil.θ_i,
                ρc_ds,
                earth_param_set,
            )
        Y.soil.ρe_int .=
            CL.Soil.volumetric_internal_energy.(
                Y.soil.θ_i,
                ρc_s,
                orog_adjusted_T,
                earth_param_set,
            )

        Y.snow.S .= FT(0)
        Y.snow.S_l .= FT(0)
        Y.snow.U .= FT(0)

        # Initialize the surface temperature so the atmosphere can compute radiation.
        @. p.T_sfc = orog_adjusted_T_surface
    end
    # Initialize the surface emissivity so the atmosphere can compute radiation.
    # Otherwise, it's initialized to 0 which causes NaNs in the radiation calculation.
    @. p.ϵ_sfc = FT(1)

    # Update cos(zenith angle) within land model every hour
    update_dt = dt isa ITime ? ITime(3600) : 3600
    updateat = [promote(tspan[1]:update_dt:(tspan[2] + dt)...)...] # add an extra time at end in case sim steps over end
    updatefunc = CL.make_update_drivers(CL.get_drivers(model))
    driver_cb = CL.DriverUpdateCallback(updateat, updatefunc)

    exp_tendency! = CL.make_exp_tendency(model)
    imp_tendency! = CL.make_imp_tendency(model)
    jacobian! = CL.make_jacobian(model)

    # set up jacobian info
    jac_kwargs = (; jac_prototype = CL.FieldMatrixWithSolver(Y), Wfact = jacobian!)

    prob = SciMLBase.ODEProblem(
        CTS.ClimaODEFunction(
            T_exp! = exp_tendency!,
            T_imp! = SciMLBase.ODEFunction(imp_tendency!; jac_kwargs...),
        ),
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
            reduction_period = :monthly,
        )
        diagnostic_handler =
            CD.DiagnosticsHandler(scheduled_diagnostics, Y, p, tspan[1]; dt = dt)
        diag_cb = CD.DiagnosticsCallback(diagnostic_handler)
    else
        output_writer = nothing
        diag_cb = nothing
    end

    # Set up time stepper and integrator
    stepper = CTS.ARS111()
    ode_algo = CTS.IMEXAlgorithm(
        stepper,
        CTS.NewtonsMethod(
            max_iters = 3,
            update_j = CTS.UpdateEvery(CTS.NewNewtonIteration),
        ),
    )
    integrator = SciMLBase.init(
        prob,
        ode_algo;
        dt,
        saveat,
        adaptive = false,
        callback = SciMLBase.CallbackSet(driver_cb, diag_cb),
    )

    return ClimaLandSimulation(model, integrator, area_fraction, output_writer)
end

###############################################################################
### Functions required by ClimaCoupler.jl for a SurfaceModelSimulation
###############################################################################

Interfacer.get_field(sim::ClimaLandSimulation, ::Val{:area_fraction}) = sim.area_fraction
Interfacer.get_field(sim::ClimaLandSimulation, ::Val{:beta}) =
    CL.surface_evaporative_scaling(sim.model, sim.integrator.u, sim.integrator.p)
Interfacer.get_field(sim::ClimaLandSimulation, ::Val{:emissivity}) = sim.integrator.p.ϵ_sfc
Interfacer.get_field(sim::ClimaLandSimulation, ::Val{:energy}) =
    CL.total_energy(sim.integrator.u, sim.integrator.p)
Interfacer.get_field(sim::ClimaLandSimulation, ::Val{:surface_direct_albedo}) =
    CL.surface_albedo(sim.model, sim.integrator.u, sim.integrator.p)
Interfacer.get_field(sim::ClimaLandSimulation, ::Val{:surface_diffuse_albedo}) =
    CL.surface_albedo(sim.model, sim.integrator.u, sim.integrator.p)
Interfacer.get_field(sim::ClimaLandSimulation, ::Val{:water}) =
    CL.total_water(sim.integrator.u, sim.integrator.p)
Interfacer.get_field(sim::ClimaLandSimulation, ::Val{:surface_temperature}) =
    sim.integrator.p.T_sfc

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
function Interfacer.update_field!(
    sim::ClimaLandSimulation,
    ::Val{:liquid_precipitation},
    field,
)
    # Arbitrarily take parameters from the soil (they are the same for all land sub-components)
    ρ_liq = (LP.ρ_cloud_liq(sim.model.soil.parameters.earth_param_set))
    Interfacer.remap!(sim.integrator.p.drivers.P_liq, field ./ ρ_liq)
end
function Interfacer.update_field!(
    sim::ClimaLandSimulation,
    ::Val{:snow_precipitation},
    field,
)
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

function Interfacer.step!(sim::ClimaLandSimulation, t)
    while float(sim.integrator.t) < float(t)
        Interfacer.step!(sim.integrator)
    end
    return nothing
end
Interfacer.close_output_writers(sim::ClimaLandSimulation) =
    isnothing(sim.output_writer) || close(sim.output_writer)

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

"""
   FieldExchanger.import_atmos_fields!(csf, sim::ClimaLandSimulation, atmos_sim)

Import non-default coupler fields from the atmosphere simulation into the coupler fields.
These include the diffuse fraction of light, shortwave and longwave downwelling radiation,
air pressure, and CO2 concentration.

The default coupler fields will be imported by the default method implemented in
FieldExchanger.jl.
"""
function FieldExchanger.import_atmos_fields!(csf, sim::ClimaLandSimulation, atmos_sim)
    Interfacer.get_field!(csf.diffuse_fraction, atmos_sim, Val(:diffuse_fraction))
    Interfacer.get_field!(csf.SW_d, atmos_sim, Val(:SW_d))
    Interfacer.get_field!(csf.LW_d, atmos_sim, Val(:LW_d))
    Interfacer.get_field!(csf.P_atmos, atmos_sim, Val(:air_pressure))
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
    land_coupler_fields = [
        :SW_d,
        :LW_d,
        :diffuse_fraction,
        :c_co2,
        :P_atmos,
        :T_atmos,
        :q_atmos,
        :P_liq,
        :P_snow,
    ]
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
        ignore = Set([
            :dss_buffer_2d,
            :dss_buffer_3d,
            :scratch1,
            :scratch2,
            :scratch3,
            :sfc_scratch,
            :subsfc_scratch,
        ]),
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
    @. coupled_atmos.thermal_state =
        TD.PhaseEquil_ρTq(thermo_params, p.scratch1, p.scratch2, p.scratch3)

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
        csf.scalar_temp1,
        canopy_dest.lhf .+ soil_dest.lhf .* (1 .- p.snow.snow_cover_fraction) .+
        p.snow.snow_cover_fraction .* snow_dest.lhf,
    )
    Interfacer.remap!(
        csf.scalar_temp2,
        canopy_dest.shf .+ soil_dest.shf .* (1 .- p.snow.snow_cover_fraction) .+
        p.snow.snow_cover_fraction .* snow_dest.shf,
    )

    # Zero out the fluxes where the area fraction is zero
    @. csf.scalar_temp1 =
        ifelse(area_fraction == 0, zero(csf.scalar_temp1), csf.scalar_temp1)
    @. csf.scalar_temp2 =
        ifelse(area_fraction == 0, zero(csf.scalar_temp2), csf.scalar_temp2)

    # Update the coupler field in-place
    @. csf.F_lh += csf.scalar_temp1 * area_fraction
    @. csf.F_sh += csf.scalar_temp2 * area_fraction

    # Combine turbulent moisture fluxes from each component of the land model
    # Note that we multiply by ρ_liq to convert from m s-1 to kg m-2 s-1
    ρ_liq = (LP.ρ_cloud_liq(sim.model.soil.parameters.earth_param_set))
    Interfacer.remap!(
        csf.scalar_temp1,
        (
            canopy_dest.transpiration .+
            (soil_dest.vapor_flux_liq .+ soil_dest.vapor_flux_ice) .*
            (1 .- p.snow.snow_cover_fraction) .+
            p.snow.snow_cover_fraction .* snow_dest.vapor_flux
        ) .* ρ_liq,
    )
    @. csf.scalar_temp1 =
        ifelse(area_fraction == 0, zero(csf.scalar_temp1), csf.scalar_temp1)
    @. csf.F_turb_moisture += csf.scalar_temp1 * area_fraction

    # Combine turbulent momentum fluxes from each component of the land model
    # Note that we exclude the canopy component here for now, since we can have nonzero momentum fluxes
    #  where there is zero LAI. This should be fixed in ClimaLand.
    Interfacer.remap!(
        csf.scalar_temp1,
        soil_dest.ρτxz .* (1 .- p.snow.snow_cover_fraction) .+
        p.snow.snow_cover_fraction .* snow_dest.ρτxz,
    )
    @. csf.scalar_temp1 =
        ifelse(area_fraction == 0, zero(csf.scalar_temp1), csf.scalar_temp1)
    @. csf.F_turb_ρτxz += csf.scalar_temp1 * area_fraction

    Interfacer.remap!(
        csf.scalar_temp1,
        soil_dest.ρτyz .* (1 .- p.snow.snow_cover_fraction) .+
        p.snow.snow_cover_fraction .* snow_dest.ρτyz,
    )
    @. csf.scalar_temp1 =
        ifelse(area_fraction == 0, zero(csf.scalar_temp1), csf.scalar_temp1)
    @. csf.F_turb_ρτyz += csf.scalar_temp1 * area_fraction

    # Combine the buoyancy flux from each component of the land model
    # Note that we exclude the canopy component here for now, since ClimaLand doesn't
    #  include its extra resistance term in the buoyancy flux calculation.
    Interfacer.remap!(
        csf.scalar_temp1,
        soil_dest.buoy_flux .* (1 .- p.snow.snow_cover_fraction) .+
        p.snow.snow_cover_fraction .* snow_dest.buoy_flux,
    )
    @. csf.scalar_temp1 =
        ifelse(area_fraction == 0, zero(csf.scalar_temp1), csf.scalar_temp1)
    @. csf.buoyancy_flux += csf.scalar_temp1 * area_fraction

    # Compute ustar from the momentum fluxes and surface air density
    #  ustar = sqrt(ρτ / ρ)
    @. csf.scalar_temp1 = sqrt(sqrt(csf.F_turb_ρτxz^2 + csf.F_turb_ρτyz^2) / csf.ρ_atmos)
    @. csf.scalar_temp1 =
        ifelse(area_fraction == 0, zero(csf.scalar_temp1), csf.scalar_temp1)
    # If ustar is zero, set it to eps to avoid division by zero in the atmosphere
    @. csf.ustar += max(csf.scalar_temp1 * area_fraction, eps(FT))

    # Compute the Monin-Obukhov length from ustar and the buoyancy flux
    #  L_MO = -u^3 / (k * buoyancy_flux)
    # Prevent dividing by zero in the case of zero buoyancy flux
    function non_zero(v::FT) where {FT}
        sign_of_v = v == 0 ? 1 : sign(v)
        return abs(v) < eps(FT) ? eps(FT) * sign_of_v : v
    end
    surface_params = LP.surface_fluxes_parameters(sim.model.soil.parameters.earth_param_set)
    @. csf.scalar_temp1 =
        -csf.ustar^3 / SFP.von_karman_const(surface_params) / non_zero(csf.buoyancy_flux)
    @. csf.scalar_temp1 =
        ifelse(area_fraction == 0, zero(csf.scalar_temp1), csf.scalar_temp1)
    # When L_MO is infinite, avoid multiplication by zero to prevent NaN
    @. csf.L_MO +=
        ifelse(isinf(csf.scalar_temp1), csf.scalar_temp1, csf.scalar_temp1 * area_fraction)

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
