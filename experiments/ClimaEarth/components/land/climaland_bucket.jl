import Dates
import SciMLBase
import Statistics
import ClimaComms
import ClimaCore as CC
import ClimaTimeSteppers as CTS
import Thermodynamics as TD
import ClimaLand as CL
import ClimaLand.Parameters as LP
import ClimaParams as CP
import ClimaDiagnostics as CD
import ClimaCoupler: Checkpointer, FluxCalculator, Interfacer, FieldExchanger, Plotting
using NCDatasets
include("climaland_helpers.jl")


###
### Functions required by ClimaCoupler.jl for a AbstractSurfaceSimulation
###
"""
    BucketSimulation{M, I, A}

The bucket model simulation object.

It contains the following objects:
- `model::M`: The `ClimaLand.Bucket.BucketModel`.
- `integrator::I`: The integrator used in timestepping this model.
- `area_fraction::A`: A ClimaCore Field on the boundary space representing the surface area fraction of this component model.
- `output_writer`: The diagnostic file writer.
"""
struct BucketSimulation{
    M <: CL.Bucket.BucketModel,
    I <: SciMLBase.AbstractODEIntegrator,
    A <: CC.Fields.Field,
    OW,
    RF,
} <: Interfacer.AbstractLandSimulation
    model::M
    integrator::I
    area_fraction::A
    output_writer::OW
    radiative_fluxes::RF
end

"""
    Interfacer.LandSimulation(::Type{FT}, ::Val{:bucket}; kwargs...)

Extension of the generic LandSimulation constructor for the bucket model.
"""
function Interfacer.LandSimulation(::Type{FT}, ::Val{:bucket}; kwargs...) where {FT}
    return BucketSimulation(FT; kwargs...)
end

"""
    BucketSimulation

Initializes the bucket model variables.
"""
function BucketSimulation(
    ::Type{FT};
    dt::TT,
    tspan::Tuple{TT, TT},
    start_date::Dates.DateTime,
    output_dir::String,
    area_fraction,
    nelements::Tuple{Int, Int} = (50, 10),
    depth::FT = FT(3.5),
    dz_tuple::Tuple{FT, FT} = FT.((1, 0.05)),
    shared_surface_space = nothing,
    saveat::Vector{TT} = [tspan[1], tspan[2]],
    surface_elevation = nothing,
    atmos_h,
    land_temperature_anomaly::String = "amip",
    use_land_diagnostics::Bool = true,
    stepper = CTS.RK4(),
    albedo_type::String = "map_static",
    bucket_initial_condition::String = "",
    coupled_param_dict = CP.create_toml_dict(FT),
    extra_kwargs...,
) where {FT, TT <: Union{Float64, ITime}}
    # Get default land parameters from ClimaLand.LandParameters
    land_toml_dict = LP.create_toml_dict(FT)
    # Override land parameters with coupled parameters
    toml_dict = CP.merge_override_default_values(coupled_param_dict, land_toml_dict)

    # Note that this does not take into account topography of the surface, which is OK for this land model.
    # But it must be taken into account when computing surface fluxes, for Δz.
    if isnothing(shared_surface_space)
        domain = make_land_domain(depth, toml_dict; nelements, dz_tuple)
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
        surface_elevation = Interfacer.remap(surface_space, surface_elevation)
    end

    if albedo_type == "map_static" # Read in albedo from static data file (default type)
        # By default, this uses a file containing bareground albedo without a time component. Snow albedo is specified separately.
        albedo = CL.Bucket.PrescribedBaregroundAlbedo(toml_dict, surface_space)
    elseif albedo_type == "map_temporal" # Read in albedo from data file containing data over time
        # By default, this uses a file containing linearly-interpolated monthly data of clear-sky albedo, generated from CERES.
        albedo = CL.Bucket.PrescribedSurfaceAlbedo{FT}(
            start_date,
            surface_space;
            albedo_file_path = CL.Artifacts.ceres_albedo_dataset_path(),
            varname = "sw_alb_clr",
        )
    elseif albedo_type == "function" # Use prescribed function of lat/lon for surface albedo
        function α_bareground(coordinate_point)
            (; lat, long) = coordinate_point
            return typeof(lat)(0.38)
        end
        α_snow = toml_dict["alpha_snow"] # snow albedo
        albedo =
            CL.Bucket.PrescribedBaregroundAlbedo{FT}(α_snow, α_bareground, surface_space)
    else
        error("invalid albedo type $albedo_type")
    end

    # This is the timescale on which snow exponentially damps to zero, in the case where all
    # the snow would melt in time `τc`. It prevents us from having to specially time step in cases where
    # all the snow melts in a single timestep.
    τc = FT(float(dt))
    params = CL.Bucket.BucketModelParameters(toml_dict; albedo, τc)

    # Interpolate atmosphere height field to surface space of land model,
    #  since that's where we compute fluxes for this land model
    atmos_h = Interfacer.remap(surface_space, atmos_h)

    args = (
        params,
        CL.CoupledAtmosphere{FT}(surface_space, atmos_h),
        CL.CoupledRadiativeFluxes{FT}(),
        domain,
    )
    model = CL.Bucket.BucketModel{FT, typeof.(args)...}(args...)

    # Initial conditions with no moisture
    Y, p, coords = CL.initialize(model)

    if land_temperature_anomaly != "nothing"
        T_functions =
            Dict("aquaplanet" => temp_anomaly_aquaplanet, "amip" => temp_anomaly_amip)
        haskey(T_functions, land_temperature_anomaly) ||
            error("land temp anomaly function $land_temperature_anomaly not supported")
        temp_anomaly = T_functions[land_temperature_anomaly]

        # Set temperature IC including anomaly, based on atmospheric setup
        # Bucket surface temperature is in `p.bucket.T_sfc` (ClimaLand.jl)
        lapse_rate = FT(6.5e-3)
        T_sfc_0 = FT(271)
        @. Y.bucket.T = T_sfc_0 + temp_anomaly(coords.subsurface)
        # `surface_elevation` is a ClimaCore.Fields.Field(`half` level)
        orog_adjusted_T_data =
            CC.Fields.field_values(Y.bucket.T) .-
            lapse_rate .* CC.Fields.field_values(surface_elevation)
        orog_adjusted_T = CC.Fields.Field(orog_adjusted_T_data, domain.space.subsurface)
        # Adjust T based on surface elevation (p.bucket.T_sfc is then set using the
        # set_initial_cache! function)
        Y.bucket.T .= orog_adjusted_T
    end

    Y.bucket.W .= 0.15
    Y.bucket.Ws .= 0.0
    Y.bucket.σS .= 0.0

    # Overwrite initial conditions with interpolated values from a netcdf file using
    # the `SpaceVaryingInputs` tool. We expect the file to contain the following variables:
    # - `W`, for subsurface water storage (2D),
    # - `Ws`, for surface water content (2D),
    # - `T`, for soil temperature (3D),
    # - `S`, for snow water equivalent (2D).

    if !isempty(bucket_initial_condition)
        @info "ClimaLand Bucket using land IC file" bucket_initial_condition
        ds = NCDataset(bucket_initial_condition)
        has_all_variables = all(key -> haskey(ds, key), ["W", "Ws", "T", "S"])
        @assert has_all_variables "The land iniital condition file is expected to contain the variables W, Ws, T, and S (read documentation about requirements)."
        close(ds)

        surface_space = domain.space.surface
        subsurface_space = domain.space.subsurface
        regridder_type = :InterpolationsRegridder
        extrapolation_bc =
            (Interpolations.Periodic(), Interpolations.Flat(), Interpolations.Flat())

        Y.bucket.W .= SpaceVaryingInput(
            bucket_initial_condition,
            "W",
            surface_space;
            regridder_type,
            regridder_kwargs = (; extrapolation_bc,),
        )
        Y.bucket.Ws .= SpaceVaryingInput(
            bucket_initial_condition,
            "Ws",
            surface_space;
            regridder_type,
            regridder_kwargs = (; extrapolation_bc,),
        )
        Y.bucket.T .= SpaceVaryingInput(
            bucket_initial_condition,
            "T",
            subsurface_space;
            regridder_type,
            regridder_kwargs = (; extrapolation_bc,),
        )
        Y.bucket.σS .= SpaceVaryingInput(
            bucket_initial_condition,
            "S",
            surface_space;
            regridder_type,
            regridder_kwargs = (; extrapolation_bc,),
        )
    end

    # Set initial aux variable values
    set_initial_cache! = CL.make_set_initial_cache(model)
    set_initial_cache!(p, Y, tspan[1])

    exp_tendency! = CL.make_exp_tendency(model)
    ode_algo = CTS.ExplicitAlgorithm(stepper)
    bucket_ode_function = CTS.ClimaODEFunction(T_exp! = exp_tendency!)
    prob = SciMLBase.ODEProblem(bucket_ode_function, Y, tspan, p)

    # Add diagnostics
    if use_land_diagnostics
        output_writer =
            CD.Writers.NetCDFWriter(domain.space.subsurface, output_dir; start_date)
        scheduled_diagnostics = CL.default_diagnostics(
            model,
            start_date,
            output_writer = output_writer,
            reduction_period = :monthly,
        )

        diagnostic_handler =
            CD.DiagnosticsHandler(scheduled_diagnostics, Y, p, tspan[1]; dt = dt)
        diag_cb = CD.DiagnosticsCallback(diagnostic_handler)
    else
        output_writer = nothing
        diag_cb = nothing
    end

    integrator = SciMLBase.init(
        prob,
        ode_algo;
        dt = dt,
        saveat = saveat,
        adaptive = false,
        callback = SciMLBase.CallbackSet(diag_cb),
    )

    # TODO move these to ClimaLand.jl as part of CoupledRadiativeFluxes
    radiative_fluxes =
        (; SW_d = CC.Fields.zeros(surface_space), LW_d = CC.Fields.zeros(surface_space))

    return BucketSimulation(
        model,
        integrator,
        area_fraction,
        output_writer,
        radiative_fluxes,
    )
end

# extensions required by Interfacer
Interfacer.get_field(sim::BucketSimulation, ::Val{:area_fraction}) = sim.area_fraction
Interfacer.get_field(sim::BucketSimulation, ::Val{:emissivity}) =
    CL.surface_emissivity(sim.model, sim.integrator.u, sim.integrator.p)
Interfacer.get_field(sim::BucketSimulation, ::Val{:roughness_buoyancy}) =
    sim.model.parameters.z_0b
Interfacer.get_field(sim::BucketSimulation, ::Val{:roughness_momentum}) =
    sim.model.parameters.z_0m
Interfacer.get_field(sim::BucketSimulation, ::Val{:surface_direct_albedo}) =
    CL.surface_albedo(sim.model, sim.integrator.u, sim.integrator.p)
Interfacer.get_field(sim::BucketSimulation, ::Val{:surface_diffuse_albedo}) =
    CL.surface_albedo(sim.model, sim.integrator.u, sim.integrator.p)
Interfacer.get_field(sim::BucketSimulation, ::Val{:surface_temperature}) =
    CL.component_temperature(sim.model, sim.integrator.u, sim.integrator.p)

"""
    Interfacer.get_field(sim::BucketSimulation, ::Val{:energy})

Extension of Interfacer.get_field that provides the total energy contained in the bucket,
computed from the temperature of the bucket and also including the latent heat of fusion
of frozen water in the snow.

This method is required by the ConservationChecker to check energy conservation.
"""
Interfacer.get_field(sim::BucketSimulation, ::Val{:energy}) =
    sim.integrator.p.bucket.total_energy

"""
    Interfacer.get_field(sim::BucketSimulation, ::Val{:water})

Extension of Interfacer.get_field that provides the total water contained in the bucket.
The total water contained in the bucket is the sum of the subsurface water storage `W`,
the snow water equivalent `σS`, and surface water content `Ws`.

This method is required by the ConservationChecker to check water conservation.
"""
function Interfacer.get_field(sim::BucketSimulation, ::Val{:water})
    ρ_cloud_liq = CL.LP.ρ_cloud_liq(sim.model.parameters.earth_param_set)
    return sim.integrator.p.bucket.total_water .* ρ_cloud_liq # kg water / m2
end

function Interfacer.update_field!(
    sim::BucketSimulation,
    ::Val{:liquid_precipitation},
    field,
)
    ρ_liq = LP.ρ_cloud_liq(sim.model.parameters.earth_param_set)
    Interfacer.remap!(sim.integrator.p.drivers.P_liq, field ./ ρ_liq)
end
function Interfacer.update_field!(sim::BucketSimulation, ::Val{:SW_d}, field)
    Interfacer.remap!(sim.radiative_fluxes.SW_d, field)
end
function Interfacer.update_field!(sim::BucketSimulation, ::Val{:LW_d}, field)
    Interfacer.remap!(sim.radiative_fluxes.LW_d, field)
end
function Interfacer.update_field!(sim::BucketSimulation, ::Val{:air_temperature}, field)
    Interfacer.remap!(sim.integrator.p.drivers.T, field)
end
function Interfacer.update_field!(sim::BucketSimulation, ::Val{:air_pressure}, field)
    Interfacer.remap!(sim.integrator.p.drivers.P, field)
end
function Interfacer.update_field!(sim::BucketSimulation, ::Val{:air_humidity}, field)
    Interfacer.remap!(sim.integrator.p.drivers.q, field)
end
function Interfacer.update_field!(sim::BucketSimulation, ::Val{:air_velocity}, u_int, v_int)
    Interfacer.remap!(sim.integrator.p.bucket.scratch1, u_int)
    Interfacer.remap!(sim.integrator.p.bucket.scratch2, v_int)
    sim.integrator.p.drivers.u .=
        StaticArrays.SVector.(
            sim.integrator.p.bucket.scratch1,
            sim.integrator.p.bucket.scratch2,
        )
end
function Interfacer.update_field!(
    sim::BucketSimulation,
    ::Val{:air_thermo_state},
    ρ_atmos,
    T_atmos,
    q_atmos,
)
    thermo_params = sim.model.parameters.earth_param_set.thermo_params
    Interfacer.remap!(sim.integrator.p.bucket.scratch1, ρ_atmos)
    Interfacer.remap!(sim.integrator.p.bucket.scratch2, T_atmos)
    Interfacer.remap!(sim.integrator.p.bucket.scratch3, q_atmos)
    sim.integrator.p.drivers.thermal_state .=
        TD.PhaseEquil_ρTq.(
            thermo_params,
            sim.integrator.p.bucket.scratch1,
            sim.integrator.p.bucket.scratch2,
            sim.integrator.p.bucket.scratch3,
        )
end
function Interfacer.update_field!(
    sim::BucketSimulation,
    ::Val{:turbulent_energy_flux},
    fields,
)
    Interfacer.remap!(sim.integrator.p.bucket.turbulent_fluxes.lhf, fields.F_lh)
    Interfacer.remap!(sim.integrator.p.bucket.turbulent_fluxes.shf, fields.F_sh)
end
function Interfacer.update_field!(sim::BucketSimulation, ::Val{:snow_precipitation}, field)
    ρ_liq = LP.ρ_cloud_liq(sim.model.parameters.earth_param_set)
    Interfacer.remap!(sim.integrator.p.drivers.P_snow, field ./ ρ_liq)
end
function Interfacer.update_field!(
    sim::BucketSimulation,
    ::Val{:turbulent_moisture_flux},
    field,
)
    ρ_liq = LP.ρ_cloud_liq(sim.model.parameters.earth_param_set)
    Interfacer.remap!(sim.integrator.p.bucket.turbulent_fluxes.vapor_flux, field ./ ρ_liq) # TODO: account for sublimation
end

Interfacer.step!(sim::BucketSimulation, t) =
    Interfacer.step!(sim.integrator, t - sim.integrator.t, true)
Interfacer.close_output_writers(sim::BucketSimulation) =
    isnothing(sim.output_writer) || close(sim.output_writer)

"""
   FieldExchanger.import_atmos_fields!(csf, sim::BucketSimulation, atmos_sim)

Import non-default coupler fields from the atmosphere simulation into the coupler fields.
These include the diffuse fraction of light, shortwave and longwave downwelling radiation,
air pressure, and CO2 concentration.

The default coupler fields will be imported by the default method implemented in
FieldExchanger.jl.
"""
function FieldExchanger.import_atmos_fields!(csf, sim::BucketSimulation, atmos_sim)
    Interfacer.get_field!(csf.P_atmos, atmos_sim, Val(:air_pressure))
    return nothing
end

"""
Extend Interfacer.add_coupler_fields! to add the fields required for BucketSimulation.

The fields added are:
- `:ρ_sfc`
- `:P_atmos`
"""
function Interfacer.add_coupler_fields!(coupler_field_names, ::BucketSimulation)
    bucket_coupler_fields = [:ρ_sfc, :P_atmos]
    push!(coupler_field_names, bucket_coupler_fields...)
end

## Extend functions for land-specific flux calculation
"""
    compute_surface_fluxes!(csf, sim::BucketSimulation, atmos_sim, thermo_params)

This function computes surface fluxes between the bucket simulation and the atmosphere.

Update the input coupler surface fields `csf` in-place with the computed fluxes
for this model. These are then summed using area-weighting across all surface
models to get the total fluxes. Fluxes where the area fraction is zero are set to zero.

Currently, this calculation is done on the land surface space, and the computed fluxes
are remapped onto the coupler boundary space as the coupler fields are updated. In the future,
we may compute fluxes in the bucket model's internal `step!` function.

# Arguments
- `csf`: [CC.Fields.Field] containing a NamedTuple of turbulent flux fields: `F_turb_ρτxz`, `F_turb_ρτyz`, `F_lh`, `F_sh`, `F_turb_moisture`.
- `sim`: [BucketSimulation] the bucket simulation to compute fluxes for.
- `atmos_sim`: [Interfacer.AbstractAtmosSimulation] the atmosphere simulation to compute fluxes with.
- `thermo_params`: [ClimaParams.ThermodynamicParameters] the thermodynamic parameters for the simulation.
"""
function FluxCalculator.compute_surface_fluxes!(
    csf,
    sim::BucketSimulation,
    atmos_sim::Interfacer.AbstractAtmosSimulation,
    thermo_params,
)
    boundary_space = axes(csf)
    FT = CC.Spaces.undertype(boundary_space)
    Y, p, t, model = sim.integrator.u, sim.integrator.p, sim.integrator.t, sim.model

    # We should change this to be on the boundary_space
    coupled_atmos = sim.model.atmos
    bucket_dest = p.bucket.turbulent_fluxes
    CL.turbulent_fluxes!(bucket_dest, coupled_atmos, model, Y, p, t)

    # Get area fraction of the land model (min = 0, max = 1)
    area_fraction = Interfacer.get_field(sim, Val(:area_fraction))

    # Combine turbulent energy fluxes from each component of the land model
    # Use temporary variables to avoid allocating
    Interfacer.remap!(csf.scalar_temp1, bucket_dest.lhf)
    Interfacer.remap!(csf.scalar_temp2, bucket_dest.shf)

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
    ρ_liq = (LP.ρ_cloud_liq(sim.model.parameters.earth_param_set))
    Interfacer.remap!(csf.scalar_temp1, bucket_dest.vapor_flux .* ρ_liq)
    @. csf.scalar_temp1 =
        ifelse(area_fraction == 0, zero(csf.scalar_temp1), csf.scalar_temp1)
    @. csf.F_turb_moisture += csf.scalar_temp1 * area_fraction

    # Combine turbulent momentum fluxes from each component of the land model
    # Note that we exclude the canopy component here for now, since we can have nonzero momentum fluxes
    #  where there is zero LAI. This should be fixed in ClimaLand.
    Interfacer.remap!(csf.scalar_temp1, bucket_dest.ρτxz)
    @. csf.scalar_temp1 =
        ifelse(area_fraction == 0, zero(csf.scalar_temp1), csf.scalar_temp1)
    @. csf.F_turb_ρτxz += csf.scalar_temp1 * area_fraction

    Interfacer.remap!(csf.scalar_temp1, bucket_dest.ρτyz)
    @. csf.scalar_temp1 =
        ifelse(area_fraction == 0, zero(csf.scalar_temp1), csf.scalar_temp1)
    @. csf.F_turb_ρτyz += csf.scalar_temp1 * area_fraction

    # Combine the buoyancy flux from each component of the land model
    # Note that we exclude the canopy component here for now, since ClimaLand doesn't
    #  include its extra resistance term in the buoyancy flux calculation.
    Interfacer.remap!(csf.scalar_temp1, bucket_dest.buoyancy_flux)
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
    surface_params = LP.surface_fluxes_parameters(sim.model.parameters.earth_param_set)
    @. csf.scalar_temp1 =
        -csf.ustar^3 / SFP.von_karman_const(surface_params) / SF.non_zero(csf.buoyancy_flux)
    @. csf.scalar_temp1 =
        ifelse(area_fraction == 0, zero(csf.scalar_temp1), csf.scalar_temp1)
    # When L_MO is infinite, avoid multiplication by zero to prevent NaN
    @. csf.L_MO +=
        ifelse(isinf(csf.scalar_temp1), csf.scalar_temp1, csf.scalar_temp1 * area_fraction)

    return nothing
end

"""
    update_sim!(sim::BucketSimulation, csf)

Updates the surface component model cache with the current coupler fields besides turbulent fluxes.

# Arguments
- `sim`: [Interfacer.AbstractSurfaceSimulation] containing a surface model simulation object.
- `csf`: [NamedTuple] containing coupler fields.
"""
function FieldExchanger.update_sim!(sim::BucketSimulation, csf)
    # radiative fluxes
    # TODO add SW_d, LW_d fields to BucketSimulation and update there instead
    Interfacer.update_field!(sim, Val(:SW_d), csf.SW_d)
    Interfacer.update_field!(sim, Val(:LW_d), csf.LW_d)

    model = sim.model
    Y = sim.integrator.u
    p = sim.integrator.p
    t = sim.integrator.t

    # Note: here we add negative signs to account for a difference in sign convention
    #  between ClimaLand.jl and ClimaCoupler.jl in SW_d and LW_d.
    σ = model.parameters.earth_param_set.Stefan
    sim.integrator.p.bucket.R_n .=
        .-(1 .- CL.surface_albedo(model, Y, p)) .* sim.radiative_fluxes.SW_d .-
        Interfacer.get_field(sim, Val(:emissivity)) .*
        (sim.radiative_fluxes.LW_d .- σ .* CL.component_temperature(model, Y, p) .^ 4)

    # precipitation
    Interfacer.update_field!(sim, Val(:liquid_precipitation), csf.P_liq)
    Interfacer.update_field!(sim, Val(:snow_precipitation), csf.P_snow)

    # Update fields for turbulent flux calculations
    # For the BucketSimulation, the coupler uses CL.turbulent_fluxes to calculate surface fluxes,
    # so we need to update the following fields in p.drivers in the bucket cache.
    Interfacer.update_field!(sim, Val(:air_velocity), csf.u_int, csf.v_int)
    Interfacer.update_field!(sim, Val(:air_temperature), csf.T_atmos)
    Interfacer.update_field!(sim, Val(:air_pressure), csf.P_atmos)
    Interfacer.update_field!(sim, Val(:air_humidity), csf.q_atmos)
    Interfacer.update_field!(
        sim,
        Val(:air_thermo_state),
        csf.ρ_atmos,
        csf.T_atmos,
        csf.q_atmos,
    )
    return nothing
end

"""
    Interfacer.set_cache!(sim::BucketSimulation, csf)

Set cache variables that cannot be initialized before the initial exchange.
This must be called after radiation, so that `p.drivers`
is filled with the initial radiation fluxes, and these can be propagated
to the rest of the cache (e.g. in canopy radative transfer).

This function does not set all the cache variables, because many are computed
as part of the tendendencies.
"""
function Interfacer.set_cache!(sim::BucketSimulation, csf)
    bucket_set_initial_cache! = CL.make_set_initial_cache(sim.model)
    bucket_set_initial_cache!(sim.integrator.p, sim.integrator.u, sim.integrator.t)

    return nothing
end

"""
    Checkpointer.get_model_prog_state(sim::BucketSimulation)

Extension of Checkpointer.get_model_prog_state to get the model state.
"""
function Checkpointer.get_model_prog_state(sim::BucketSimulation)
    return sim.integrator.u.bucket
end

function Checkpointer.get_model_cache(sim::BucketSimulation)
    return sim.integrator.p
end

function Checkpointer.restore_cache!(sim::BucketSimulation, new_cache)
    old_cache = Checkpointer.get_model_cache(sim)
    comms_ctx = ClimaComms.context(sim.model)
    Checkpointer.restore!(
        old_cache,
        new_cache,
        comms_ctx,
        ignore = Set([:rc, :params, :dss_buffer_2d, :dss_buffer_3d, :graph_context]),
    )
end

# Additional BucketSimulation getter methods for plotting debug fields
Interfacer.get_field(sim::BucketSimulation, ::Val{:σS}) = sim.integrator.u.bucket.σS
Interfacer.get_field(sim::BucketSimulation, ::Val{:Ws}) = sim.integrator.u.bucket.Ws
Interfacer.get_field(sim::BucketSimulation, ::Val{:W}) = sim.integrator.u.bucket.W

Plotting.debug_plot_fields(sim::BucketSimulation) =
    (:area_fraction, :surface_temperature, :σS, :Ws, :W)
