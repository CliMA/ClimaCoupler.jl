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
import ClimaCoupler: Checkpointer, FluxCalculator, Interfacer, FieldExchanger
using NCDatasets
include("climaland_helpers.jl")

include("../shared/restore.jl")

###
### Functions required by ClimaCoupler.jl for a SurfaceModelSimulation
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
} <: Interfacer.LandModelSimulation
    model::M
    integrator::I
    area_fraction::A
    output_writer::OW
    radiative_fluxes::RF
end

"""
    bucket_init

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
    land_temperature_anomaly::String = "amip",
    use_land_diagnostics::Bool = true,
    stepper = CTS.RK4(),
    albedo_type::String = "map_static",
    bucket_initial_condition::String = "",
    coupled_param_dict = CP.create_toml_dict(FT),
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

    args = (params, CL.CoupledAtmosphere{FT}(), CL.CoupledRadiativeFluxes{FT}(), domain)
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
            reduction_period = :hourly,
            reduction_type = :instantaneous,
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
Interfacer.get_field(sim::BucketSimulation, ::Val{:beta}) =
    CL.surface_evaporative_scaling(sim.model, sim.integrator.u, sim.integrator.p)
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
    CL.surface_temperature(sim.model, sim.integrator.u, sim.integrator.p, sim.integrator.t)

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

function Interfacer.update_field!(sim::BucketSimulation, ::Val{:air_density}, field)
    Interfacer.remap!(sim.integrator.p.bucket.ρ_sfc, field)
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
Extend Interfacer.add_coupler_fields! to add the fields required for BucketSimulation.

The fields added are:
- `:ρ_sfc`
"""
function Interfacer.add_coupler_fields!(coupler_field_names, ::BucketSimulation)
    bucket_coupler_fields = [:ρ_sfc]
    push!(coupler_field_names, bucket_coupler_fields...)
end

# extensions required by FluxCalculator
function FluxCalculator.update_turbulent_fluxes!(sim::BucketSimulation, fields::NamedTuple)
    (; F_lh, F_sh, F_turb_moisture) = fields
    Interfacer.update_field!(sim, Val(:turbulent_energy_flux), (; F_lh, F_sh))
    Interfacer.update_field!(sim, Val(:turbulent_moisture_flux), F_turb_moisture)
    return nothing
end

"""
    update_sim!(sim::BucketSimulation, csf)

Updates the surface component model cache with the current coupler fields besides turbulent fluxes.

# Arguments
- `sim`: [Interfacer.SurfaceModelSimulation] containing a surface model simulation object.
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
        (sim.radiative_fluxes.LW_d .- σ .* CL.surface_temperature(model, Y, p, t) .^ 4)

    # precipitation
    Interfacer.update_field!(sim, Val(:liquid_precipitation), csf.P_liq)
    Interfacer.update_field!(sim, Val(:snow_precipitation), csf.P_snow)
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
    restore!(
        old_cache,
        new_cache,
        comms_ctx,
        ignore = Set([:rc, :params, :dss_buffer_2d, :dss_buffer_3d, :graph_context]),
    )
end
