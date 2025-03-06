import Dates
import SciMLBase
import Statistics
import ClimaComms
import ClimaCore as CC
import ClimaTimeSteppers as CTS
import Thermodynamics as TD
import ClimaLand as CL
import ClimaLand.Parameters as LP
import ClimaDiagnostics as CD
import ClimaCoupler: Checkpointer, FluxCalculator, Interfacer
using NCDatasets

include("../shared/restore.jl")

###
### Functions required by ClimaCoupler.jl for a SurfaceModelSimulation
###
"""
    BucketSimulation{M, D, I}

The bucket model simulation object.
"""
struct BucketSimulation{M, D, I, A} <: Interfacer.LandModelSimulation
    model::M
    domain::D
    integrator::I
    area_fraction::A
end
Interfacer.name(::BucketSimulation) = "BucketSimulation"

"""
    get_new_cache(p, Y, energy_check)

Returns a new `p` with an updated field to store e_per_area if energy conservation
    checks are turned on.
"""
function get_new_cache(p, Y, energy_check)
    if energy_check
        e_per_area_field = CC.Fields.zeros(axes(Y.bucket.W))
        return merge(p, (; e_per_area = e_per_area_field))
    else
        return p
    end
end

"""
    bucket_init

Initializes the bucket model variables.
"""
function BucketSimulation(
    ::Type{FT},
    tspan::Tuple{TT, TT},
    config::String,
    albedo_type::String,
    land_initial_condition::String,
    land_temperature_anomaly::String,
    output_dir::String;
    space,
    dt::TT,
    saveat::Vector{TT},
    area_fraction,
    stepper = CTS.RK4(),
    date_ref::Dates.DateTime,
    t_start::TT,
    energy_check::Bool,
    surface_elevation,
    use_land_diagnostics::Bool,
) where {FT, TT <: Union{Float64, ITime}}
    if config != "sphere"
        println(
            "Currently only spherical shell domains are supported; single column set-up will be addressed in future PR.",
        )
        @assert config == "sphere"
    end

    α_snow = FT(0.8) # snow albedo
    if albedo_type == "map_static" # Read in albedo from static data file (default type)
        # By default, this uses a file containing bareground albedo without a time component. Snow albedo is specified separately.
        albedo = CL.Bucket.PrescribedBaregroundAlbedo{FT}(α_snow, space)
    elseif albedo_type == "map_temporal" # Read in albedo from data file containing data over time
        # By default, this uses a file containing linearly-interpolated monthly data of clear-sky albedo, generated from CERES.
        if pkgversion(CL) < v"0.15"
            albedo = CL.Bucket.PrescribedSurfaceAlbedo{FT}(date_ref, t_start, space)
        else
            albedo = CL.Bucket.PrescribedSurfaceAlbedo{FT}(
                date_ref,
                space;
                albedo_file_path = CL.Artifacts.ceres_albedo_dataset_path(),
                varname = "sw_alb_clr",
            )
        end
    elseif albedo_type == "function" # Use prescribed function of lat/lon for surface albedo
        function α_bareground(coordinate_point)
            (; lat, long) = coordinate_point
            return typeof(lat)(0.38)
        end
        albedo = CL.Bucket.PrescribedBaregroundAlbedo{FT}(α_snow, α_bareground, space)
    else
        error("invalid albedo type $albedo_type")
    end

    d_soil = FT(3.5) # soil depth
    z_0m = FT(1e-3) # roughness length for momentum over smooth bare soil
    z_0b = FT(1e-3) # roughness length for tracers over smooth bare soil
    τc = FT(float(dt)) # This is the timescale on which snow exponentially damps to zero, in the case where all
    # the snow would melt in time `τc`. It prevents us from having to specially time step in cases where
    # all the snow melts in a single timestep.
    σS_c = FT(0.2) # critical snow water equivalent
    W_f = FT(0.2) # bucket capacity
    κ_soil = FT(1.5) # soil conductivity
    ρc_soil = FT(2e6) # soil volumetric heat capacity

    params = CL.Bucket.BucketModelParameters(FT; albedo, z_0m, z_0b, τc, σS_c, W_f, κ_soil, ρc_soil)

    n_vertical_elements = 7
    # Note that this does not take into account topography of the surface, which is OK for this land model.
    # But it must be taken into account when computing surface fluxes, for Δz.
    domain = make_land_domain(space, (-d_soil, FT(0.0)), n_vertical_elements)
    args = (params, CL.CoupledAtmosphere{FT}(), CL.CoupledRadiativeFluxes{FT}(), domain)
    model = CL.Bucket.BucketModel{FT, typeof.(args)...}(args...)

    # Initial conditions with no moisture
    Y, p, coords = CL.initialize(model)
    p = get_new_cache(p, Y, energy_check)

    # Get temperature anomaly function
    T_functions = Dict("aquaplanet" => temp_anomaly_aquaplanet, "amip" => temp_anomaly_amip)
    haskey(T_functions, land_temperature_anomaly) ||
        error("land temp anomaly function $land_temperature_anomaly not supported")
    temp_anomaly = T_functions[land_temperature_anomaly]

    # Set temperature IC including anomaly, based on atmospheric setup
    # Bucket surface temperature is in `p.bucket.T_sfc` (ClimaLand.jl)
    lapse_rate = FT(6.5e-3)
    T_sfc_0 = FT(271)
    @. Y.bucket.T = T_sfc_0 + temp_anomaly(coords.subsurface)
    # `surface_elevation` is a ClimaCore.Fields.Field(`half` level)
    orog_adjusted_T = CC.Fields.field_values(Y.bucket.T) .- lapse_rate .* CC.Fields.field_values(surface_elevation)
    # Adjust T based on surface elevation (p.bucket.T_sfc is then set using the
    # set_initial_cache! function)
    parent(Y.bucket.T) .= parent(orog_adjusted_T)

    Y.bucket.W .= 0.15
    Y.bucket.Ws .= 0.0
    Y.bucket.σS .= 0.0

    # Overwrite initial conditions with interpolated values from a netcdf file using
    # the `SpaceVaryingInputs` tool. We expect the file to contain the following variables:
    # - `W`, for subsurface water storage (2D),
    # - `Ws`, for surface water content (2D),
    # - `T`, for soil temperature (3D),
    # - `S`, for snow water equivalent (2D).

    if !isempty(land_initial_condition)
        ds = NCDataset(land_initial_condition)
        has_all_variables = all(key -> haskey(ds, key), ["W", "Ws", "T", "S"])
        @assert has_all_variables "The land iniital condition file is expected to contain the variables W, Ws, T, and S (read documentation about requirements)."
        close(ds)

        surface_space = domain.space.surface
        subsurface_space = domain.space.subsurface
        regridder_type = :InterpolationsRegridder
        extrapolation_bc = (Interpolations.Periodic(), Interpolations.Flat(), Interpolations.Flat())
        Y.bucket.W .= SpaceVaryingInput(
            land_initial_condition,
            "W",
            surface_space;
            regridder_type,
            regridder_kwargs = (; extrapolation_bc,),
        )
        Y.bucket.Ws .= SpaceVaryingInput(
            land_initial_condition,
            "Ws",
            surface_space;
            regridder_type,
            regridder_kwargs = (; extrapolation_bc,),
        )
        Y.bucket.T .= SpaceVaryingInput(
            land_initial_condition,
            "T",
            subsurface_space;
            regridder_type,
            regridder_kwargs = (; extrapolation_bc,),
        )
        Y.bucket.σS .= SpaceVaryingInput(
            land_initial_condition,
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
    bucket_ode_function = CTS.ClimaODEFunction(T_exp! = exp_tendency!, dss! = CL.dss!)
    prob = SciMLBase.ODEProblem(bucket_ode_function, Y, tspan, p)

    # Add diagnostics
    if use_land_diagnostics
        netcdf_writer = CD.Writers.NetCDFWriter(domain.space.subsurface, output_dir)
        scheduled_diagnostics =
            CL.default_diagnostics(model, date_ref, output_writer = netcdf_writer, average_period = :monthly)

        diagnostic_handler = CD.DiagnosticsHandler(scheduled_diagnostics, Y, p, t_start; dt = dt)
        diag_cb = CD.DiagnosticsCallback(diagnostic_handler)
    else
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

    sim = BucketSimulation(model, (; domain = domain, soil_depth = d_soil), integrator, area_fraction)

    # DSS state to ensure we have continuous fields
    dss_state!(sim)
    return sim
end

# extensions required by Interfacer
Interfacer.get_field(sim::BucketSimulation, ::Val{:area_fraction}) = sim.area_fraction
Interfacer.get_field(sim::BucketSimulation, ::Val{:beta}) =
    CL.surface_evaporative_scaling(sim.model, sim.integrator.u, sim.integrator.p)
Interfacer.get_field(sim::BucketSimulation, ::Val{:roughness_buoyancy}) = sim.model.parameters.z_0b
Interfacer.get_field(sim::BucketSimulation, ::Val{:roughness_momentum}) = sim.model.parameters.z_0m
Interfacer.get_field(sim::BucketSimulation, ::Val{:surface_direct_albedo}) =
    CL.surface_albedo(sim.model, sim.integrator.u, sim.integrator.p)
Interfacer.get_field(sim::BucketSimulation, ::Val{:surface_diffuse_albedo}) =
    CL.surface_albedo(sim.model, sim.integrator.u, sim.integrator.p)
Interfacer.get_field(sim::BucketSimulation, ::Val{:surface_humidity}) =
    CL.surface_specific_humidity(sim.model, sim.integrator.u, sim.integrator.p, sim.integrator.t)
Interfacer.get_field(sim::BucketSimulation, ::Val{:surface_temperature}) =
    CL.surface_temperature(sim.model, sim.integrator.u, sim.integrator.p, sim.integrator.t)

"""
    Interfacer.get_field(bucket_sim::BucketSimulation, ::Val{:energy})

Extension of Interfacer.get_field that provides the total energy contained in the bucket, including the latent heat due to snow melt.
"""
function Interfacer.get_field(bucket_sim::BucketSimulation, ::Val{:energy})
    # required by ConservationChecker
    e_per_area = bucket_sim.integrator.p.e_per_area .= 0
    CC.Operators.column_integral_definite!(
        e_per_area,
        bucket_sim.model.parameters.ρc_soil .* bucket_sim.integrator.u.bucket.T,
    )

    e_per_area .+=
        -LP.LH_f0(bucket_sim.model.parameters.earth_param_set) .*
        LP.ρ_cloud_liq(bucket_sim.model.parameters.earth_param_set) .* bucket_sim.integrator.u.bucket.σS

    return e_per_area
end

"""
    Interfacer.get_field(bucket_sim::BucketSimulation, ::Val{:water})

Extension of Interfacer.get_field that provides the total water contained in the bucket, including the liquid water in snow.
"""
function Interfacer.get_field(bucket_sim::BucketSimulation, ::Val{:water})
    ρ_cloud_liq = CL.LP.ρ_cloud_liq(bucket_sim.model.parameters.earth_param_set)
    return
    @. (bucket_sim.integrator.u.bucket.σS + bucket_sim.integrator.u.bucket.W + bucket_sim.integrator.u.bucket.Ws) *
       ρ_cloud_liq  # kg water / m2
end

function Interfacer.update_field!(sim::BucketSimulation, ::Val{:air_density}, field)
    parent(sim.integrator.p.bucket.ρ_sfc) .= parent(field)
end
function Interfacer.update_field!(sim::BucketSimulation, ::Val{:liquid_precipitation}, field)
    ρ_liq = (LP.ρ_cloud_liq(sim.model.parameters.earth_param_set))
    parent(sim.integrator.p.drivers.P_liq) .= parent(field ./ ρ_liq)
end
function Interfacer.update_field!(sim::BucketSimulation, ::Val{:radiative_energy_flux_sfc}, field)
    parent(sim.integrator.p.bucket.R_n) .= parent(field)
end
function Interfacer.update_field!(sim::BucketSimulation, ::Val{:turbulent_energy_flux}, field)
    parent(sim.integrator.p.bucket.turbulent_fluxes.shf) .= parent(field)
end
function Interfacer.update_field!(sim::BucketSimulation, ::Val{:snow_precipitation}, field)
    ρ_liq = (LP.ρ_cloud_liq(sim.model.parameters.earth_param_set))
    parent(sim.integrator.p.drivers.P_snow) .= parent(field ./ ρ_liq)
end
function Interfacer.update_field!(sim::BucketSimulation, ::Val{:turbulent_moisture_flux}, field)
    ρ_liq = (LP.ρ_cloud_liq(sim.model.parameters.earth_param_set))
    parent(sim.integrator.p.bucket.turbulent_fluxes.vapor_flux) .= parent(field ./ ρ_liq) # TODO: account for sublimation
end

# extensions required by FieldExchanger
Interfacer.step!(sim::BucketSimulation, t) = Interfacer.step!(sim.integrator, t - sim.integrator.t, true)
Interfacer.reinit!(sim::BucketSimulation) = Interfacer.reinit!(sim.integrator)

"""
Extend Interfacer.add_coupler_fields! to add the fields required for BucketSimulation.

The fields added are:
- `:ρ_sfc`
- `:F_radiative` (for radiation input)
- `:P_liq` (for precipitation input)
- `:P_snow` (for precipitation input)
"""
function Interfacer.add_coupler_fields!(coupler_field_names, ::BucketSimulation)
    bucket_coupler_fields = [:ρ_sfc, :F_radiative, :P_liq, :P_snow]
    push!(coupler_field_names, bucket_coupler_fields...)
end

# extensions required by FluxCalculator (partitioned fluxes)
function FluxCalculator.update_turbulent_fluxes!(sim::BucketSimulation, fields::NamedTuple)
    (; F_turb_energy, F_turb_moisture) = fields
    turbulent_fluxes = sim.integrator.p.bucket.turbulent_fluxes
    turbulent_fluxes.shf .= F_turb_energy
    earth_params = sim.model.parameters.earth_param_set
    turbulent_fluxes.vapor_flux .= F_turb_moisture ./ LP.ρ_cloud_liq(earth_params)
    return nothing
end

# extension of FluxCalculator.FluxCalculator.surface_thermo_state, overriding the saturated-surface default
function FluxCalculator.surface_thermo_state(
    sim::BucketSimulation,
    thermo_params::TD.Parameters.ThermodynamicsParameters,
    thermo_state_int,
)

    T_sfc = Interfacer.get_field(sim, Val(:surface_temperature))
    # Note that the surface air density, ρ_sfc, is computed using the atmospheric state at the first level and making ideal gas
    # and hydrostatic balance assumptions. The land model does not compute the surface air density so this is
    # a reasonable stand-in.
    ρ_sfc = FluxCalculator.extrapolate_ρ_to_sfc.(thermo_params, thermo_state_int, T_sfc) # ideally the # calculate elsewhere, here just getter...
    q_sfc = Interfacer.get_field(sim, Val(:surface_humidity)) # already calculated in rhs! (cache)
    @. TD.PhaseEquil_ρTq.(thermo_params, ρ_sfc, T_sfc, q_sfc)
end

"""
    Checkpointer.get_model_prog_state(sim::BucketSimulation)

Extension of Checkpointer.get_model_prog_state to get the model state.
"""
function Checkpointer.get_model_prog_state(sim::BucketSimulation)
    return sim.integrator.u.bucket
end

"""
    temp_anomaly_aquaplanet(coord)

Introduce a temperature IC anomaly for the aquaplanet case.
The values for this case follow the moist Held-Suarez setup of Thatcher &
Jablonowski (2016, eq. 6), consistent with ClimaAtmos aquaplanet.
"""
temp_anomaly_aquaplanet(coord) = 29 * exp(-coord.lat^2 / (2 * 26^2))

"""
    temp_anomaly_amip(coord)

Introduce a temperature IC anomaly for the AMIP case.
The values used in this case have been tuned to align with observed temperature
and result in stable simulations.
"""
temp_anomaly_amip(coord) = 40 * cosd(coord.lat)^4

"""
    make_land_domain(
        atmos_boundary_space::CC.Spaces.SpectralElementSpace2D,
        zlim::Tuple{FT, FT},
        nelements_vert::Int,) where {FT}

Creates the land model domain from the horizontal space of the atmosphere, and information
about the number of elements and extent of the vertical domain.
"""
function make_land_domain(
    atmos_boundary_space::CC.Spaces.SpectralElementSpace2D,
    zlim::Tuple{FT, FT},
    nelements_vert::Int,
) where {FT}
    @assert zlim[1] < zlim[2]
    depth = zlim[2] - zlim[1]

    mesh = CC.Spaces.topology(atmos_boundary_space).mesh

    radius = mesh.domain.radius
    nelements_horz = mesh.ne
    npolynomial = CC.Spaces.Quadratures.polynomial_degree(CC.Spaces.quadrature_style(atmos_boundary_space))
    nelements = (nelements_horz, nelements_vert)
    vertdomain = CC.Domains.IntervalDomain(
        CC.Geometry.ZPoint(FT(zlim[1])),
        CC.Geometry.ZPoint(FT(zlim[2]));
        boundary_names = (:bottom, :top),
    )

    vertmesh = CC.Meshes.IntervalMesh(vertdomain, CC.Meshes.Uniform(), nelems = nelements[2])
    verttopology = CC.Topologies.IntervalTopology(vertmesh)
    vert_center_space = CC.Spaces.CenterFiniteDifferenceSpace(verttopology)
    subsurface_space = CC.Spaces.ExtrudedFiniteDifferenceSpace(atmos_boundary_space, vert_center_space)
    space = (; surface = atmos_boundary_space, subsurface = subsurface_space)

    fields = CL.Domains.get_additional_domain_fields(subsurface_space)

    return CL.Domains.SphericalShell{FT}(radius, depth, nothing, nelements, npolynomial, space, fields)
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

"""
    dss_state!(sim::BucketSimulation)

Perform DSS on the state of a component simulation, intended to be used
before the initial step of a run. This method acts on bucket land simulations.
The `dss!` function of ClimaLand must be called because it uses either the 2D
or 3D dss buffer stored in the cache depending on space of each variable in
`sim.integrator.u`.
"""
function dss_state!(sim::BucketSimulation)
    CL.dss!(sim.integrator.u, sim.integrator.p, sim.integrator.t)
end
