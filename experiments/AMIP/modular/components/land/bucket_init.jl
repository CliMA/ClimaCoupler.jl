# slab_rhs!
using ClimaCore
using ClimaLSM
import ClimaLSM
import ClimaTimeSteppers as CTS
import Thermodynamics as TD
using Dates: DateTime

include(joinpath(pkgdir(ClimaLSM), "parameters", "create_parameters.jl"))
using ClimaLSM.Bucket: BucketModel, BucketModelParameters, AbstractAtmosphericDrivers, AbstractRadiativeDrivers
using ClimaComms: AbstractCommsContext

import ClimaLSM.Bucket:
    surface_fluxes,
    net_radiation,
    surface_air_density,
    liquid_precipitation,
    BulkAlbedoTemporal,
    BulkAlbedoStatic,
    BulkAlbedoFunction,
    surface_albedo,
    snow_precipitation

using ClimaLSM:
    make_exp_tendency, initialize, obtain_surface_space, make_set_initial_aux_state, surface_evaporative_scaling

import ClimaCoupler.Interfacer: LandModelSimulation, get_field, update_field!, name
import ClimaCoupler.FieldExchanger: step!, reinit!
import ClimaCoupler.FluxCalculator: update_turbulent_fluxes_point!, surface_thermo_state

"""
    BucketSimulation{M, Y, D, I}

The bucket model simulation object.
"""
struct BucketSimulation{M, Y, D, I, A} <: LandModelSimulation
    model::M
    Y_init::Y
    domain::D
    integrator::I
    area_fraction::A
end
name(::BucketSimulation) = "BucketSimulation"

include("./bucket_utils.jl")

"""
    CoupledRadiativeFluxes{FT} <: AbstractRadiativeDrivers{FT}

To be used when coupling to an atmosphere model; internally, used for
multiple dispatch on `surface_fluxes`.
"""
struct CoupledRadiativeFluxes{FT} <: AbstractRadiativeDrivers{FT} end

"""
    CoupledAtmosphere{FT} <: AbstractAtmosphericDrivers{FT}

To be used when coupling to an atmosphere model; internally, used for
multiple dispatch on `surface_fluxes`.
"""
struct CoupledAtmosphere{FT} <: AbstractAtmosphericDrivers{FT} end

"""
    surface_fluxes(atmos::CoupledAtmosphere{FT},
                    model::BucketModel{FT},
                    Y,
                    p,
                    ) where {FT <: AbstractFloat}

Computes the turbulent surface fluxes terms at the ground for a coupled simulation.
Note that `Ch` is not used with the current implementation of the bucket model,
but will be used once the canopy is incorporated.

The turbulent energy flux is currently not split up between latent and sensible
heat fluxes. This will be fixed once `lhf` and `shf` are added to the bucket's
cache.
"""
function ClimaLSM.Bucket.surface_fluxes(
    atmos::CoupledAtmosphere{FT},
    model::BucketModel{FT},
    Y,
    p,
    _...,
) where {FT <: AbstractFloat}
    space = model.domain.space.surface
    return (
        lhf = ClimaCore.Fields.zeros(space),
        shf = p.bucket.turbulent_energy_flux,
        vapor_flux = p.bucket.evaporation,
        Ch = ClimaCore.Fields.similar(p.bucket.evaporation),
    )
end

"""
    net_radiation(radiation::CoupledRadiativeFluxes{FT},
                    model::BucketModel{FT},
                    Y,
                    p,
                    _...,
                    ) where {FT <: AbstractFloat}

Computes the net radiative flux at the ground for a coupled simulation.
"""
function ClimaLSM.Bucket.net_radiation(
    radiation::CoupledRadiativeFluxes{FT},
    model::BucketModel{FT},
    Y::ClimaCore.Fields.FieldVector,
    p::NamedTuple,
    _...,
) where {FT <: AbstractFloat}
    # coupler has done its thing behind the scenes already
    return p.bucket.R_n
end


"""
    ClimaLSM.surface_air_density(
                    atmos::CoupledAtmosphere{FT},
                    model::BucketModel{FT},
                    Y,
                    p,
                    _...,
                )
an extension of the bucket model method which returns the surface air
density in the case of a coupled simulation.
"""
function surface_air_density(
    atmos::CoupledAtmosphere{FT},
    model::BucketModel{FT},
    Y,
    p,
    _...,
) where {FT <: AbstractFloat}
    return p.bucket.ρ_sfc
end

"""
    ClimaLSM.Bucket.liquid_precipitation(atmos::CoupledAtmosphere, p, t)

an extension of the bucket model method which returns the precipitation
(m/s) in the case of a coupled simulation.
"""
function liquid_precipitation(atmos::CoupledAtmosphere, p, t)
    # coupler has filled this in
    return p.bucket.P_liq
end

"""
    ClimaLSM.Bucket.snow_precipitation(atmos::CoupledAtmosphere, p, t)

an extension of the bucket model method which returns the precipitation
(m/s) in the case of a coupled simulation.
"""
function snow_precipitation(atmos::CoupledAtmosphere, p, t)
    # coupler has filled this in
    return p.bucket.P_snow
end

"""
    bucket_init

Initializes the bucket model variables.
"""
function bucket_init(
    ::Type{FT},
    tspan::Tuple{FT, FT},
    config::String,
    albedo_type::String,
    land_temperature_anomaly::String,
    comms_ctx::AbstractCommsContext,
    regrid_dirpath::String;
    space,
    dt::FT,
    saveat::FT,
    area_fraction,
    stepper = CTS.RK4(),
    date_ref::DateTime,
    t_start::FT,
) where {FT}
    @show "start of bucket_init"
    if config != "sphere"
        println(
            "Currently only spherical shell domains are supported; single column set-up will be addressed in future PR.",
        )
        @assert config == "sphere"
    end

    earth_param_set = create_lsm_parameters(FT)

    α_snow = FT(0.8) # snow albedo
    if albedo_type == "map_static" # Read in albedo from static data file (default type)
        # By default, this uses a file containing bareground albedo without a time component. Snow albedo is specified separately.
        if ClimaComms.iamroot(comms_ctx)
            @show "before BulkAlbedoStatic call in iamroot"
            albedo = BulkAlbedoStatic{FT}(regrid_dirpath, α_snow = α_snow)
            @show "after BulkAlbedoStatic call in iamroot"
        end
        @show "right before barrier"
        ClimaComms.barrier(comms_ctx)
        @show "right after barrier"
        albedo = BulkAlbedoStatic{FT}(regrid_dirpath, α_snow = α_snow)
        @show "after BulkAlbedoStatic all process call"
    elseif albedo_type == "map_temporal" # Read in albedo from data file containing data over time
        # By default, this uses a file containing linearly-interpolated monthly data of total albedo, generated by CESM2's land model (CLM).
        if ClimaComms.iamroot(comms_ctx)
            albedo = BulkAlbedoTemporal{FT}(regrid_dirpath, date_ref, t_start, space)
        end
        ClimaComms.barrier(comms_ctx)
        albedo = BulkAlbedoTemporal{FT}(regrid_dirpath, date_ref, t_start, space)
    elseif albedo_type == "function" # Use prescribed function of lat/lon for surface albedo
        function α_sfc(coordinate_point)
            (; lat, long) = coordinate_point
            return typeof(lat)(0.38)
        end
        albedo = BulkAlbedoFunction{FT}(α_snow, α_sfc)
    else
        error("invalid albedo type $albedo_type")
    end

    σS_c = FT(0.2)
    W_f = FT(10)
    d_soil = FT(3.5) # soil depth
    z_0m = FT(1e-2)
    z_0b = FT(1e-3)
    κ_soil = FT(0.7)
    ρc_soil = FT(2e8)
    t_crit = dt # This is the timescale on which snow exponentially damps to zero, in the case where all
    # the snow would melt in time t_crit. It prevents us from having to specially time step in cases where
    # all the snow melts in a single timestep.
    params = BucketModelParameters(κ_soil, ρc_soil, albedo, σS_c, W_f, z_0m, z_0b, t_crit, earth_param_set)
    n_vertical_elements = 7
    # Note that this does not take into account topography of the surface, which is OK for this land model.
    # But it must be taken into account when computing surface fluxes, for Δz.
    domain = make_lsm_domain(space, (-d_soil, FT(0.0)), n_vertical_elements)
    args = (params, CoupledAtmosphere{FT}(), CoupledRadiativeFluxes{FT}(), domain)
    model = BucketModel{FT, typeof.(args)...}(args...)

    # Initial conditions with no moisture
    @show "before bucket initialize"
    Y, p, coords = initialize(model)
    @show "after bucket initialize"
    Y.bucket.T = map(coords.subsurface) do coord
        T_sfc_0 = FT(271.0)
        radlat = coord.lat / FT(180) * pi
        ΔT = FT(0)
        if land_temperature_anomaly == "zonally_asymmetric"
            anom_ampl = FT(0)# this is zero, no anomaly
            lat_0 = FT(60) / FT(180) * pi
            lon_0 = FT(-90) / FT(180) * pi
            radlon = coord.long / FT(180) * pi
            stdev = FT(5) / FT(180) * pi
            ΔT = anom_ampl * exp(-((radlat - lat_0)^2 / 2stdev^2 + (radlon - lon_0)^2 / 2stdev^2))
        elseif land_temperature_anomaly == "aquaplanet"
            ΔT = FT(29) * exp(-coord.lat^2 / (2 * 26^2))
        elseif land_temperature_anomaly == "amip"
            ΔT = FT(40 * cos(radlat)^4)
        else
            ΔT = FT(0)
        end
        T_sfc_0 + ΔT
    end

    Y.bucket.W .= 10#0.14
    Y.bucket.Ws .= 0.0
    Y.bucket.σS .= 0.0
    P_liq = zeros(axes(Y.bucket.W)) .+ FT(0.0)
    P_snow = zeros(axes(Y.bucket.W)) .+ FT(0.0)
    variable_names = (propertynames(p.bucket)..., :P_liq, :P_snow)
    orig_fields = map(x -> getproperty(p.bucket, x), propertynames(p.bucket))
    fields = (orig_fields..., P_liq, P_snow)
    p_new = (;
        :bucket => (; zip(variable_names, fields)...),
        :dss_buffer_2d => p.dss_buffer_2d,
        :dss_buffer_3d => p.dss_buffer_3d,
    )

    # Set initial aux variable values
    @show "before set_initial_aux_state"
    set_initial_aux_state! = make_set_initial_aux_state(model)
    set_initial_aux_state!(p_new, Y, tspan[1])
    @show "after set_initial_aux_state"

    exp_tendency! = make_exp_tendency(model)
    ode_algo = CTS.ExplicitAlgorithm(stepper)
    bucket_ode_function = CTS.ClimaODEFunction(T_exp! = exp_tendency!, dss! = ClimaLSM.dss!)
    prob = ODEProblem(bucket_ode_function, Y, tspan, p_new)
    integrator = init(prob, ode_algo; dt = dt, saveat = saveat, adaptive = false)
    @show "after bucket integrator init"

    sim = BucketSimulation(model, Y, (; domain = domain, soil_depth = d_soil), integrator, area_fraction)

    # DSS state to ensure we have continuous fields
    dss_state!(sim)
    @show "end of bucket_init"
    return sim
end
