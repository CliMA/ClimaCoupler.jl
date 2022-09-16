# slab_rhs!
using ClimaCore
using ClimaLSM
import ClimaLSM
include(joinpath(pkgdir(ClimaLSM), "parameters", "create_parameters.jl"))
using ClimaLSM.Bucket: BucketModel, BucketModelParameters, AbstractAtmosphericDrivers, AbstractRadiativeDrivers

import ClimaLSM.Bucket:
    surface_fluxes,
    net_radiation,
    surface_air_density,
    liquid_precipitation,
    BulkAlbedo,
    surface_albedo,
    snow_precipitation

using ClimaLSM: make_ode_function, initialize, obtain_surface_space

"""
    BucketSimulation{P, Y, D, I}

The bucket model simulation object.
"""
struct BucketSimulation{P, Y, D, I}
    params::P
    Y_init::Y
    domain::D
    integrator::I
end

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
    surface_fluxes(Y, p, 
                    t::FT,
                    parameters::P,
                    atmos::PA,
                    ) where {FT <: AbstractFloat, P <: BucketModelParameters{FT},  PA <: CoupledAtmosphere{FT}}

Computes the turbulent surface fluxes terms at the ground for a coupled simulation.
"""
function ClimaLSM.Bucket.surface_fluxes(
    Y::ClimaCore.Fields.FieldVector,
    p::ClimaCore.Fields.FieldVector,
    t::FT,
    parameters::BucketModelParameters{FT, P},
    atmos::CoupledAtmosphere{FT},
) where {FT <: AbstractFloat, P}
    # coupler has done its thing behind the scenes already
    return (turbulent_energy_flux = p.bucket.turbulent_energy_flux, evaporation = p.bucket.evaporation)
end



"""
    net_radiation(Y, p, 
                    t::FT,
                    parameters::P,
                    radiation::PR, p,
                    ) where {FT <: AbstractFloat, P <: BucketModelParameters{FT}, PR <: CoupledRadiativeFluxes{FT}}

Computes the net radiative flux at the ground for a coupled simulation.
"""
function ClimaLSM.Bucket.net_radiation(
    Y::ClimaCore.Fields.FieldVector,
    p::ClimaCore.Fields.FieldVector,
    t::FT,
    parameters::BucketModelParameters{FT, P},
    radiation::CoupledRadiativeFluxes{FT},
) where {FT <: AbstractFloat, P}
    # coupler has done its thing behind the scenes already
    return p.bucket.R_n
end


"""
    ClimaLSM.Bucket.surface_air_density(p, atmos::CoupledAtmosphere)

an extension of the bucket model method which returns the surface air 
density in the case of a coupled simulation.
"""
function ClimaLSM.Bucket.surface_air_density(p, atmos::CoupledAtmosphere)
    # coupler has filled this in
    return p.bucket.ρ_sfc
end

"""
    ClimaLSM.Bucket.liquid_precipitation(p, atmos::CoupledAtmosphere, t)

an extension of the bucket model method which returns the precipitation
(m/s) in the case of a coupled simulation.
"""
function liquid_precipitation(p, atmos::CoupledAtmosphere, t)
    # coupler has filled this in
    return p.bucket.P_liq
end

"""
    ClimaLSM.Bucket.snow_precipitation(p, atmos::CoupledAtmosphere, t)

an extension of the bucket model method which returns the precipitation
(m/s) in the case of a coupled simulation.
"""
function snow_precipitation(p, atmos::CoupledAtmosphere, t)
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
    config::String;
    space,
    dt::FT,
    saveat::FT,
    stepper = Euler(),
) where {FT}
    if config != "sphere"
        println(
            "Currently only spherical shell domains are supported; single column set-up will be addressed in future PR.",
        )
        @assert config == "sphere"
    end

    earth_param_set = create_lsm_parameters(FT)
    function α_soil(coordinate_point)
        (; lat, long) = coordinate_point
        return typeof(lat)(0.2)
    end

    α_snow = FT(0.8) # snow albedo
    albedo = BulkAlbedo{FT}(α_snow, α_soil)
    σS_c = FT(0.2)
    W_f = FT(0.15)
    d_soil = FT(3.5) # soil depth
    z_0m = FT(1e-2)
    z_0b = FT(1e-3)
    κ_soil = FT(0.7)
    ρc_soil = FT(2e6)
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
    Y, p, coords = initialize(model)
    anomaly = true
    hs_sfc = false
    Y.bucket.T = map(coords.subsurface) do coord
        T_sfc_0 = FT(280.0)
        radlat = coord.lat / FT(180) * pi
        ΔT = FT(0)
        if anomaly == true
            anom_ampl = FT(0)# this is zero, no anomaly
            lat_0 = FT(60) / FT(180) * pi
            lon_0 = FT(-90) / FT(180) * pi
            radlon = coord.long / FT(180) * pi
            stdev = FT(5) / FT(180) * pi
            ΔT = anom_ampl * exp(-((radlat - lat_0)^2 / 2stdev^2 + (radlon - lon_0)^2 / 2stdev^2))
        elseif hs_sfc == true
            ΔT = -FT(60) * sin(radlat)^2
        end
        T_sfc_0 + ΔT
    end

    Y.bucket.W .= 0.14
    Y.bucket.Ws .= 0.0
    Y.bucket.σS .= 0.0
    ρ_sfc = zeros(axes(Y.bucket.W)) .+ FT(1.1)
    P_liq = zeros(axes(Y.bucket.W)) .+ FT(0.0)
    P_snow = zeros(axes(Y.bucket.W)) .+ FT(0.0)
    variable_names = (propertynames(p.bucket)..., :ρ_sfc, :P_liq, :P_snow)
    orig_fields = map(x -> getproperty(p.bucket, x), propertynames(p.bucket))
    fields = (orig_fields..., ρ_sfc, P_liq, P_snow)
    p_new = ClimaCore.Fields.FieldVector(; :bucket => (; zip(variable_names, fields)...))

    ode_function! = make_ode_function(model)

    prob = ODEProblem(ode_function!, Y, tspan, p_new)
    integrator = init(prob, stepper; dt = dt, saveat = saveat)

    BucketSimulation(params, Y, (; domain = domain, soil_depth = d_soil), integrator)
end
