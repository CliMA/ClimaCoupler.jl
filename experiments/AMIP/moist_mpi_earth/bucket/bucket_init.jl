# slab_rhs!
using ClimaCore
using ClimaLSM
import ClimaLSM
include(joinpath(pkgdir(ClimaLSM), "parameters", "create_parameters.jl"))
using ClimaLSM.Bucket: BucketModel, BucketModelParameters, AbstractAtmosphericDrivers, AbstractRadiativeDrivers

import ClimaLSM.Bucket: surface_fluxes, surface_air_density, liquid_precipitation, BulkAlbedo, surface_albedo, snow_precipitation

using ClimaLSM: make_ode_function, initialize_prognostic, initialize_auxiliary

import ClimaLSM: initialize
import ClimaLSM.Parameters as LSMP
include(joinpath(pkgdir(ClimaLSM), "parameters", "create_parameters.jl"))
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


"""
    initialize(model::BucketModel, space::ClimaCore.Spaces.AbstractSpace)

A method for initializing land model variables from a predefined space.
"""
function initialize(model::BucketModel, space::ClimaCore.Spaces.AbstractSpace)
    coords = ClimaCore.Fields.coordinate_field(space)
    Y = initialize_prognostic(model, coords)
    p = initialize_auxiliary(model, coords)
    return Y, p, coords
end


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
                    radiation::PR, p,
                    ) where {FT <: AbstractFloat, P <: BucketModelParameters{FT},  PA <: CoupledAtmosphere{FT}, PR <: CoupledRadiativeFluxes{FT}}

Computes the surface flux terms at the ground for a coupled simulation:
net radiation,  SHF,  LHF,
as well as the water vapor flux (in units of m^3/m^2/s of water).
Positive fluxes indicate flow from the ground to the atmosphere.
Currently, we only support soil covered surfaces.
"""
function ClimaLSM.Bucket.surface_fluxes(
    Y::ClimaCore.Fields.FieldVector,
    p::ClimaCore.Fields.FieldVector,
    t::FT,
    parameters::BucketModelParameters{FT, P},
    atmos::CoupledAtmosphere{FT},
    radiation::CoupledRadiativeFluxes{FT},
) where {FT <: AbstractFloat, P}
    # coupler has done its thing behind the scenes already
    return (R_n = p.bucket.R_n, LHF = p.bucket.LHF, SHF = p.bucket.SHF, E = p.bucket.E)
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
    return p.bucket.P_snow
end

"""
    get_land_temp(slab_sim::BucketSimulation)

Returns the surface temperature of the earth; 
a method for the bucket model 
when used as the land model.
"""
function get_land_temp(slab_sim::BucketSimulation)
    return slab_sim.integrator.u.bucket.T_sfc
end

"""
    get_land_roughness(slab_sim::BucketSimulation)

Returns the roughness length parameters of the bucket; 
a method for the bucket model 
when used as the land model.
"""
function get_land_roughness(slab_sim::BucketSimulation)
    return slab_sim.params.z_0m, slab_sim.params.z_0b
end

"""
   land_albedo(slab_sim::BucketSimulation)

Returns the surface albedo of the earth; 
a method for the bucket model 
when used as the land model.
"""
function land_albedo(slab_sim::BucketSimulation)
    coords = ClimaCore.Fields.coordinate_field(axes(slab_sim.integrator.u.bucket.σS))
    α_land = surface_albedo.(Ref(slab_sim.params.albedo), coords, slab_sim.integrator.u.bucket.σS, slab_sim.params.σS_c)
    return parent(α_land)
end

"""
    get_land_q(slab_sim::Bucketimulation, _...)

Returns the surface specific humidity of the earth; 
a method for the bucket 
when used as the land model.
"""
function get_land_q(slab_sim::BucketSimulation, _...)
    return slab_sim.integrator.p.bucket.q_sfc
end

"""
    get_land_energy(bucket_sim::BucketSimulation, boundary_space)

Returns the internal energy per unit area of the bucket land model.
"""
function get_land_energy(bucket_sim::BucketSimulation, boundary_space)
    σS = swap_space!(bucket_sim.integrator.u.bucket.σS, boundary_space)
    T_sfc = swap_space!(bucket_sim.integrator.u.bucket.T_sfc, boundary_space)
    return bucket_sim.params.ρc_soil .* T_sfc .* bucket_sim.params.d_soil .-
        LSMP.LH_f0(bucket_sim.params.earth_param_set) .* σS
end
get_land_energy(bucket_sim::SlabSimulation, boundary_space) = get_slab_energy(bucket_sim, boundary_space)



"""
    bucket_init

Initializes the bucket model variables.
"""
function bucket_init(::Type{FT}, tspan::Tuple{FT, FT}; space, dt::FT, saveat::FT, stepper = Euler()) where {FT}

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
    T0 = FT(280.0)
    z_0m = FT(1e-2)
    z_0b = FT(1e-3)
    κ_soil = FT(0.0)# setting this to zero allows us to test energy conservation; zero flux in soil column at bottom
    ρc_soil = FT(2e6)
    τc = FT(5)*dt
    params = BucketModelParameters(d_soil, T0, κ_soil, ρc_soil, albedo, σS_c, W_f, z_0m, z_0b, τc, earth_param_set)

    args = (params, CoupledAtmosphere{FT}(), CoupledRadiativeFluxes{FT}(), nothing)
    model = BucketModel{FT, typeof.(args)...}(args...)

    # Initial conditions with no moisture
    Y, p, coords = initialize(model, space)
    anomaly = true
    hs_sfc = false
    Y.bucket.T_sfc = map(coords) do coord
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
    # add ρ_sfc to cache
    # this needs to be initialized!!! Turbulent surface fluxes need this set to be computed.
    ρ_sfc = zeros(space) .+ FT(1.1)
    P_liq = zeros(space) .+ FT(0.0)
    P_snow = zeros(space) .+ FT(0.0)
    variable_names = (propertynames(p.bucket)..., :ρ_sfc, :P_liq, :P_snow)
    orig_fields = map(x -> getproperty(p.bucket, x), propertynames(p.bucket))
    fields = (orig_fields..., ρ_sfc, P_liq, P_snow)
    p_new = ClimaCore.Fields.FieldVector(; :bucket => (; zip(variable_names, fields)...))

    ode_function! = make_ode_function(model)

    prob = ODEProblem(ode_function!, Y, tspan, p_new)
    integrator = init(prob, stepper; dt = dt, saveat = saveat)

    BucketSimulation(params, Y, space, integrator)
end
