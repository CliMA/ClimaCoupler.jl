# slab_rhs!
using ClimaCore
using CLIMAParameters: AbstractEarthParameterSet
using ClimaLSM
using ClimaLSM.Bucket: BucketModel, BucketModelParameters, AbstractAtmosphericDrivers, AbstractRadiativeDrivers

import ClimaLSM.Bucket: surface_fluxes, surface_air_density, liquid_precipitation

using ClimaLSM: make_ode_function, initialize_prognostic, initialize_auxiliary

import ClimaLSM: initialize

struct EarthParameterSet <: AbstractEarthParameterSet end

struct BucketSimulation{P, Y, D, I}
    params::P
    Y_init::Y
    domain::D
    integrator::I
end

# init simulation
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
    parameters::BucketModelParameters{FT,P},
    atmos::CoupledAtmosphere{FT},
    radiation::CoupledRadiativeFluxes{FT}
) where {
    FT <: AbstractFloat,
    P
}
    # coupler has done its thing behind the scenes already
    return (
        R_n = p.bucket.R_n,
        LHF = p.bucket.LHF,
        SHF = p.bucket.SHF,
        E = p.bucket.E,
    )
end


function surface_air_density(p, atmos::CoupledAtmosphere)
    # coupler has filled this in
    return p.bucket.ρ_sfc
end

function liquid_precipitation(p, atmos::CoupledAtmosphere, t)
    # coupler has filled this in
    return p.bucket.P_liq
end

get_bucket_energy(bucket_sim, T_sfc) = bucket_sim.params.ρc_soil .* T_sfc .* bucket_sim.params.d_soil

function bucket_init(
    ::Type{FT},
    tspan::Tuple{FT,FT};
    space,
    dt::FT,
    saveat::FT,
    stepper = Euler(),
) where {FT}

    earth_param_set = EarthParameterSet()
    α_soil = FT(0.2) # soil albedo
    α_snow = FT(0.8) # snow albedo
    S_c = FT(0.2)
    W_f = FT(0.15)
    d_soil = FT(100.0) # soil depth
    T0 = FT(280.0)
    z_0m = FT(1e-2)
    z_0b = FT(1e-3)
    κ_soil = FT(0.0)# setting this to zero allows us to test energy conservation
    ρc_soil = FT(2e6)
    params = BucketModelParameters(d_soil,T0,κ_soil,ρc_soil,α_soil,α_snow,S_c,W_f,z_0m,z_0b,earth_param_set)

    args  = (params, CoupledAtmosphere{FT}(), CoupledRadiativeFluxes{FT}(), nothing)
    model= BucketModel{FT, typeof.(args)...}(args...)

    # Initial conditions with no moisture
    Y, p, coords = initialize(model, space)
    anomaly = false
    hs_sfc = false
    Y.bucket.T_sfc = map(coords) do coord
        T_sfc_0 = FT(280.0)
        radlat = coord.lat / FT(180) * pi
        ΔT = FT(0)
        if anomaly == true
            anom_ampl = FT(0)
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
    Y.bucket.S .= 0.0 # no snow
    # add ρ_sfc to cache
    # this needs to be initialized!!!
    ρ_sfc = zeros(space) .+ FT(1.1)
    P_liq = zeros(space) .+ FT(0.0)
    variable_names = (propertynames(p.bucket)..., :ρ_sfc, :P_liq)
    orig_fields = map(x -> getproperty(p.bucket,x), propertynames(p.bucket))
    fields = (orig_fields..., ρ_sfc, P_liq)
    p_new = ClimaCore.Fields.FieldVector(; :bucket => (;zip(variable_names, fields)...))
    ode_function! = make_ode_function(model)
    
    prob = ODEProblem(ode_function!, Y, tspan, p_new)
    integrator = init(prob, stepper; dt = dt, saveat = saveat)

    BucketSimulation(params, Y, space, integrator)
end
