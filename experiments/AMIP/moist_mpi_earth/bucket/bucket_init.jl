# slab_rhs!
using ClimaCore
using CLIMAParameters: AbstractEarthParameterSet
using ClimaLSM
using ClimaLSM.Bucket:
    BucketModel,
    BucketModelParameters,
    CoupledAtmosphere,
    CoupledRadiativeFluxes
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

get_slab_energy(slab_sim, T_sfc) = slab_sim.params.ρ .* slab_sim.params.c .* T_sfc .* slab_sim.params.h

function bucket_init(
    ::Type{FT},
    tspan;
    stepper = Euler(),
    nelements = 6,
    npolynomial = 4,
    dt = 0.02,
    saveat = 1.0e10,
    space = nothing,
    mask = nothing,
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
    κ_soil = FT(1.5)
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
    ρ_sfc = similar(p.bucket.E)
    names = (propertynames(p.bucket)..., :ρ_sfc)
    orig_fields = map(x -> getproperty(p.bucket,x), propertynames(p.bucket))
    fields = (orig_fields..., ρ_sfc)
    p_new = ClimaCore.Fields.FieldVector(; :bucket => (;zip(names, fields)...))
    ode_function! = make_ode_function(model)
    
    prob = ODEProblem(ode_function!, Y, tspan, p_new)
    integrator = init(prob, stepper; dt = dt, saveat = saveat)

    BucketSimulation(params, Y, space, integrator)
end
