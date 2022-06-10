# slab_ice

# sea-ice parameters
struct IceSlabParameters# <: CLIMAParameters.AbstractEarthParameterSet{F} 
    h::FT
    ρ::FT
    c::FT
    α::FT # albedo
    T_freeze::FT
end

# init simulation
function slab_ice_space_init(::Type{FT}, space, p) where {FT}

    # prognostic variable
    Y = Fields.FieldVector(T_sfc = ones(space) .* p.T_freeze,)

    return Y
end
function slab_ice_rhs!(dY, Y, Ya, t)
    (; params, F_aero, F_rad) = Ya
    dY.T_sfc .= @. (F_aero + F_rad) / (params.h * params.ρ * params.c)
end

function slab_ice_init(
    ::Type{FT},
    tspan;
    stepper = Euler(),
    dt,
    saveat,
    space,
    ice_mask,
) where {FT}

    params = IceSlabParameters(
        FT(2),
        FT(1500.0),
        FT(800.0),
        FT(0.6),
        FT(273.16),
    ) # TODO: better interface, use CLIMAParameters



    Y = slab_ice_space_init(FT, space, params)
    Ya = (
        params = params,
        F_aero = ClimaCore.Fields.zeros(space),
        F_rad = ClimaCore.Fields.zeros(space),
        ice_mask = ice_mask,
    )
    problem = OrdinaryDiffEq.ODEProblem(slab_ice_rhs!, Y, tspan, Ya)
    integrator = OrdinaryDiffEq.init(problem, stepper, dt = dt, saveat = saveat)

    SlabSimulation(params, Y, space, integrator)
end

get_ice_mask(h_ice) = h_ice > FT(0) ? FT(1) : FT(0)
