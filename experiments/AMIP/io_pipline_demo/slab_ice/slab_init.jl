# slab_ice

# sea-ice parameters
struct IceSlabParameters{FT <: AbstractFloat}
    h::FT # sea ice height
    ρ::FT
    c::FT
    T_base::FT
    z0m::FT
    z0b::FT
    T_freeze::FT # temperature at freezing point [K] 
    k_ice::FT   # thermal conductivity of ice [W / m / K]
    α::FT # albedo
end

# init simulation
function slab_ice_space_init(::Type{FT}, space, p) where {FT}
    Y = Fields.FieldVector(T_sfc = ones(space) .* p.T_freeze)
    return Y
end


function ice_rhs!(du, u, p, t)
    dY = du
    Y = u
    FT = eltype(dY)

    params = p.Ya.params
    F_aero = p.Ya.F_aero
    F_rad = p.Ya.F_rad
    ice_mask = p.Ya.ice_mask

    F_conductive = @. params.k_ice / (params.h) * (params.T_base - Y.T_sfc)
    rhs = @. (-F_aero - F_rad + F_conductive) / (params.h * params.ρ * params.c)
    parent(dY.T_sfc) .= apply_mask.(parent(ice_mask), >, parent(rhs), FT(0))
end

function slab_ice_init(::Type{FT}; tspan, saveat, dt, space, ice_mask, stepper = Euler()) where {FT}

    params = IceSlabParameters(
        FT(2),
        FT(1500.0),
        FT(800.0),
        FT(280.0),
        FT(1e-3),
        FT(1e-5),
        FT(273.15),
        FT(0.0),# k_ice
        FT(0.38), # albedo
    )

    Y = slab_ice_space_init(FT, space, params)
    Ya = (;
        params = params,
        F_aero = ClimaCore.Fields.zeros(space),
        F_rad = ClimaCore.Fields.zeros(space),
        ice_mask = ice_mask,
    )

    problem = OrdinaryDiffEq.ODEProblem(ice_rhs!, Y, tspan, (; Ya = Ya))
    integrator = OrdinaryDiffEq.init(problem, stepper, dt = dt, saveat = saveat)


    SlabSimulation(params, Y, space, integrator)
end

get_ice_mask(h_ice, FT) = h_ice > FT(0) ? FT(1) : FT(0)
