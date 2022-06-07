# slab_ice

# sea-ice parameters
struct IceSlabParameters# <: CLIMAParameters.AbstractEarthParameterSet{F} 
    h::FT
    ρ::FT
    c::FT
    T_init::FT
    z0m::FT
    z0b::FT
    ρc_ml::FT    # density times heat transfer coefficient for mixed layer [J / m2 / K ]
    F0_base::FT  # ice base transfer coefficient [W / m2 / K]
    T_base::FT   # ice base temperature [K] 
    L_ice::FT    # latent heat coefficient for ice [J / m3]
    h_ml::FT     # mixed layer depth [m] 
    T_freeze::FT # temperature at freezing point [K] 
    k_ice::FT   # thermal conductivity of ice [W / m / K]
    α::FT # albedo
    σ::FT # Stefan Boltzmann constant
end

# init simulation
function slab_ice_space_init(::Type{FT}, space, p) where {FT}

    # prognostic variable
    Y = Fields.FieldVector(T_sfc = ones(space) .* p.T_freeze, h_ice = zeros(space), T_ml = ones(space) .* p.T_freeze)

    return Y, space
end

"""
    solve_ice!(dT_sfc, T_sfc, (parameters, F_accumulated), t)

slab RHS with an implicit solve ice and explicit (forward Euler) solve for ocean

"""
function solve_ice!(integ, Δt)

    Y = integ.u
    Ya = integ.p.Ya
    p = integ.p.p

    FT = eltype(Y)
    ocean_qflux = FT(0)

    # prognostic
    T_ml = Y.T_ml
    T_sfc = Y.T_sfc
    h_ice = Y.h_ice

    # auxiliary
    F_atm = Ya.F_aero + Ya.F_rad

    ∂F_atm∂T_sfc = get_∂F_atm∂T_sfc(p, T_sfc, Ya) # this will be passed from atmos/SF.jl

    ΔT_ml = similar(F_atm)
    Δh_ice = similar(F_atm)
    F_conductive = similar(F_atm)
    ΔT_sfc = similar(F_atm)

    # ice thickness and mixed layer temperature changes due to atmosphereic and ocean fluxes 
    if h_ice[1] > 0 # ice-covered 
        @. ΔT_ml = -(p.F0_base * (T_ml - p.T_base) + ocean_qflux) * Δt / (p.h_ml * p.ρc_ml)
        @. Δh_ice = (F_atm - p.F0_base * (T_ml - p.T_base)) * Δt / p.L_ice
    else # ice-free
        @. ΔT_ml = -(F_atm + ocean_qflux) * Δt / (p.h_ml * p.ρc_ml)
        @. Δh_ice = 0
    end

    # adjust if transition to ice-covered 
    if (T_ml[1] + ΔT_ml[1] < p.T_freeze)
        @. Δh_ice = Δh_ice - (T_ml + ΔT_ml - p.T_freeze) * (p.h_ml * p.ρc_ml) / p.L_ice
        @. ΔT_ml = p.T_freeze - T_ml
    end

    # adjust if transition to ice-free
    if ((h_ice[1] > 0) & (h_ice[1] + Δh_ice[1] <= 0))
        @. ΔT_ml = ΔT_ml - (h_ice + Δh_ice) * p.L_ice / (p.h_ml * p.ρc_ml)
        @. Δh_ice = -h_ice
    end

    # solve for T_sfc
    if (h_ice[1] + Δh_ice[1] > 0) #  surface is ice-covered
        # if ice covered, solve implicity (for now one Newton iteration: ΔT_s = - F(T_s) / dF(T_s)/dT_s )
        @. F_conductive = p.k_ice / (h_ice + Δh_ice) * (p.T_base - T_sfc)
        @. ΔT_sfc = (-F_atm + F_conductive) / (p.k_ice / (h_ice + Δh_ice) + ∂F_atm∂T_sfc)
        if (T_sfc[1] + ΔT_sfc[1] > p.T_freeze)
            @. ΔT_sfc = p.T_freeze - T_sfc
        end
        # surface is ice-covered, so update T_sfc as ice surface temperature
        @. Y.T_sfc += ΔT_sfc
    else # ice-free, so update T_sfc as mixed layer temperature
        @. Y.T_sfc = T_ml + ΔT_ml
    end

    # update state
    @. Y.T_ml += ΔT_ml
    @. Y.h_ice += Δh_ice

    @. Ya.ice_mask = get_ice_mask(h_ice)

    return nothing
end

get_∂F_atm∂T_sfc(p, T_sfc, Ya) = FT(4) * (FT(1) - p.α) * p.σ * T_sfc^3 + Ya.∂F_aero∂T_sfc

function ∑tendencies_ice_stub(du, u, p, t)
    dY = du
    Y = u
    FT = eltype(dY)

    if p.prescribed_sic !== nothing # Prognostic eqn for T_sfc

        @show "pres SIC"
        p, F_aero, F_rad, ice_mask = Ya

        F_conductive = p.k_ice / (p.h_ice) * (p.T_base - Y.T_sfc)
        rhs = @. (F_aero + F_rad + F_conductive) * p.h_ice / p.k_ice
        parent(dY.T_sfc) .= apply_mask.(parent(ice_mask), >, parent(rhs), FT(0), value = FT(0))
    else
        solve_ice!((; u = u, p = p), p.Δt) # timestepping outside of DeffEq (but DeffEq still used here for saving vars in `integ.sol`)

        @. dY.T_ml = FT(0)
        @. dY.h_ice = FT(0)
        @. dY.T_sfc = FT(0)
    end
end

function slab_ice_init(
    ::Type{FT},
    tspan;
    stepper = Euler(),
    nelements = 6,
    npolynomial = 4,
    dt = 0.02,
    saveat = 1.0e10,
    space = nothing,
    mask = nothing,
    prescribed_sic = nothing,
    ocean_params = (; ρ = FT(1e3), c = FT(4e3), h = FT(1) )
) where {FT}

    params = IceSlabParameters(
        FT(2),
        FT(1500.0),
        FT(800.0),
        FT(280.0),
        FT(1e-3),
        FT(1e-5),
        FT(ocean_params.ρ * ocean_params.c), #rho c
        FT(120),
        FT(273.16),
        FT(3e8),
        FT(ocean_params.h), # h_ml
        FT(273.16),
        FT(2),
        FT(0.38),
        FT(5.67e-8),
    ) # TODO: better interface, use CLIMAParameters

    ice_mask = prescribed_sic !== nothing ? get_ice_mask.(prescribed_sic .- FT(25)) : ClimaCore.Fields.zeros(space) # here 25% and lower is considered ice free # TODO: generalize to a smaoot function of ice fraction

    Y, space = slab_ice_space_init(FT, space, params)
    Ya = (
        params = params,
        F_aero = ClimaCore.Fields.zeros(space),
        ∂F_aero∂T_sfc = ClimaCore.Fields.zeros(space),
        F_rad = ClimaCore.Fields.zeros(space),
        mask = mask,
        ice_mask = ice_mask,
        Δt = dt,
        prescribed_sic = prescribed_sic,
    ) #auxiliary
    problem = OrdinaryDiffEq.ODEProblem(slab_ocean_rhs!, Y, tspan, Ya)
    integrator = OrdinaryDiffEq.init(problem, stepper, dt = dt, saveat = saveat)

    SlabSimulation(params, Y, space, integrator)
end

get_ice_mask(h_ice) = h_ice > FT(0) ? FT(1) : FT(0)
