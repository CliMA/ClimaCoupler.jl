import ClimateMachine.NumericalFluxes: numerical_flux_first_order!

roe_average(ρ⁻, ρ⁺, var⁻, var⁺) =
    (sqrt(ρ⁻) * var⁻ + sqrt(ρ⁺) * var⁺) / (sqrt(ρ⁻) + sqrt(ρ⁺))

function numerical_flux_first_order!(
    ::RoeNumericalFlux,
    model::CNSE3D,
    fluxᵀn::Vars{S},
    n⁻::SVector,
    state⁻::Vars{S},
    aux⁻::Vars{A},
    state⁺::Vars{S},
    aux⁺::Vars{A},
    t,
    direction,
) where {S, A}
    numerical_flux_first_order!(
        CentralNumericalFluxFirstOrder(),
        model,
        fluxᵀn,
        n⁻,
        state⁻,
        aux⁻,
        state⁺,
        aux⁺,
        t,
        direction,
    )

    FT = eltype(fluxᵀn)

    # - states
    ρ⁻ = state⁻.ρ
    ρu⁻ = state⁻.ρu
    ρθ⁻ = state⁻.ρθ

    # constructed states
    u⁻ = ρu⁻ / ρ⁻
    θ⁻ = ρθ⁻ / ρ⁻
    uₙ⁻ = u⁻' * n⁻

    # in general thermodynamics
    p⁻ = calc_pressure(model.eos, state⁻)
    c⁻ = calc_sound_speed(model.eos, state⁻)

    # + states
    ρ⁺ = state⁺.ρ
    ρu⁺ = state⁺.ρu
    ρθ⁺ = state⁺.ρθ

    # constructed states
    u⁺ = ρu⁺ / ρ⁺
    θ⁺ = ρθ⁺ / ρ⁺
    uₙ⁺ = u⁺' * n⁻

    # in general thermodynamics
    p⁺ = calc_pressure(model.eos, state⁺)
    c⁺ = calc_sound_speed(model.eos, state⁺)

    # construct roe averges
    ρ = sqrt(ρ⁻ * ρ⁺)
    u = roe_average(ρ⁻, ρ⁺, u⁻, u⁺)
    θ = roe_average(ρ⁻, ρ⁺, θ⁻, θ⁺)
    c = roe_average(ρ⁻, ρ⁺, c⁻, c⁺)

    # construct normal velocity
    uₙ = u' * n⁻

    # differences
    Δρ = ρ⁺ - ρ⁻
    Δp = p⁺ - p⁻
    Δu = u⁺ - u⁻
    Δρθ = ρθ⁺ - ρθ⁻
    Δuₙ = Δu' * n⁻

    # constructed values
    c⁻² = 1 / c^2
    w1 = abs(uₙ - c) * (Δp - ρ * c * Δuₙ) * 0.5 * c⁻²
    w2 = abs(uₙ + c) * (Δp + ρ * c * Δuₙ) * 0.5 * c⁻²
    w3 = abs(uₙ) * (Δρ - Δp * c⁻²)
    w4 = abs(uₙ) * ρ
    w5 = abs(uₙ) * (Δρθ - θ * Δp * c⁻²)

    # fluxes!!!
    
    fluxᵀn.ρ -= (w1 + w2 + w3) * 0.5
    fluxᵀn.ρu -=
        (
            w1 * (u - c * n⁻) +
            w2 * (u + c * n⁻) +
            w3 * u +
            w4 * (Δu - Δuₙ * n⁻)
        ) * 0.5
    fluxᵀn.ρθ -= ((w1 + w2) * θ + w5) * 0.5
    
    #=
    max_wavespeed = 10.0 # max(c⁻, c⁺)
    fluxᵀn.ρ -= max_wavespeed * 0.5 * Δρ * 1.0
    fluxᵀn.ρu -=
        (
            max_wavespeed * Δu
        ) * 0.5
    fluxᵀn.ρθ -= max_wavespeed * Δρθ * 1.0
    =#

    return nothing
end
