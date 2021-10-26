########
# Set up parameters
########
parameters = (
    ρₒ = 1.0,  # reference density
    cₛ = 1e-2, # sound speed
    ℓᵐ = 10,   # jet thickness, (larger is thinner)
    ℓ = 20,   # perturbation thickness, (larger is thinner)
    m = 2,    # number of sign changes on equator for perturbation
    ϕᵖ = π / 2 * 0.05, # of centerdness of perturbation
    ϵ = 0.3,  # perturbation amplitude
    vˢ = 5e-4, # velocity scale
    α = 2e-4,
    Ω = 1e-3,
    λ_coupler = (500 / 60 / 86400), #(L_airsea / τ_airsea)
)


########
# Set up inital condition
########
uᵐ(p, λ, ϕ, r) = p.ℓᵐ * sech(p.ℓᵐ * ϕ)^2
vᵐ(p, λ, ϕ, r) = 0.0
hᵐ(p, λ, ϕ, r) = 0.0

u1(p, λ, ϕ, r) = p.ℓ * 2 * (ϕ - p.ϕᵖ) * exp(-p.ℓ * (ϕ - p.ϕᵖ)^2) * cos(ϕ) * cos(2 * (ϕ - p.ϕᵖ)) * sin(p.m * λ)
u2(p, λ, ϕ, r) = exp(-p.ℓ * (ϕ - p.ϕᵖ)^2) * sin(ϕ) * cos(2 * (ϕ - p.ϕᵖ)) * sin(p.m * λ)
u3(p, λ, ϕ, r) = 2 * exp(-p.ℓ * (ϕ - p.ϕᵖ)^2) * cos(ϕ) * sin(2 * (ϕ - p.ϕᵖ)) * sin(p.m * λ)
uᵖ(p, λ, ϕ, r) = u1(p, λ, ϕ, r) + u2(p, λ, ϕ, r) + u3(p, λ, ϕ, r)
vᵖ(p, λ, ϕ, r) = p.m * exp(-p.ℓ * (ϕ - p.ϕᵖ)^2) * cos(2 * (ϕ - p.ϕᵖ)) * cos(p.m * λ)
hᵖ(p, λ, ϕ, r) = 0.0

ρ₀(p, λ, ϕ, r) = p.ρₒ
ρuʳᵃᵈ(p, λ, ϕ, r) = 0
ρuˡᵃᵗ(p, λ, ϕ, r) = p.vˢ * ρ₀(p, λ, ϕ, r) * (p.ϵ * vᵖ(p, λ, ϕ, r))
ρuˡᵒⁿ(p, λ, ϕ, r) = p.vˢ * ρ₀(p, λ, ϕ, r) * (uᵐ(p, λ, ϕ, r) + p.ϵ * uᵖ(p, λ, ϕ, r))
ρθ₀B(p, λ, ϕ, r) = ρ₀(p, λ, ϕ, r) * tanh(p.ℓᵐ * ϕ)
ρθ₀A(p, λ, ϕ, r) = 0

# Cartesian Representation for initialization (boiler plate really)
ρ₀ᶜᵃʳᵗ(p, x...) = ρ₀(p, lon(x...), lat(x...), rad(x...))
ρu⃗₀ᶜᵃʳᵗ(p, x...) = (
    ρuʳᵃᵈ(p, lon(x...), lat(x...), rad(x...)) * r̂(x...) +
    ρuˡᵃᵗ(p, lon(x...), lat(x...), rad(x...)) * ϕ̂(x...) +
    ρuˡᵒⁿ(p, lon(x...), lat(x...), rad(x...)) * λ̂(x...)
)
ρθ₀ᶜᵃʳᵗA(p, x...) = ρθ₀A(p, lon(x...), lat(x...), rad(x...))
ρθ₀ᶜᵃʳᵗB(p, x...) = ρθ₀B(p, lon(x...), lat(x...), rad(x...))
