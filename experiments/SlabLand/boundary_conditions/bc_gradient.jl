"""
    calc_boundary_state!(::NumericalFluxGradient, ::Impenetrable{FreeSlip}, ::ModelSetup)
apply free slip boundary condition for velocity
sets non-reflective ghost point
"""
function calc_boundary_state!(
    ::NumericalFluxGradient,
    ::Impenetrable{FreeSlip},
    ::Union{ModelSetup},
    ::Union{AbstractDiffusion, Nothing},
    state⁺,
    aux⁺,
    n⁻,
    state⁻,
    aux⁻,
    t,
    _...,
)
    state⁺.ρ = state⁻.ρ

    ρu⁻ = state⁻.ρu
    state⁺.ρu = ρu⁻ - n⁻ ⋅ ρu⁻ .* SVector(n⁻)

    return nothing
end

"""
    calc_boundary_state!(::NumericalFluxGradient, ::Impenetrable{NoSlip}, ::ModelSetup)
apply no slip boundary condition for velocity
set numerical flux to zero for U
"""
@inline function calc_boundary_state!(
    ::NumericalFluxGradient,
    ::Impenetrable{NoSlip},
    ::Union{ModelSetup},
    ::Union{AbstractDiffusion, Nothing},
    state⁺,
    aux⁺,
    n⁻,
    state⁻,
    aux⁻,
    t,
    _...,
)
    FT = eltype(state⁺)
    state⁺.ρu = @SVector zeros(FT, 3)

    return nothing
end

"""
    calc_boundary_state!(::NumericalFluxGradient, ::Insulating, ::HBModel)

apply insulating boundary condition for temperature
sets transmissive ghost point
"""
@inline function calc_boundary_state!(
    ::Union{NumericalFluxGradient},
    ::Union{Insulating, CoupledPrimaryBoundary},
    ::Union{ModelSetup},
    ::Union{AbstractDiffusion, Nothing},
    state⁺,
    aux⁺,
    n⁻,
    state⁻,
    aux⁻,
    t,
    _...,
)
    state⁺.ρθ = state⁻.ρθ

    return nothing
end

calc_boundary_state!(
    ::Union{NumericalFluxGradient},
    ::Union{Insulating, CoupledSecondaryBoundary},
    ::Union{SlabLandModelSetup},
    _...,) = nothing

