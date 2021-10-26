abstract type AbstractBoundaryCondition end

struct DefaultBC <: AbstractBoundaryCondition end

Base.@kwdef struct BulkFormulaTemperature{𝒯, 𝒰, 𝒱} <: AbstractBoundaryCondition
    drag_coef_temperature::𝒯
    drag_coef_moisture::𝒰
    surface_temperature::𝒱
end

Base.@kwdef struct CoupledPrimarySlabOceanBC{𝒯, 𝒰} <: AbstractBoundaryCondition
    drag_coef_temperature::𝒯
    drag_coef_moisture::𝒰
end

abstract type TemperatureBC end
struct Insulating <: TemperatureBC end
struct CoupledSecondaryAtmosModelBC <: TemperatureBC end

function numerical_boundary_flux_first_order!(
    numerical_flux::NumericalFluxFirstOrder,
    ::DefaultBC,
    balance_law::DryAtmosModel,
    fluxᵀn::Vars{S},
    n̂::SVector,
    state⁻::Vars{S},
    aux⁻::Vars{A},
    state⁺::Vars{S},
    aux⁺::Vars{A},
    t,
    direction,
    state1⁻::Vars{S},
    aux1⁻::Vars{A},
) where {S, A}
    state⁺.ρ = state⁻.ρ
    state⁺.ρe = state⁻.ρe
    state⁺.ρq = state⁻.ρq

    # project and reflect for impenetrable condition, but 
    # leave tangential component untouched
    ρu⁻ = state⁻.ρu
    state⁺.ρu = ρu⁻ - 2n̂ ⋅ ρu⁻ .* SVector(n̂)

    numerical_flux_first_order!(numerical_flux, balance_law, fluxᵀn, n̂, state⁻, aux⁻, state⁺, aux⁺, t, direction)
end

function numerical_boundary_flux_first_order!(
    numerical_flux::NumericalFluxFirstOrder,
    bctype::BulkFormulaTemperature,
    model::DryAtmosModel,
    fluxᵀn::Vars{S},
    n̂::SVector,
    state⁻::Vars{S},
    aux⁻::Vars{A},
    state⁺::Vars{S},
    aux⁺::Vars{A},
    t,
    direction,
    state1⁻::Vars{S},
    aux1⁻::Vars{A},
) where {S, A}
    # Impenetrable free-slip condition to reflect and project momentum 
    # at the boundary
    numerical_boundary_flux_first_order!(
        numerical_flux,
        DefaultBC(),
        model,
        fluxᵀn,
        n̂,
        state⁻,
        aux⁻,
        state⁺,
        aux⁺,
        t,
        direction,
        state1⁻,
        aux1⁻,
    )

    E, H = calc_ocean_sfc_fluxes(model.physics, bctype, state⁻, aux⁻)
    LH_v0 = model.physics.parameters.LH_v0

    fluxᵀn.ρ -= E / LH_v0
    fluxᵀn.ρe -= E + H
    fluxᵀn.ρq -= E / LH_v0
end

function numerical_boundary_flux_first_order!(
    numerical_flux::NumericalFluxFirstOrder,
    bctype::CoupledPrimarySlabOceanBC,
    model::DryAtmosModel,
    fluxᵀn::Vars{S},
    n̂::SVector,
    state⁻::Vars{S},
    aux⁻::Vars{A},
    state⁺::Vars{S},
    aux⁺::Vars{A},
    t,
    direction,
    state1⁻::Vars{S},
    aux1⁻::Vars{A},
) where {S, A}
    # Impenetrable free-slip condition to reflect and project momentum 
    # at the boundary

    numerical_boundary_flux_first_order!(
        numerical_flux,
        DefaultBC(),
        model,
        fluxᵀn,
        n̂,
        state⁻,
        aux⁻,
        state⁺,
        aux⁺,
        t,
        direction,
        state1⁻,
        aux1⁻,
    )

    # The following will be moved the the second-order kernel (as a Neumann BC) once available
    LH_v0 = model.physics.parameters.LH_v0
    E, H = calc_ocean_sfc_fluxes(model.physics, bctype, state⁻, aux⁻) #[W/m^2]

    fluxᵀn.ρ -= E / LH_v0
    fluxᵀn.ρe -= E + H
    fluxᵀn.ρq -= E / LH_v0

end
