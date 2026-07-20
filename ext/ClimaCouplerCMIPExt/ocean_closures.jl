# Ported from ClimaOcean v0.10.0, src/OceanConfigurations/{OceanConfigurations,one_degree_tripolar}.jl,
# and NumericalEarth v0.5.8, src/Oceans/ocean_simulation.jl

import Oceananigans.Units: days
import Oceananigans.TurbulenceClosures.TKEBasedVerticalDiffusivities:
    CATKEVerticalDiffusivity, CATKEMixingLength, CATKEEquation

@inline νhb(i, j, k, grid, ℓx, ℓy, ℓz, clock, fields, λ) =
    OC.Operators.Az(i, j, k, grid, ℓx, ℓy, ℓz)^2 / λ

# Background tracer diffusivity following Henyey et al. (1986), as implemented by Harrison and Hallberg (2008)
# and used in GFDL OM4 (Adcroft et al. 2019).
@inline henyey_diffusivity(x, y, z, t) = max(2e-6, 3e-5 * abs(sind(y)))

@inline ν_step_simple(x, y, z, t) = ifelse(z >= -100, 1e-2, 1e-4)
@inline κ_step_simple(x, y, z, t) = z >= -10 ? 5e-2 : z >= -100 ? 1e-2 : 1e-5

function default_ocean_closure(FT = OC.defaults.FloatType)
    mixing_length = CATKEMixingLength(Cᵇ = 0.01)
    turbulent_kinetic_energy_equation = CATKEEquation(Cᵂϵ = 1.0)
    return CATKEVerticalDiffusivity(
        OC.VerticallyImplicitTimeDiscretization(),
        FT;
        mixing_length,
        turbulent_kinetic_energy_equation,
    )
end

"""
    simplified_ocean_closure(FT=Oceananigans.defaults.FloatType)

A minimal closure suitable for testing on memory-limited GPUs (e.g. P100). Uses `ConvectiveAdjustmentVerticalDiffusivity`with a background
vertical diffusivity and viscosity, avoiding the large parameter space of CATKE + Gent-McWilliams + biharmonic closures.
"""
function simplified_ocean_closure(FT = OC.defaults.FloatType)
    horizontal_viscosity = OC.HorizontalScalarBiharmonicDiffusivity(
        ν = νhb,
        discrete_form = true,
        parameters = 10days,
    )
    convective_mixing = OC.ConvectiveAdjustmentVerticalDiffusivity(
        FT;
        convective_κz = 1.0,
        convective_νz = 1.0,
    )
    vertical_mixing = OC.VerticalScalarDiffusivity(
        OC.VerticallyImplicitTimeDiscretization(),
        FT;
        ν = ν_step_simple,
        κ = κ_step_simple,
    )
    return (horizontal_viscosity, convective_mixing, vertical_mixing)
end

function default_one_degree_closure(;
    κ_skew = 500,
    κ_symmetric = 200,
    biharmonic_timescale = 15days,
    background_κ = henyey_diffusivity,
    background_ν = 1e-5,
)
    catke = default_ocean_closure()
    eddy = OC.IsopycnalSkewSymmetricDiffusivity(; κ_skew, κ_symmetric)
    horizontal_viscosity = OC.HorizontalScalarBiharmonicDiffusivity(
        ν = νhb,
        discrete_form = true,
        parameters = biharmonic_timescale,
    )
    vertical_diffusivity = OC.VerticalScalarDiffusivity(ν = background_ν, κ = background_κ)
    return (catke, eddy, horizontal_viscosity, vertical_diffusivity)
end
