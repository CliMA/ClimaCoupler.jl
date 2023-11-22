# from ClimaAtmos:test/test_cases/run_3d_baroclinic_wave.jl, removing `step!`

using OrdinaryDiffEq: SSPRK33
using ClimaCorePlots, Plots

import ClimaCore
import ClimaCore: Fields, Geometry, Operators
import ClimaCore.Geometry: ⊗

using ClimaAtmos.Utils.InitialConditions: init_3d_baroclinic_wave
using ClimaAtmos.Domains
using ClimaAtmos.BoundaryConditions
using ClimaAtmos.Models: ConstantViscosity
using ClimaAtmos.Models.Nonhydrostatic3DModels
using ClimaAtmos.Simulations

using CLIMAParameters
struct DryBaroclinicWaveParameters <: CLIMAParameters.AbstractEarthParameterSet end

# initiate the lower Atmos boundary Field
function atmos_bcfield_init(domain)
    center_space, face_space = make_function_space(domain)
    ClimaCore.Fields.level(ClimaCore.Fields.zeros(center_space), 1)
end

function atmos_bcvectorfield_init(domain)
    center_space, face_space = make_function_space(domain)
    zero_field = ClimaCore.Fields.zeros(center_space)
    uv_local = Geometry.UVVector.(zero_field, zero_field)
    ClimaCore.Fields.level(Geometry.Covariant3Vector.(zero_field) .⊗ Geometry.Covariant12Vector.(uv_local), 1)
end

function atmos_init(::Type{FT}, tspan; stepper = SSPRK33(), nelements = (6, 10), npolynomial = 4, dt = 0.02) where {FT}
    params = DryBaroclinicWaveParameters()

    domain = SphericalShell(
        FT,
        radius = CLIMAParameters.Planet.planet_radius(params),
        height = FT(30.0e3),
        nelements = nelements,
        npolynomial = npolynomial,
    )

    # initiate boundary conditions
    boundary_conditions = (;
        ρe_tot = (
            top = NoFlux(),
            bottom = CouplerEnergyFlux(
                FT(0.001), # transfer coefficient, default: 1e-3 [dimenisonless]
                atmos_bcfield_init(domain),
            ), # surface specific enthalpy ~ (T_sfc * c_p)
        ),
        uh = (top = NoVectorFlux(), bottom = CouplerVectorFlux(FT(0.001), atmos_bcvectorfield_init(domain))),
    )

    model = Nonhydrostatic3DModel(
        domain = domain,
        boundary_conditions = boundary_conditions,
        parameters = params,
        hyperdiffusivity = FT(1e16),
        vertical_diffusion = ConstantViscosity(ν = FT(5)), #default: 5 m^2/s
    )

    # execute (using the regression test from ClimaAtmos)
    simulation = Simulation(model, stepper, dt = dt, tspan = tspan)

    # test set function
    (; ρ, uh, w, ρe_tot) = init_3d_baroclinic_wave(FT, params)
    set!(simulation, :base, ρ = ρ, uh = uh, w = w)
    set!(simulation, :thermodynamics, ρe_tot = ρe_tot)

    simulation
end
