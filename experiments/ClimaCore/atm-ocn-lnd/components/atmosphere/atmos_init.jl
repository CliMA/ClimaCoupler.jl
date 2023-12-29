using ClimaCore.Geometry
using ClimaCore: DataLayouts, Operators, Geometry
using ClimaCore: Fields, Domains, Topologies, Meshes, Spaces

using IntervalSets
import SciMLBase: ODEProblem, solve, step!, init, reinit!
import ClimaTimeSteppers as CTS
using Logging: global_logger
using TerminalLoggers: TerminalLogger

using RecursiveArrayTools
using LinearAlgebra
using Random

# define model equations:
include("atmos_rhs.jl")

parameters = (
    # domain
    zmin_atm = FT(0.0), # height of atm stack bottom [m]
    zmax_atm = FT(1.0), # height of atm stack top [m]
    n = 15,  # number of elements in atm stack

    # atmos physics (overwriting CliMAParameters defaults?)
    MSLP = 1e5, # mean sea level pressure [Pa]
    grav = 9.8, # gravitational acceleration [m /s^2]
    R_d = 287.058, # R dry (gas constant / mol mass dry air)  [J / K / kg]
    C_p = 287.058 * 1.4 / (1.4 - 1), # heat capacity at constant pressure [J / K / kg]
    C_v = 287.058 / (1.4 - 1), # heat capacity at constant volume [J / K / kg]
    R_m = 87.058, # moist R, assumed to be dry [J / K / kg]
    f = 7.29e-5, # Coriolis parameters [1/s]
    ν = 0.1, #0.01 # viscosity, diffusivity
    Ch = 0.0015, # bulk transfer coefficient for sensible heat
    Cd = 0.01 / (2e2 / 30.0), #drag coeff
    ug = 1.0,
    vg = 0.0,
    d = sqrt(2.0 * 0.01 / 5e-5), #?

    # radiation parameters for DryBulkFormulaWithRadiation() SurfaceFluxType
    τ = 0.9,     # atmospheric transmissivity
    α = 0.5,     # surface albedo
    σ = 5.67e-8, # Steffan Boltzmann constant [kg / s^3 / K^4]
    g_a = 0.06,    # aerodynamic conductance for heat transfer [kg / m^2 / s]
    ϵ = 0.98,    # broadband emissivity / absorptivity of the surface
    F_a = 0.0,     # downward LW flux from the atmosphere [W / m^2]
    F_sol = 1361,   # incoming solar TOA radiation [W / m^2]
    τ_d = 10,     # idealized daily cycle period [s]

    # surface fluxes
    λ = FT(0.01),#FT(1e-5)    # coupling transfer coefficient for LinearRelaxation() SurfaceFluxType
)


function atmos_init(T_sfc_init;
    Lz,
    Nz,
    minimum_reference_temperature = 230.0, # K
    start_time = 0.0,
    stop_time = 1.0,
    Δt_min = 0.02,
)

    # Get surface temperature from the land simulation
    surface_temperature = T_sfc_init

    ########
    # Set up atmos domain
    ########
    domain_atm = Domains.IntervalDomain(
        Geometry.ZPoint{FT}(parameters.zmin_atm),
        Geometry.ZPoint{FT}(parameters.zmax_atm);
        boundary_tags = (:bottom, :top),
    );
    mesh_atm = Meshes.IntervalMesh(domain_atm, nelems = parameters.n); # struct, allocates face boundaries to 5,6: atmos
    center_space_atm = Spaces.CenterFiniteDifferenceSpace(mesh_atm); # collection of the above, discretises space into FD and provides coords
    face_space_atm = Spaces.FaceFiniteDifferenceSpace(center_space_atm)

    """ Initialize fields located at cell centers in the vertical. """
    function init_centers(zc, parameters)
        (; grav, C_p, MSLP, R_d) = parameters

        # temperature
        Γ = grav / C_p
        T = max(surface_temperature - Γ * zc, minimum_reference_temperature)

        # pressure
        p = MSLP * (T / surface_temperature)^(grav / (R_d * Γ))

        if T == minimum_reference_temperature
            z_top = (surface_temperature - minimum_reference_temperature) / Γ
            H_min = R_d * minimum_reference_temperature / grav
            p *= exp(-(zc - z_top) / H_min)
        end

        # potential temperature
        θ = surface_temperature

        # density
        ρ = p / (R_d * θ * (p / MSLP)^(R_d / C_p))

        # velocties
        u = 1.0
        v = 0.0

        return (ρ = ρ, u = u, v = v, ρθ = ρ * θ)
    end

    """ Initialize fields located at cell interfaces in the vertical. """
    function init_faces(zf, parameters)
        return (; w = 0.0 .* zf)
    end

    # Initialize the atmospheric states Yc and Yf
    z_centers = Fields.coordinate_field(center_space_atm).z
    z_faces = Fields.coordinate_field(face_space_atm).z
    Yc = init_centers.(z_centers, Ref(parameters))
    Yf = init_faces.(z_faces, Ref(parameters))
    T_atm_0 = (Yc = Yc, Yf = Yf)

    # Put all prognostic variable arrays into a vector and ensure that solve can partition them
    Y_atm = ArrayPartition((T_atm_0.Yc, T_atm_0.Yf, zeros(3)))

    ode_algo = CTS.ExplicitAlgorithm(CTS.RK4())
    ode_function = CTS.ClimaODEFunction(T_exp! = ∑tendencies_atm!)

    problem = ODEProblem(ode_function, Y_atm, (start_time, stop_time), (parameters, [surface_temperature]))
    simulation = init(problem, ode_algo, dt = Δt_min, saveat = 1 * Δt_min, adaptive = false)

    return simulation
end
