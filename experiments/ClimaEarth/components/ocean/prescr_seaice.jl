import SciMLBase
import ClimaCore as CC
import ClimaTimeSteppers as CTS
import ClimaUtilities.TimeVaryingInputs: TimeVaryingInput, evaluate!
import ClimaUtilities.ClimaArtifacts: @clima_artifact
import Interpolations # triggers InterpolationsExt in ClimaUtilities
import Thermodynamics as TD
import ClimaCoupler: Checkpointer, FluxCalculator, Interfacer, Utilities

###
### Functions required by ClimaCoupler.jl for a SurfaceModelSimulation
###
"""
    PrescribedIceSimulation{P, Y, D, I}

Ice concentration is prescribed, and we solve the following energy equation:

    (h * ρ * c) d T_sfc dt = -(F_turb_energy + F_radiative) + F_conductive

    with
    F_conductive = k_ice (T_base - T_sfc) / (h)

    The surface temperature (`T_sfc`) is the prognostic variable which is being
    modified by turbulent aerodynamic (`F_turb_energy`) and radiative (`F_turb_energy`) fluxes,
    as well as a conductive flux that depends on the temperature difference
    across the ice layer (with `T_base` being prescribed).

In the current version, the sea ice has a prescribed thickness.
"""
struct PrescribedIceSimulation{P, Y, D, I} <: Interfacer.SeaIceModelSimulation
    params::P
    Y_init::Y
    domain::D
    integrator::I
end
Interfacer.name(::PrescribedIceSimulation) = "PrescribedIceSimulation"

# sea-ice parameters
Base.@kwdef struct IceSlabParameters{FT <: AbstractFloat}
    h::FT = 2               # ice thickness [m]
    ρ::FT = 900             # density of sea ice [kg / m3]
    c::FT = 2100            # specific heat of sea ice [J / kg / K]
    T_base::FT = 271.2      # temperature of sea water at the ice base
    z0m::FT = 1e-4          # roughness length for momentum [m]
    z0b::FT = 1e-4          # roughness length for tracers [m]
    T_freeze::FT = 271.2    # freezing temperature of sea water [K]
    k_ice::FT = 2           # thermal conductivity of sea ice [W / m / K] (less in HM71)
    α::FT = 0.65            # albedo of sea ice, roughly tuned to match observations
end

Interfacer.name(::IceSlabParameters) = "IceSlabParameters"

# init simulation
function slab_ice_space_init(::Type{FT}, space, params) where {FT}
    Y = CC.Fields.FieldVector(T_sfc = ones(space) .* params.T_freeze)
    return Y
end

"""
    PrescribedIceSimulation(
        ::Type{FT};
        tspan,
        dt,
        saveat,
        space,
        thermo_params,
        comms_ctx,
        date0,
        mono_surface,
        land_fraction,
        stepper = CTS.RK4()
    ) where {FT}

Initializes the `DiffEq` problem, and creates a Simulation-type object
containing the necessary information for `Interfacer.step!` in the coupling loop.

This model reads in prescribed sea ice concentration data and solves the energy equation
for the surface temperature of the sea ice. The sea ice concentration is updated
at each timestep.

"""
function PrescribedIceSimulation(
    ::Type{FT};
    tspan,
    dt,
    saveat,
    space,
    thermo_params,
    comms_ctx,
    date0,
    mono_surface,
    land_fraction,
    stepper = CTS.RK4(),
) where {FT}
    # Set up prescribed sea ice concentration object
    sic_data = try
        joinpath(@clima_artifact("historical_sst_sic", comms_ctx), "MODEL.ICE.HAD187001-198110.OI198111-202206.nc")
    catch error
        @warn "Using lowres SIC. If you want the higher resolution version, you have to obtain it from ClimaArtifacts"
        joinpath(
            @clima_artifact("historical_sst_sic_lowres", comms_ctx),
            "MODEL.ICE.HAD187001-198110.OI198111-202206_lowres.nc",
        )
    end

    SIC_timevaryinginput = TimeVaryingInput(
        sic_data,
        "SEAICE",
        space,
        reference_date = date0,
        file_reader_kwargs = (; preprocess_func = (data) -> data / 100,), ## convert to fraction
    )

    # Get initial SIC values and use them to calculate ice fraction
    SIC_init = CC.Fields.zeros(space)
    evaluate!(SIC_init, SIC_timevaryinginput, tspan[1])
    ice_fraction = get_ice_fraction.(SIC_init, mono_surface)

    # Overwrite ice fraction with the static land area fraction anywhere we have nonzero land area
    #  max needed to avoid Float32 errors (see issue #271; Heisenbug on HPC)
    @. ice_fraction = max(min(ice_fraction, FT(1) - land_fraction), FT(0))

    params = IceSlabParameters{FT}()

    Y = slab_ice_space_init(FT, space, params)
    cache = (;
        F_turb_energy = CC.Fields.zeros(space),
        F_radiative = CC.Fields.zeros(space),
        q_sfc = CC.Fields.zeros(space),
        ρ_sfc = CC.Fields.zeros(space),
        area_fraction = ice_fraction,
        SIC_timevaryinginput = SIC_timevaryinginput,
        land_fraction = land_fraction,
        dt = dt,
        thermo_params = thermo_params,
        # add dss_buffer to cache to avoid runtime dss allocation
        dss_buffer = CC.Spaces.create_dss_buffer(Y),
    )

    ode_algo = CTS.ExplicitAlgorithm(stepper)
    ode_function = CTS.ClimaODEFunction(T_exp! = ice_rhs!, dss! = (Y, p, t) -> CC.Spaces.weighted_dss!(Y, p.dss_buffer))

    problem = SciMLBase.ODEProblem(ode_function, Y, Float64.(tspan), (; cache..., params = params))
    integrator = SciMLBase.init(problem, ode_algo, dt = Float64(dt), saveat = Float64(saveat), adaptive = false)

    sim = PrescribedIceSimulation(params, Y, space, integrator)

    # DSS state to ensure we have continuous fields
    dss_state!(sim)
    return sim
end

# extensions required by Interfacer
Interfacer.get_field(sim::PrescribedIceSimulation, ::Val{:air_density}) = sim.integrator.p.ρ_sfc
Interfacer.get_field(sim::PrescribedIceSimulation, ::Val{:area_fraction}) = sim.integrator.p.area_fraction
Interfacer.get_field(sim::PrescribedIceSimulation, ::Val{:beta}) = convert(eltype(sim.integrator.u), 1.0)
Interfacer.get_field(sim::PrescribedIceSimulation, ::Val{:roughness_buoyancy}) = sim.integrator.p.params.z0b
Interfacer.get_field(sim::PrescribedIceSimulation, ::Val{:roughness_momentum}) = sim.integrator.p.params.z0m
Interfacer.get_field(sim::PrescribedIceSimulation, ::Union{Val{:surface_direct_albedo}, Val{:surface_diffuse_albedo}}) =
    sim.integrator.p.params.α
Interfacer.get_field(sim::PrescribedIceSimulation, ::Val{:surface_humidity}) = sim.integrator.p.q_sfc
Interfacer.get_field(sim::PrescribedIceSimulation, ::Val{:surface_temperature}) = sim.integrator.u.T_sfc
Interfacer.get_field(sim::PrescribedIceSimulation, ::Val{:water}) = nothing

"""
    Interfacer.get_field(sim::PrescribedIceSimulation, ::Val{:energy})

Extension of Interfacer.get_field to get the energy of the ocean.
It multiplies the the slab temperature by the heat capacity, density, and depth.
"""
Interfacer.get_field(sim::PrescribedIceSimulation, ::Val{:energy}) =
    sim.integrator.p.params.ρ .* sim.integrator.p.params.c .* sim.integrator.u.T_sfc .* sim.integrator.p.params.h

function Interfacer.update_field!(sim::PrescribedIceSimulation, ::Val{:air_density}, field)
    parent(sim.integrator.p.ρ_sfc) .= parent(field)
end
function Interfacer.update_field!(sim::PrescribedIceSimulation, ::Val{:area_fraction}, field::CC.Fields.Field)
    sim.integrator.p.area_fraction .= field
end
function Interfacer.update_field!(sim::PrescribedIceSimulation, ::Val{:radiative_energy_flux_sfc}, field)
    parent(sim.integrator.p.F_radiative) .= parent(field)
end
function Interfacer.update_field!(sim::PrescribedIceSimulation, ::Val{:turbulent_energy_flux}, field)
    parent(sim.integrator.p.F_turb_energy) .= parent(field)
end
Interfacer.update_field!(sim::PrescribedIceSimulation, ::Val{:turbulent_moisture_flux}, field) = nothing

# extensions required by FieldExchanger
Interfacer.step!(sim::PrescribedIceSimulation, t) = Interfacer.step!(sim.integrator, t - sim.integrator.t, true)
Interfacer.reinit!(sim::PrescribedIceSimulation) = Interfacer.reinit!(sim.integrator)

# extensions required by FluxCalculator (partitioned fluxes)
function FluxCalculator.update_turbulent_fluxes!(sim::PrescribedIceSimulation, fields::NamedTuple)
    (; F_turb_energy) = fields
    @. sim.integrator.p.F_turb_energy = F_turb_energy
end

"""
    Checkpointer.get_model_prog_state(sim::PrescribedIceSimulation)

Extension of Checkpointer.get_model_prog_state to get the model state.
"""
function Checkpointer.get_model_prog_state(sim::PrescribedIceSimulation)
    return sim.integrator.u
end

###
### Sea ice model-specific functions (not explicitly required by ClimaCoupler.jl)
###
# setting that SIC < 0.5 is counted as ocean if binary remapping.
get_ice_fraction(h_ice::FT, mono::Bool, threshold = 0.5) where {FT} =
    mono ? h_ice : Utilities.binary_mask(h_ice, threshold)

"""
    ice_rhs!(dY, Y, p, t)

Rhs method in the form as required by `ClimeTimeSteppers`, with the tendency vector `dY`,
the state vector `Y`, the parameter vector `p`, and the simulation time `t` as input arguments.

This sea-ice energy formulation follows [Holloway and Manabe 1971](https://journals.ametsoc.org/view/journals/mwre/99/5/1520-0493_1971_099_0335_socbag_2_3_co_2.xml?tab_body=pdf),
where sea-ice concentrations and thicknes are prescribed, and the model solves
for temperature (curbed at the freezing point).
"""
function ice_rhs!(dY, Y, p, t)
    FT = eltype(Y)
    params = p.params

    # Update the cached area fraction with the current SIC
    evaluate!(p.area_fraction, p.SIC_timevaryinginput, t)

    # Overwrite ice fraction with the static land area fraction anywhere we have nonzero land area
    #  max needed to avoid Float32 errors (see issue #271; Heisenbug on HPC)
    @. p.area_fraction = max(min(p.area_fraction, FT(1) - p.land_fraction), FT(0))

    F_conductive = @. params.k_ice / (params.h) * (params.T_base - Y.T_sfc) # fluxes are defined to be positive when upward
    rhs = @. (-p.F_turb_energy - p.F_radiative + F_conductive) / (params.h * params.ρ * params.c)
    # If tendencies lead to temperature above freezing, set temperature to freezing
    @. rhs = min(rhs, (params.T_freeze - Y.T_sfc) / p.dt)
    # mask out no-ice areas
    area_mask = Utilities.binary_mask.(p.area_fraction)
    dY.T_sfc .= rhs .* area_mask

    @. p.q_sfc = TD.q_vap_saturation_generic.(p.thermo_params, Y.T_sfc, p.ρ_sfc, TD.Ice())
end

"""
    dss_state!(sim::PrescribedIceSimulation)

Perform DSS on the state of a component simulation, intended to be used
before the initial step of a run. This method acts on prescribed ice simulations.
"""
dss_state!(sim::PrescribedIceSimulation) = CC.Spaces.weighted_dss!(sim.integrator.u, sim.integrator.p.dss_buffer)
