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
    PrescribedIceSimulation{P, I}

Ice concentration is prescribed, and we solve the following energy equation:

    (h * ρ * c) d T_bulk dt = (-F_turb_energy + (1 - α) * SW_d + LW_d - LW_u) + F_conductive

    with
    F_conductive = k_ice (T_base - T_sfc) / (h)

    The bulk temperature (`T_bulk`) is the prognostic variable which is being
    modified by turbulent aerodynamic (`F_turb_energy`) and radiative (`F_turb_energy`) fluxes,
    as well as a conductive flux that depends on the temperature difference
    across the ice layer (with `T_base` being prescribed).

In the current version, the sea ice has a prescribed thickness. T_sfc is extrapolated from T_base
and T_bulk assuming the ice temperature varies linearly between the ice surface and the base.
"""
struct PrescribedIceSimulation{P, I} <: Interfacer.SeaIceModelSimulation
    params::P
    integrator::I
end

# sea-ice parameters
Base.@kwdef struct IceSlabParameters{FT <: AbstractFloat}
    h::FT                   # ice thickness [m]
    ρ::FT                   # density of sea ice [kg / m3]
    c::FT                   # specific heat of sea ice [J / kg / K]
    T_base::FT              # temperature of sea water at the ice base
    z0m::FT                 # roughness length for momentum [m]
    z0b::FT                 # roughness length for tracers [m]
    T_freeze::FT            # freezing temperature of sea water [K]
    k_ice::FT               # thermal conductivity of sea ice [W / m / K] (less in HM71)
    α::FT                   # albedo of sea ice, roughly tuned to match observations
    ϵ::FT                   # emissivity of sea ice
    σ::FT                   # Stefan-Boltzmann constant [W / m2 / K4]
end

"""
    IceSlabParameters{FT}(coupled_param_dict; h = FT(2), ρ = FT(900), c = FT(2100),
                          T_base = FT(271.2), z0m = FT(1e-4), z0b = FT(1e-4),
                          T_freeze = FT(271.2), k_ice = FT(2), α = FT(0.65), ϵ = FT(1))

Initialize the `IceSlabParameters` object with the coupled parameters.

# Arguments
- `coupled_param_dict`: a dictionary of coupled parameters (required)
- `h`: ice thickness [m] (default: 2)
- `ρ`: density of sea ice [kg / m3] (default: 900)
- `c`: specific heat of sea ice [J / kg / K] (default: 2100)
- `T_base`: temperature of sea water at the ice base [K] (default: 271.2)
- `z0m`: roughness length for momentum [m] (default: 1e-4)
- `z0b`: roughness length for tracers [m] (default: 1e-4)
- `T_freeze`: freezing temperature of sea water [K] (default: 271.2)
- `k_ice`: thermal conductivity of sea ice [W / m / K] (default: 2)
- `α`: albedo of sea ice (default: 0.65)
- `ϵ`: emissivity of sea ice (default: 1)

# Returns
- `IceSlabParameters{FT}`: an `IceSlabParameters` object
"""
function IceSlabParameters{FT}(
    coupled_param_dict;
    h = FT(2),
    ρ = FT(900),
    c = FT(2100),
    T_base = FT(271.2),
    z0m = FT(1e-4),
    z0b = FT(1e-4),
    T_freeze = FT(271.2),
    k_ice = FT(2),
    α = FT(0.65),
    ϵ = FT(1),
) where {FT}
    return IceSlabParameters{FT}(;
        h,
        ρ,
        c,
        T_base,
        z0m,
        z0b,
        T_freeze,
        k_ice,
        α,
        ϵ,
        σ = coupled_param_dict["stefan_boltzmann_constant"],
    )
end

# init simulation
function slab_ice_space_init(::Type{FT}, space, params) where {FT}
    # bulk temperatures commonly 10-20 K below freezing for sea ice 2m thick in winter,
    # and closer to freezing in summer and when melting.
    Y = CC.Fields.FieldVector(T_bulk = ones(space) .* params.T_freeze .- FT(5.0))
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
        start_date,
        land_fraction,
        stepper = CTS.RK4(),
        sic_path::Union{Nothing, String} = nothing,
    ) where {FT}

Initializes the `DiffEq` problem, and creates a Simulation-type object
containing the necessary information for `Interfacer.step!` in the coupling loop.

This model reads in prescribed sea ice concentration data and solves the energy equation
for the surface temperature of the sea ice. The sea ice concentration is updated
at each timestep.

The sea ice concentration is read from the file specified by `sic_path`. If `sic_path` is `nothing`,
the model will use the default path from `ClimaArtifacts`.

"""
function PrescribedIceSimulation(
    ::Type{FT};
    tspan,
    dt,
    saveat,
    space,
    coupled_param_dict,
    thermo_params,
    comms_ctx,
    start_date,
    land_fraction,
    stepper = CTS.RK4(),
    sic_path::Union{Nothing, String} = nothing,
) where {FT}
    # Set up prescribed sea ice concentration object
    sic_data =
        isnothing(sic_path) ?
        try
            joinpath(
                @clima_artifact("historical_sst_sic", comms_ctx),
                "MODEL.ICE.HAD187001-198110.OI198111-202206.nc",
            )
        catch error
            @warn "Using lowres SIC. If you want the higher resolution version, you have to obtain it from ClimaArtifacts"
            joinpath(
                @clima_artifact("historical_sst_sic_lowres", comms_ctx),
                "MODEL.ICE.HAD187001-198110.OI198111-202206_lowres.nc",
            )
        end : sic_path
    @info "PrescribedIce: using SIC file" sic_data

    SIC_timevaryinginput = TimeVaryingInput(
        sic_data,
        "SEAICE",
        space,
        reference_date = start_date,
        file_reader_kwargs = (; preprocess_func = (data) -> data / 100,), ## convert to fraction
    )

    # Get initial SIC values and use them to calculate ice fraction
    ice_fraction = CC.Fields.zeros(space)
    evaluate!(ice_fraction, SIC_timevaryinginput, tspan[1])

    # Make ice fraction binary rather than fractional
    ice_fraction = ifelse.(ice_fraction .> FT(0.5), FT(1), FT(0))

    params = IceSlabParameters{FT}(coupled_param_dict)

    Y = slab_ice_space_init(FT, space, params)
    cache = (;
        F_turb_energy = CC.Fields.zeros(space),
        SW_d = CC.Fields.zeros(space),
        LW_d = CC.Fields.zeros(space),
        area_fraction = ice_fraction,
        SIC_timevaryinginput = SIC_timevaryinginput,
        land_fraction = land_fraction,
        dt = dt,
        thermo_params = thermo_params,
        # add dss_buffer to cache to avoid runtime dss allocation
        dss_buffer = CC.Spaces.create_dss_buffer(Y),
    )

    ode_algo = CTS.ExplicitAlgorithm(stepper)
    ode_function = CTS.ClimaODEFunction(
        T_exp! = ice_rhs!,
        dss! = (Y, p, t) -> CC.Spaces.weighted_dss!(Y, p.dss_buffer),
    )
    if typeof(dt) isa Number
        dt = Float64(dt)
        tspan = Float64.(tspan)
        saveat = Float64.(saveat)
    end
    problem = SciMLBase.ODEProblem(ode_function, Y, tspan, (; cache..., params = params))
    integrator =
        SciMLBase.init(problem, ode_algo, dt = dt, saveat = saveat, adaptive = false)

    sim = PrescribedIceSimulation(params, integrator)

    # DSS state to ensure we have continuous fields
    dss_state!(sim)
    return sim
end

# extensions required by Interfacer
Interfacer.get_field(sim::PrescribedIceSimulation, ::Val{:area_fraction}) =
    sim.integrator.p.area_fraction
Interfacer.get_field(sim::PrescribedIceSimulation, ::Val{:ice_concentration}) =
    sim.integrator.p.area_fraction
Interfacer.get_field(sim::PrescribedIceSimulation, ::Val{:emissivity}) =
    sim.integrator.p.params.ϵ
Interfacer.get_field(sim::PrescribedIceSimulation, ::Val{:roughness_buoyancy}) =
    sim.integrator.p.params.z0b
Interfacer.get_field(sim::PrescribedIceSimulation, ::Val{:roughness_momentum}) =
    sim.integrator.p.params.z0m
Interfacer.get_field(
    sim::PrescribedIceSimulation,
    ::Union{Val{:surface_direct_albedo}, Val{:surface_diffuse_albedo}},
) = sim.integrator.p.params.α
Interfacer.get_field(sim::PrescribedIceSimulation, ::Val{:surface_temperature}) =
    ice_surface_temperature.(sim.integrator.u.T_bulk, sim.integrator.p.params.T_base)

function Interfacer.get_field(sim::PrescribedIceSimulation, ::Val{:beta})
    # assume no LHF over sea ice
    FT = eltype(sim.integrator.u)
    return FT(0)
end

# Approximates the surface temperature of the sea ice assuming
# the ice temperature varies linearly between the ice surface and the base
ice_surface_temperature(T_bulk, T_base) = 2 * T_bulk - T_base

"""
    Interfacer.get_field(sim::PrescribedIceSimulation, ::Val{:energy})

Extension of Interfacer.get_field to get the energy of the ocean.
It multiplies the the slab temperature by the heat capacity, density, and depth.
"""
Interfacer.get_field(sim::PrescribedIceSimulation, ::Val{:energy}) =
    sim.integrator.p.params.ρ .* sim.integrator.p.params.c .* sim.integrator.u.T_bulk .*
    sim.integrator.p.params.h

function Interfacer.update_field!(
    sim::PrescribedIceSimulation,
    ::Val{:area_fraction},
    field::CC.Fields.Field,
)
    sim.integrator.p.area_fraction .= field
    return nothing
end
function Interfacer.update_field!(sim::PrescribedIceSimulation, ::Val{:SW_d}, field)
    Interfacer.remap!(sim.integrator.p.SW_d, field)
end
function Interfacer.update_field!(sim::PrescribedIceSimulation, ::Val{:LW_d}, field)
    Interfacer.remap!(sim.integrator.p.LW_d, field)
end
function Interfacer.update_field!(
    sim::PrescribedIceSimulation,
    ::Val{:turbulent_energy_flux},
    field,
)
    Interfacer.remap!(sim.integrator.p.F_turb_energy, field)
end
Interfacer.update_field!(
    sim::PrescribedIceSimulation,
    ::Val{:turbulent_moisture_flux},
    field,
) = nothing

# extensions required by FieldExchanger
Interfacer.step!(sim::PrescribedIceSimulation, t) =
    Interfacer.step!(sim.integrator, t - sim.integrator.t, true)

function FluxCalculator.update_turbulent_fluxes!(
    sim::PrescribedIceSimulation,
    fields::NamedTuple,
)
    Interfacer.update_field!(sim, Val(:turbulent_energy_flux), fields.F_lh .+ fields.F_sh)
    return nothing
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
    (; k_ice, h, T_base, ρ, c, ϵ, α, T_freeze, σ) = p.params

    # Update the cached area fraction with the current SIC
    evaluate!(p.area_fraction, p.SIC_timevaryinginput, t)
    @. p.area_fraction = ifelse(p.area_fraction > FT(0.5), FT(1), FT(0))

    # Overwrite ice fraction with the static land area fraction anywhere we have nonzero land area
    #  max needed to avoid Float32 errors (see issue #271; Heisenbug on HPC)
    @. p.area_fraction = max(min(p.area_fraction, FT(1) - p.land_fraction), FT(0))

    # Calculate the conductive flux, and set it to zero if the area fraction is zero
    F_conductive = @. k_ice / (h) * (T_base - ice_surface_temperature(Y.T_bulk, T_base)) # fluxes are defined to be positive when upward
    rhs = @. (
        -p.F_turb_energy +
        (1 - α) * p.SW_d +
        ϵ * (p.LW_d - σ * ice_surface_temperature(Y.T_bulk, T_base)^4) +
        F_conductive
    ) / (h * ρ * c)
    # Zero out tendencies where there is no ice, so that ice temperature remains constant there
    @. rhs = ifelse(p.area_fraction ≈ 0, zero(rhs), rhs)

    # If tendencies lead to temperature above freezing, set temperature to freezing
    @. dY.T_bulk = min(rhs, (T_freeze - Y.T_bulk) / float(p.dt))
end

"""
    dss_state!(sim::PrescribedIceSimulation)

Perform DSS on the state of a component simulation, intended to be used
before the initial step of a run. This method acts on prescribed ice simulations.
"""
dss_state!(sim::PrescribedIceSimulation) =
    CC.Spaces.weighted_dss!(sim.integrator.u, sim.integrator.p.dss_buffer)

function Checkpointer.get_model_cache(sim::PrescribedIceSimulation)
    return sim.integrator.p
end

function Checkpointer.restore_cache!(sim::PrescribedIceSimulation, new_cache)
    old_cache = Checkpointer.get_model_cache(sim)
    for p in propertynames(old_cache)
        if getproperty(old_cache, p) isa Field
            ArrayType = ClimaComms.array_type(getproperty(old_cache, p))
            parent(getproperty(old_cache, p)) .=
                ArrayType(parent(getproperty(new_cache, p)))
        end
    end
end
