import ClimaCore as CC
import ClimaTimeSteppers as CTS
import ClimaUtilities.TimeManager: ITime
import ClimaUtilities.TimeVaryingInputs: TimeVaryingInput, evaluate!
import ClimaUtilities.ClimaArtifacts: @clima_artifact
import Interpolations # triggers InterpolationsExt in ClimaUtilities
import Thermodynamics as TD
import ..Checkpointer, ..FluxCalculator, ..Interfacer, ..Utilities, ..TimeManager
import ClimaComms
import NCDatasets

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
struct PrescribedIceSimulation{P, I} <: Interfacer.AbstractSeaIceSimulation
    params::P
    integrator::I
end

# sea-ice parameters
Base.@kwdef struct IceSlabParameters{FT <: AbstractFloat}
    h_nh::FT                # ice thickness in Northern Hemisphere [m]
    h_sh::FT                # ice thickness in Southern Hemisphere [m]
    ρ::FT                   # density of sea ice [kg / m3]
    c::FT                   # specific heat of sea ice [J / kg / K]
    T_base::FT              # temperature of sea water at the ice base [K]
    z0m::FT                 # roughness length for momentum [m]
    z0b::FT                 # roughness length for tracers [m]
    T_freeze::FT            # freezing/melting temperature of sea-ice surface [K]
    k_ice::FT               # thermal conductivity of sea ice [W / m / K] (less in HM71)
    α_cold::FT              # cold snow-covered sea-ice albedo
    α_melt::FT              # sea-ice albedo with melt ponds at T_freeze
    ΔT_melt::FT             # temperature range below T_freeze over which albedo ramps [K]
    ϵ::FT                   # emissivity of sea ice
    σ::FT                   # Stefan-Boltzmann constant [W / m2 / K4]
    T_sfc_min::FT           # minimum allowed surface temperature [K]
end

Base.broadcastable(p::IceSlabParameters) = tuple(p)

"""
    IceSlabParameters{FT}(coupled_param_dict; h_nh = FT(2), h_sh = FT(1),
                          ρ = FT(900), c = FT(2100), T_base = FT(271.2),
                          z0m = FT(1e-4), z0b = FT(1e-4), T_freeze = FT(273.15),
                          k_ice = FT(2), α_cold = FT(0.78), α_melt = FT(0.55),
                          ΔT_melt = FT(3), ϵ = FT(1), T_sfc_min = FT(200))

Initialize the `IceSlabParameters` object with the coupled parameters.

# Arguments
- `coupled_param_dict`: a dictionary of coupled parameters (required)
- `h_nh`: ice thickness in the Northern Hemisphere (`lat ≥ 0`) [m] (default: 2)
- `h_sh`: ice thickness in the Southern Hemisphere (`lat < 0`) [m] (default: 1)
- `ρ`: density of sea ice [kg / m3] (default: 900)
- `c`: specific heat of sea ice [J / kg / K] (default: 2100)
- `T_base`: temperature of sea water at the ice base [K] (default: 271.2)
- `z0m`: roughness length for momentum [m] (default: 1e-4)
- `z0b`: roughness length for tracers [m] (default: 1e-4)
- `T_freeze`: melting temperature of sea-ice surface [K] (default: 273.15)
- `k_ice`: thermal conductivity of sea ice [W / m / K] (default: 2)
- `α_cold`: cold snow-covered sea-ice albedo for `T_sfc ≤ T_freeze - ΔT_melt`
  (default: 0.78)
- `α_melt`: sea-ice albedo with melt ponds at `T_sfc = T_freeze` (default: 0.55)
- `ΔT_melt`: melt-ramp width below `T_freeze` [K] (default: 3)
- `ϵ`: emissivity of sea ice (default: 1)
- `T_sfc_min`: minimum allowed surface temperature [K] (default: 200)

# Returns
- `IceSlabParameters{FT}`: an `IceSlabParameters` object
"""
function IceSlabParameters{FT}(
    coupled_param_dict;
    h_nh = FT(2),
    h_sh = FT(1),
    ρ = FT(900),
    c = FT(2100),
    T_base = FT(271.2),
    z0m = FT(1e-4),
    z0b = FT(1e-4),
    T_freeze = FT(273.15),
    k_ice = FT(2),
    α_cold = FT(0.78),
    α_melt = FT(0.55),
    ΔT_melt = FT(3),
    ϵ = FT(1),
    T_sfc_min = FT(200),
) where {FT}
    return IceSlabParameters{FT}(;
        h_nh,
        h_sh,
        ρ,
        c,
        T_base,
        z0m,
        z0b,
        T_freeze,
        k_ice,
        α_cold,
        α_melt,
        ΔT_melt,
        ϵ,
        σ = coupled_param_dict["stefan_boltzmann_constant"],
        T_sfc_min,
    )
end

"""
    ice_thickness(lat, params::IceSlabParameters)

Prescribed sea-ice slab thickness [m] at latitude `lat` [degrees]: `params.h_nh`
in the Northern Hemisphere (`lat ≥ 0`) and `params.h_sh` in the Southern
Hemisphere (`lat < 0`).
"""
ice_thickness(lat::FT, params::IceSlabParameters{FT}) where {FT} =
    ifelse(lat < zero(FT), params.h_sh, params.h_nh)

"""
    slab_ice_space_init(::Type{FT}, space, params, sic_data, start_date, t_start, ice_fraction)

Initialize ice bulk temperature from the ISTL1 variable in `sic_data` if available.
ISTL1 from ERA5 (ice surface temperature layer 1, 0-7cm) is used as a proxy for surface temperature,
then back-calculates T_bulk assuming a linear temperature profile: T_bulk = (T_sfc + T_base) / 2.

Only uses ISTL1 temperatures where there is ice (ice_fraction > 0). In ice-free regions,
uses the fallback temperature to avoid initializing with values at/above freezing which
can cause numerical issues.

Falls back to constant initialization (T_freeze - 5K) if ISTL1 is not present or data is not finite.
"""
function slab_ice_space_init(
    ::Type{FT},
    space,
    params,
    sic_data,
    start_date,
    t_start,
    ice_fraction,
) where {FT}
    (; T_base, T_freeze) = params
    T_fallback = T_base - FT(5.0)

    # Initialize ice temperature if available
    has_ISTL1 = NCDatasets.NCDataset(sic_data, "r") do ds
        haskey(ds, "ISTL1")
    end

    if has_ISTL1
        # Use ERA5 ISTL1 (0-7cm, near-surface layer) as a proxy for surface temperature
        ISTL1_input = TimeVaryingInput(sic_data, "ISTL1", space; start_date)
        T_sfc_data = CC.Fields.zeros(space)
        evaluate!(T_sfc_data, ISTL1_input, t_start)

        # Back-calculate T_bulk from T_sfc assuming linear temperature profile
        T_bulk = @. (T_sfc_data + T_base) / FT(2)

        # Only use ISTL1 temperatures where there is ice; use fallback elsewhere
        @. T_bulk = ifelse(ice_fraction > 0, T_bulk, T_fallback)

        # Fall back to constant where data is not finite
        @. T_bulk = ifelse(isfinite(T_bulk), T_bulk, T_fallback)

        @info "PrescribedIce: initialized T_bulk from ISTL1"
        T_bulk_init = T_bulk
    else
        @info "PrescribedIce: ISTL1 not in file, using constant initialization" file =
            sic_data
        T_bulk_init = ones(space) .* T_fallback
    end

    @. T_bulk_init = ifelse(isnan(T_bulk_init), T_fallback, T_bulk_init)

    Y = CC.Fields.FieldVector(T_bulk = T_bulk_init)
    return Y
end

"""
    slab_ice_space_init(::Type{FT}, space, params, ::Nothing, args...)

Fallback initialization when no SIC data file is provided.
Initializes bulk temperature to 5K below freezing.
"""
function slab_ice_space_init(::Type{FT}, space, params, ::Nothing, args...) where {FT}
    # bulk temperatures commonly 10-20 K below freezing for sea ice 2m thick in winter,
    # and closer to freezing in summer and when melting.
    Y = CC.Fields.FieldVector(T_bulk = ones(space) .* params.T_freeze .- FT(5.0))
    return Y
end

"""
    Interfacer.SeaIceSimulation(::Type{FT}, ::Val{:prescribed}; kwargs...)

Extension of the generic SeaIceSimulation constructor for the prescribed sea ice model.
"""
function Interfacer.SeaIceSimulation(::Type{FT}, ::Val{:prescribed}; kwargs...) where {FT}
    return PrescribedIceSimulation(FT; kwargs...)
end

"""
    PrescribedIceSimulation(
        ::Type{FT};
        tspan,
        dt,
        saveat,
        boundary_space,
        thermo_params,
        comms_ctx,
        start_date,
        land_fraction,
        stepper = CTS.RK4(),
        sic_path::Union{Nothing, String} = nothing,
        binary_area_fraction::Bool = true,
        extra_kwargs...,
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
    boundary_space,
    coupled_param_dict,
    thermo_params,
    comms_ctx,
    start_date,
    land_fraction,
    stepper = CTS.RK4(),
    sic_path::Union{Nothing, String} = nothing,
    binary_area_fraction::Bool = true,
    domain_type::String = "global",
    extra_kwargs...,
) where {FT}
    # Set up prescribed sea ice concentration object
    sic_data =
        isnothing(sic_path) ?
        try
            joinpath(
                @clima_artifact("historical_sst_sic", comms_ctx),
                "MODEL.ICE.HAD187001-198110.OI198111-202206.nc",
            )
        catch
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
        boundary_space;
        start_date,
        file_reader_kwargs = (; preprocess_func = (data) -> data / 100,), ## convert to fraction
    )

    # Get initial SIC values and use them to calculate ice fraction
    ice_fraction = CC.Fields.zeros(boundary_space)
    if domain_type == "column"
        ice_fraction .= FT(1)
    else
        evaluate!(ice_fraction, SIC_timevaryinginput, tspan[1])

        # Ensure ice fraction is finite and not NaN
        @. ice_fraction = ifelse(isfinite(ice_fraction), ice_fraction, FT(0))
        # binary ice fraction (threshold at 0.5)
        if binary_area_fraction
            @. ice_fraction = ifelse(ice_fraction > FT(0.5), FT(1), FT(0))
        end
        # max/min needed to avoid Float32 errors (see issue #271; Heisenbug on HPC)
        @. ice_fraction = max(min(ice_fraction, FT(1) - land_fraction), FT(0))
    end

    params = IceSlabParameters{FT}(coupled_param_dict)

    # Initialize ice temperature: use ERA5 ice layer data if available, otherwise constant
    Y = slab_ice_space_init(
        FT,
        boundary_space,
        params,
        sic_data,
        start_date,
        tspan[1],
        ice_fraction,
    )

    cache = (;
        F_turb_energy = CC.Fields.zeros(boundary_space),
        SW_d = CC.Fields.zeros(boundary_space),
        LW_d = CC.Fields.zeros(boundary_space),
        α_sfc = CC.Fields.zeros(boundary_space),
        h_ice = ice_thickness.(CC.Fields.coordinate_field(boundary_space).lat, params),
        area_fraction = ice_fraction,
        SIC_timevaryinginput = SIC_timevaryinginput,
        land_fraction = land_fraction,
        binary_area_fraction = binary_area_fraction,
        domain_type = domain_type,
        dt = dt,
        thermo_params = thermo_params,
        # add dss_buffer to cache to avoid runtime dss allocation
        dss_buffer = Utilities.init_dss_buffer(Y),
    )

    ode_algo = CTS.ExplicitAlgorithm(stepper)
    ode_function = CTS.ClimaODEFunction(
        T_exp! = ice_rhs!,
        dss! = (Y, p, t) -> Utilities.apply_dss!(Y, p.dss_buffer),
    )
    if dt isa Number
        dt = Float64(dt)
        tspan = Float64.(tspan)
        saveat = Float64.(saveat)
    end

    # Create a schedule to update the SIC daily. The data provides monthly
    # averages, so it doesn't make sense to update the SIC more frequently than daily.
    SIC_schedule = TimeManager.calendar_dt_schedule("1days", start_date, tspan[1])
    SIC_update_cb =
        CTS.DiscreteCallback((_, _, integrator) -> SIC_schedule(integrator), read_sic_data!)

    problem = CTS.ODEProblem(ode_function, Y, tspan, (; cache..., params = params))
    integrator = CTS.init(
        problem,
        ode_algo,
        dt = dt,
        saveat = saveat,
        adaptive = false,
        callback = CTS.CallbackSet(SIC_update_cb),
    )

    sim = PrescribedIceSimulation(params, integrator)

    # DSS state to ensure we have continuous fields
    dss_state!(sim)
    update_ice_albedo!(sim.integrator.u, sim.integrator.p)
    return sim
end

"""
    read_sic_data!(integrator)

Read in the sea ice concentration data at the current time.
This function is intended to be used within a callback.
"""
function read_sic_data!(integrator)
    evaluate!(integrator.p.area_fraction, integrator.p.SIC_timevaryinginput, integrator.t)
end

"""
    ice_surface_albedo(T_sfc, T_freeze, ΔT_melt, α_cold, α_melt)

Sea-ice broadband surface albedo as a continuous linear ramp between cold
snow-covered sea ice (`α_cold` for `T_sfc ≤ T_freeze - ΔT_melt`) and summer
ice with melt ponds (`α_melt` at `T_sfc = T_freeze`).

The piecewise-linear temperature ramp follows the standard single-slab GCM
parameterization of [Ingram et al. (1989)](https://doi.org/10.1029/JD094iD06p08609)
and [Briegleb et al. (2004)](https://doi.org/10.5065/D65Q4T4Z), with default
effective pack-ice parameters (`α_cold = 0.78`, `α_melt = 0.55`, `ΔT_melt = 3 K`)
calibrated against CERES EBAF surface radiative fluxes.
"""
function ice_surface_albedo(
    T_sfc::FT,
    T_freeze::FT,
    ΔT_melt::FT,
    α_cold::FT,
    α_melt::FT,
) where {FT}
    return α_melt +
           (α_cold - α_melt) * clamp((T_freeze - T_sfc) / ΔT_melt, zero(FT), one(FT))
end

"""
    update_ice_albedo!(u, p)

Update the cached sea-ice surface albedo `p.α_sfc` in-place from the current
prognostic bulk temperature `u.T_bulk`.
"""
function update_ice_albedo!(u, p)
    (; T_base, T_sfc_min, T_freeze, ΔT_melt, α_cold, α_melt) = p.params
    @. p.α_sfc = ice_surface_albedo(
        ice_surface_temperature(u.T_bulk, T_base, T_sfc_min, T_freeze),
        T_freeze,
        ΔT_melt,
        α_cold,
        α_melt,
    )
    return nothing
end

"""
    Interfacer.step!(sim::PrescribedIceSimulation, t)

Advance `sim.integrator` to coupling time `t` using the cached surface albedo
`p.α_sfc`, then refresh `p.α_sfc` once at the coupling step so the atmosphere
and sea ice see the same absorbed shortwave flux.
"""
function Interfacer.step!(sim::PrescribedIceSimulation, t::Float64)
    @invoke Interfacer.step!(sim::Interfacer.AbstractComponentSimulation, t::Float64)
    update_ice_albedo!(sim.integrator.u, sim.integrator.p)
    return nothing
end
function Interfacer.step!(sim::PrescribedIceSimulation, t::ITime)
    @invoke Interfacer.step!(sim::Interfacer.AbstractComponentSimulation, t::ITime)
    update_ice_albedo!(sim.integrator.u, sim.integrator.p)
    return nothing
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
) = sim.integrator.p.α_sfc
Interfacer.get_field(sim::PrescribedIceSimulation, ::Val{:roughness_model}) = :constant
Interfacer.get_field(sim::PrescribedIceSimulation, ::Val{:surface_temperature}) =
    ice_surface_temperature.(
        sim.integrator.u.T_bulk,
        sim.integrator.p.params.T_base,
        sim.integrator.p.params.T_sfc_min,
        sim.integrator.p.params.T_freeze,
    )

"""
    ice_surface_temperature(T_bulk, T_base, T_sfc_min, T_freeze)

Approximates the surface temperature of the sea ice assuming
the ice temperature varies linearly between the ice surface and the base.
Clamps the result to be between `T_sfc_min` and `T_freeze`.
"""
ice_surface_temperature(T_bulk, T_base, T_sfc_min, T_freeze) =
    clamp(2 * T_bulk - T_base, T_sfc_min, T_freeze)

"""
    Interfacer.get_field(sim::PrescribedIceSimulation, ::Val{:energy})

Extension of Interfacer.get_field to get the energy of the sea-ice slab.
It multiplies the slab temperature by the heat capacity, density, and
hemispherically prescribed slab thickness.
"""
Interfacer.get_field(sim::PrescribedIceSimulation, ::Val{:energy}) =
    sim.integrator.p.params.ρ .* sim.integrator.p.params.c .* sim.integrator.u.T_bulk .*
    sim.integrator.p.h_ice

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

function Checkpointer.restart_model_state!(
    sim::PrescribedIceSimulation,
    input_file,
    comms_ctx,
)
    @invoke Checkpointer.restart_model_state!(
        sim::Interfacer.AbstractComponentSimulation,
        input_file,
        comms_ctx,
    )
    update_ice_albedo!(sim.integrator.u, sim.integrator.p)
    return nothing
end

###
### Sea ice model-specific functions (not explicitly required by ClimaCoupler.jl)
###

"""
    ice_energy_balance_tendency(
        T_bulk, area_fraction, F_turb_energy, SW_d, LW_d, α, h, params, dt,
    )

Pointwise slab sea-ice bulk temperature tendency [K/s] from the surface and
basal energy balance using cached surface albedo `α` and prescribed slab
thickness `h`, masked to zero in ice-free grid cells and bounded so `T_bulk`
does not exceed `(T_base + T_freeze) / 2` over `dt`, i.e. so the extrapolated
surface temperature `T_sfc = 2 T_bulk - T_base` does not exceed `T_freeze`.
Because ice thickness is prescribed, this cap discards surplus melt energy.
Note that F_conductive can be negative when T_base < T_sfc (during melt season).
"""
function ice_energy_balance_tendency(
    T_bulk::FT,
    area_fraction::FT,
    F_turb_energy::FT,
    SW_d::FT,
    LW_d::FT,
    α::FT,
    h::FT,
    params::IceSlabParameters{FT},
    dt::FT,
) where {FT}
    (; k_ice, T_base, ρ, c, ϵ, T_freeze, σ, T_sfc_min) = params
    T_sfc = ice_surface_temperature(T_bulk, T_base, T_sfc_min, T_freeze)
    F_conductive = (k_ice / h) * (T_base - T_sfc)
    rhs =
        (-F_turb_energy + (one(FT) - α) * SW_d + ϵ * (LW_d - σ * T_sfc^4) + F_conductive) /
        (h * ρ * c)
    # zero tendency where there is no ice
    rhs_masked = ifelse(area_fraction ≈ zero(FT), zero(FT), rhs)
    # cap where T_sfc reaches T_freeze
    T_bulk_max = (T_base + T_freeze) / 2
    return min(rhs_masked, (T_bulk_max - T_bulk) / dt)
end

"""
    ice_rhs!(dY, Y, p, t)

Rhs method in the form as required by `ClimeTimeSteppers`, with the tendency vector `dY`,
the state vector `Y`, the parameter vector `p`, and the simulation time `t` as input arguments.

This sea-ice energy formulation follows [Holloway and Manabe 1971](https://journals.ametsoc.org/view/journals/mwre/99/5/1520-0493_1971_099_0335_socbag_2_3_co_2.xml?tab_body=pdf),
where sea-ice concentrations and thickness are prescribed, and the model solves
for temperature (curbed at the freezing point, discarding surplus melt energy).
"""
function ice_rhs!(dY, Y, p, t)
    FT = eltype(Y)
    params = p.params

    # Update the cached area fraction with the current SIC
    if p.domain_type == "column"
        @. p.area_fraction = FT(1)
    else
        @. p.area_fraction = ifelse(isfinite(p.area_fraction), p.area_fraction, FT(0))
        # binary ice fraction (threshold at 0.5)
        if p.binary_area_fraction
            @. p.area_fraction = ifelse(p.area_fraction > FT(0.5), FT(1), FT(0))
        end

        # Overwrite ice fraction with the static land area fraction anywhere we have nonzero land area
        #  max needed to avoid Float32 errors (see issue #271; Heisenbug on HPC)
        @. p.area_fraction = max(min(p.area_fraction, FT(1) - p.land_fraction), FT(0))
    end

    dt_ft = FT(p.dt)
    @. dY.T_bulk = ice_energy_balance_tendency(
        Y.T_bulk,
        p.area_fraction,
        p.F_turb_energy,
        p.SW_d,
        p.LW_d,
        p.α_sfc,
        p.h_ice,
        params,
        dt_ft,
    )
    return nothing
end

"""
    dss_state!(sim::PrescribedIceSimulation)

Perform DSS on the state of a component simulation, intended to be used
before the initial step of a run. This method acts on prescribed ice simulations.
"""
dss_state!(sim::PrescribedIceSimulation) =
    Utilities.apply_dss!(sim.integrator.u, sim.integrator.p.dss_buffer)

function Checkpointer.get_model_cache(sim::PrescribedIceSimulation)
    return sim.integrator.p
end

function Checkpointer.restore_cache!(sim::PrescribedIceSimulation, new_cache)
    old_cache = Checkpointer.get_model_cache(sim)
    for p in propertynames(old_cache)
        if getproperty(old_cache, p) isa CC.Fields.Field
            ArrayType = ClimaComms.array_type(getproperty(old_cache, p))
            parent(getproperty(old_cache, p)) .=
                ArrayType(parent(getproperty(new_cache, p)))
        end
    end
end
