import SciMLBase
import ClimaCore as CC
import ClimaTimeSteppers as CTS
import ClimaUtilities.TimeVaryingInputs: TimeVaryingInput, evaluate!
import ClimaUtilities.ClimaArtifacts: @clima_artifact
import Interpolations # triggers InterpolationsExt in ClimaUtilities
import Thermodynamics as TD
import Statistics: mean
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
    α = FT(0.7),
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
    Y = CC.Fields.FieldVector(T_bulk = ones(space) .* params.T_freeze .- FT(15.0))
    return Y
end

"""
    slab_ice_space_init_from_file(::Type{FT}, space, params, sic_data, start_date, t_start, ice_fraction) where {FT}

Initialize ice temperature from ERA5 ice layer temperatures (ISTL1-4) if available.
Uses the average of the 4 IFS sea-ice slab layers:
- Layer 1: 0-7cm
- Layer 2: 7-28cm  
- Layer 3: 28-100cm
- Layer 4: 100-150cm

Only uses ERA5 temperatures where there is ice (ice_fraction > 0). In ice-free regions,
uses the fallback temperature to avoid initializing with values at/above freezing which
can cause numerical issues.

Falls back to constant initialization if the data is not available or not finite.

Note: Uses TimeVaryingInput because the ISTL data has a time dimension in the NetCDF file.
We evaluate at t_start to get the initial values.
"""
function slab_ice_space_init_from_file(
    ::Type{FT},
    space,
    params,
    sic_data,
    start_date,
    t_start,
    ice_fraction,
) where {FT}
    T_fallback = params.T_freeze - FT(15.0)

    T_bulk_init = try
        # Read the 4 ice temperature layers from ERA5 using TimeVaryingInput
        # TimeVaryingInput is required because the ISTL data has a time dimension
        @info "PrescribedIce: attempting to read ice temperature from ERA5 layers (ISTL1-4)"
        @info "PrescribedIce: reference_date=$start_date, t_start=$t_start"

        ISTL1_input = TimeVaryingInput(sic_data, "ISTL1", space; reference_date = start_date)
        ISTL2_input = TimeVaryingInput(sic_data, "ISTL2", space; reference_date = start_date)
        ISTL3_input = TimeVaryingInput(sic_data, "ISTL3", space; reference_date = start_date)
        ISTL4_input = TimeVaryingInput(sic_data, "ISTL4", space; reference_date = start_date)

        ISTL1 = CC.Fields.zeros(space)
        ISTL2 = CC.Fields.zeros(space)
        ISTL3 = CC.Fields.zeros(space)
        ISTL4 = CC.Fields.zeros(space)

        evaluate!(ISTL1, ISTL1_input, t_start)
        evaluate!(ISTL2, ISTL2_input, t_start)
        evaluate!(ISTL3, ISTL3_input, t_start)
        evaluate!(ISTL4, ISTL4_input, t_start)

        # Debug: log individual layer statistics
        for (name, field) in [("ISTL1", ISTL1), ("ISTL2", ISTL2), ("ISTL3", ISTL3), ("ISTL4", ISTL4)]
            vals = parent(field)
            n_nan = sum(isnan.(vals))
            n_inf = sum(isinf.(vals))
            @info "PrescribedIce: $name stats" minimum(vals) maximum(vals) mean(vals) n_nan n_inf
        end

        # Average the 4 layers for bulk temperature
        T_avg = @. (ISTL1 + ISTL2 + ISTL3 + ISTL4) / 4

        # Debug: log statistics of raw ERA5 ice temperatures
        T_avg_vals = parent(T_avg)
        n_nan_raw = sum(isnan.(T_avg_vals))
        n_inf_raw = sum(isinf.(T_avg_vals))
        @info "PrescribedIce: ERA5 ice temp stats (raw)" minimum(T_avg_vals) maximum(T_avg_vals) mean(T_avg_vals) n_nan_raw n_inf_raw

        # Only use ERA5 temperatures where there is ice; use fallback elsewhere
        # This prevents using warm (~273K) values from ice-free ocean regions
        @. T_avg = ifelse(ice_fraction > 0, T_avg, T_fallback)

        # Fall back to constant where data is not finite
        @. T_avg = ifelse(isfinite(T_avg), T_avg, T_fallback)

        # Clamp temperature to be at most T_freeze to prevent numerical issues
        # (ice temperature should never exceed the freezing point)
        @. T_avg = min(T_avg, params.T_freeze)

        # Debug: log statistics after processing
        T_avg_vals_processed = parent(T_avg)
        n_nan_processed = sum(isnan.(T_avg_vals_processed))
        @info "PrescribedIce: ice temp stats (after masking & clamping)" minimum(T_avg_vals_processed) maximum(T_avg_vals_processed) mean(T_avg_vals_processed) n_nan_processed

        @info "PrescribedIce: initialized ice temperature from ERA5 data (average of ISTL1-4)"
        T_avg
    catch e
        @warn "PrescribedIce: could not read ice temperature from file, using constant initialization" exception =
            (e, catch_backtrace())
        ones(space) .* T_fallback
    end

    # Final check for NaNs
    if any(isnan, parent(T_bulk_init))
        @warn "PrescribedIce: NaNs detected in ice temperature after initialization, replacing with fallback"
        @. T_bulk_init = ifelse(isnan(T_bulk_init), T_fallback, T_bulk_init)
    end

    Y = CC.Fields.FieldVector(T_bulk = T_bulk_init)
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
    binary_area_fraction::Bool = true,
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

    # Ensure ice fraction is finite and not NaN
    @. ice_fraction = ifelse(isfinite(ice_fraction), ice_fraction, FT(0))
    # binary ice fraction (threshold at 0.5)
    if binary_area_fraction
        @. ice_fraction = ifelse(ice_fraction > FT(0.5), FT(1), FT(0))
    end
    # max/min needed to avoid Float32 errors (see issue #271; Heisenbug on HPC)
    @. ice_fraction = max(min(ice_fraction, FT(1) - land_fraction), FT(0))

    params = IceSlabParameters{FT}(coupled_param_dict)

    # Initialize ice temperature: use ERA5 ice layer data if available, otherwise constant
    Y = if !isnothing(sic_path)
        # Custom SIC file provided - try to read ERA5 ice temperatures
        # Pass ice_fraction so we only use ERA5 temps where there's actually ice
        slab_ice_space_init_from_file(
            FT,
            space,
            params,
            sic_data,
            start_date,
            tspan[1],
            ice_fraction,
        )
    else
        # Using default artifact - use constant initialization
        slab_ice_space_init(FT, space, params)
    end

    # Post-initialization verification
    T_bulk_vals = parent(Y.T_bulk)
    @info "PrescribedIce: POST-INIT Y.T_bulk verification" minimum(T_bulk_vals) maximum(T_bulk_vals) mean(T_bulk_vals) n_nan=sum(isnan.(T_bulk_vals)) n_inf=sum(isinf.(T_bulk_vals))
    
    # Extra safety: ensure no NaNs in Y after initialization
    if any(isnan, T_bulk_vals)
        T_fallback = params.T_freeze - FT(15.0)
        @warn "PrescribedIce: POST-INIT found NaNs in Y.T_bulk, replacing with fallback=$T_fallback"
        @. Y.T_bulk = ifelse(isnan(Y.T_bulk), T_fallback, Y.T_bulk)
    end
    
    # Check ice_fraction consistency with Y.T_bulk initialization
    ice_vals = parent(ice_fraction)
    @info "PrescribedIce: ice_fraction at init" minimum(ice_vals) maximum(ice_vals) n_ice_points=sum(ice_vals .> 0)
    cache = (;
        F_turb_energy = CC.Fields.zeros(space),
        SW_d = CC.Fields.zeros(space),
        LW_d = CC.Fields.zeros(space),
        area_fraction = ice_fraction,
        SIC_timevaryinginput = SIC_timevaryinginput,
        land_fraction = land_fraction,
        binary_area_fraction = binary_area_fraction,
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
    if dt isa Number
        dt = Float64(dt)
        tspan = Float64.(tspan)
        saveat = Float64.(saveat)
    end
    @info "PrescribedIce: dt type and value" typeof(dt) dt
    problem = SciMLBase.ODEProblem(ode_function, Y, tspan, (; cache..., params = params))
    integrator =
        SciMLBase.init(problem, ode_algo, dt = dt, saveat = saveat, adaptive = false)

    sim = PrescribedIceSimulation(params, integrator)

    # DSS state to ensure we have continuous fields
    dss_state!(sim)
    
    # Post-DSS verification
    T_bulk_post_dss = parent(sim.integrator.u.T_bulk)
    @info "PrescribedIce: POST-DSS Y.T_bulk verification" minimum(T_bulk_post_dss) maximum(T_bulk_post_dss) mean(T_bulk_post_dss) n_nan=sum(isnan.(T_bulk_post_dss))
    if any(isnan, T_bulk_post_dss)
        @error "PrescribedIce: DSS introduced NaNs into Y.T_bulk!"
    end
    
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
function Interfacer.get_field(sim::PrescribedIceSimulation, ::Val{:surface_temperature})
    FT = eltype(sim.integrator.u)
    T_sfc = ice_surface_temperature.(sim.integrator.u.T_bulk, sim.integrator.p.params.T_base)
    
    # Clamp T_sfc to a minimum of 200K to prevent numerical issues in atmosphere
    # Very cold surface temps (< 200K) can cause issues in radiation calculations
    T_sfc_min_allowed = FT(210)
    T_sfc_vals = parent(T_sfc)
    n_too_cold = sum(T_sfc_vals .< T_sfc_min_allowed)
    if n_too_cold > 0
        @warn "PrescribedIce: clamping $n_too_cold T_sfc values below $T_sfc_min_allowed K" T_sfc_min=minimum(T_sfc_vals)
        @. T_sfc = max(T_sfc, T_sfc_min_allowed)
    end
    
    T_sfc_vals = parent(T_sfc)
    if any(isnan, T_sfc_vals)
        T_bulk_vals = parent(sim.integrator.u.T_bulk)
        @warn "PrescribedIce: get_field(:surface_temperature) has NaNs" n_nan=sum(isnan.(T_sfc_vals)) T_sfc_min=minimum(T_sfc_vals) T_sfc_max=maximum(T_sfc_vals) T_bulk_min=minimum(T_bulk_vals) T_bulk_max=maximum(T_bulk_vals) T_base=sim.integrator.p.params.T_base
    end
    return T_sfc
end

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
    # Debug: check for NaNs in incoming field
    SW_vals = parent(sim.integrator.p.SW_d)
    if any(isnan, SW_vals)
        @warn "PrescribedIce: update_field!(:SW_d) received NaN values" n_nan=sum(isnan.(SW_vals))
    end
end
function Interfacer.update_field!(sim::PrescribedIceSimulation, ::Val{:LW_d}, field)
    Interfacer.remap!(sim.integrator.p.LW_d, field)
    # Debug: check for NaNs in incoming field
    LW_vals = parent(sim.integrator.p.LW_d)
    if any(isnan, LW_vals)
        @warn "PrescribedIce: update_field!(:LW_d) received NaN values" n_nan=sum(isnan.(LW_vals))
    end
end
function Interfacer.update_field!(
    sim::PrescribedIceSimulation,
    ::Val{:turbulent_energy_flux},
    field,
)
    Interfacer.remap!(sim.integrator.p.F_turb_energy, field)
    # Debug: check for NaNs in incoming field
    F_turb_vals = parent(sim.integrator.p.F_turb_energy)
    if any(isnan, F_turb_vals)
        @warn "PrescribedIce: update_field!(:turbulent_energy_flux) received NaN values" n_nan=sum(isnan.(F_turb_vals))
    end
end
Interfacer.update_field!(
    sim::PrescribedIceSimulation,
    ::Val{:turbulent_moisture_flux},
    field,
) = nothing

# extensions required by FieldExchanger
function Interfacer.step!(sim::PrescribedIceSimulation, t)
    FT = eltype(sim.integrator.u)
    T_freeze = sim.integrator.p.params.T_freeze
    
    # Check state before stepping
    T_bulk_pre = parent(sim.integrator.u.T_bulk)
    if any(isnan, T_bulk_pre)
        @error "PrescribedIce: NaN detected BEFORE step!" t current_t=sim.integrator.t n_nan=sum(isnan.(T_bulk_pre))
    end
    
    result = Interfacer.step!(sim.integrator, t - sim.integrator.t, true)
    
    # Clamp T_bulk to physically reasonable bounds after step:
    # - Cannot exceed T_freeze (ice would melt)
    # - Minimum of 200K (very cold but physically possible)
    T_bulk_min = FT(210)
    @. sim.integrator.u.T_bulk = clamp(sim.integrator.u.T_bulk, T_bulk_min, T_freeze)
    
    # Check state after stepping
    T_bulk_post = parent(sim.integrator.u.T_bulk)
    if any(isnan, T_bulk_post)
        @error "PrescribedIce: NaN detected AFTER step!" t n_nan=sum(isnan.(T_bulk_post)) minimum(T_bulk_post) maximum(T_bulk_post)
    end
    
    return result
end

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

# Global counter for debug logging (only log first few timesteps)
const ICE_RHS_DEBUG_COUNTER = Ref(0)
const ICE_RHS_DEBUG_MAX = 5  # Only log first N calls

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

    # Increment debug counter
    ICE_RHS_DEBUG_COUNTER[] += 1
    do_debug = ICE_RHS_DEBUG_COUNTER[] <= ICE_RHS_DEBUG_MAX

    # Helper function to check for NaNs and log
    function check_nans(field, name)
        vals = parent(field)
        n_nan = sum(isnan.(vals))
        n_inf = sum(isinf.(vals))
        if n_nan > 0 || n_inf > 0
            @warn "ice_rhs! NaN/Inf detected in $name" t n_nan n_inf minimum(vals) maximum(vals)
            return true
        end
        return false
    end

    # Debug: Check input state Y.T_bulk
    if do_debug || any(isnan, parent(Y.T_bulk))
        T_bulk_vals = parent(Y.T_bulk)
        n_nan = sum(isnan.(T_bulk_vals))
        n_inf = sum(isinf.(T_bulk_vals))
        @info "ice_rhs! [$(ICE_RHS_DEBUG_COUNTER[])] INPUT Y.T_bulk" t minimum(T_bulk_vals) maximum(T_bulk_vals) mean(T_bulk_vals) n_nan n_inf
        if n_nan > 0
            @error "ice_rhs! NaN in INPUT Y.T_bulk - this is the source!"
        end
    end

    # Debug: Check cache fields before computation
    if do_debug
        F_turb_vals = parent(p.F_turb_energy)
        SW_vals = parent(p.SW_d)
        LW_vals = parent(p.LW_d)
        @info "ice_rhs! [$(ICE_RHS_DEBUG_COUNTER[])] CACHE" F_turb_min=minimum(F_turb_vals) F_turb_max=maximum(F_turb_vals) SW_min=minimum(SW_vals) SW_max=maximum(SW_vals) LW_min=minimum(LW_vals) LW_max=maximum(LW_vals)
        check_nans(p.F_turb_energy, "F_turb_energy")
        check_nans(p.SW_d, "SW_d")
        check_nans(p.LW_d, "LW_d")
    end

    # Update the cached area fraction with the current SIC
    evaluate!(p.area_fraction, p.SIC_timevaryinginput, t)
    @. p.area_fraction = ifelse(isfinite(p.area_fraction), p.area_fraction, FT(0))
    # binary ice fraction (threshold at 0.5)
    if p.binary_area_fraction
        @. p.area_fraction = ifelse(p.area_fraction > FT(0.5), FT(1), FT(0))
    end

    # Overwrite ice fraction with the static land area fraction anywhere we have nonzero land area
    #  max needed to avoid Float32 errors (see issue #271; Heisenbug on HPC)
    @. p.area_fraction = max(min(p.area_fraction, FT(1) - p.land_fraction), FT(0))

    if do_debug
        af_vals = parent(p.area_fraction)
        @info "ice_rhs! [$(ICE_RHS_DEBUG_COUNTER[])] area_fraction" minimum(af_vals) maximum(af_vals) n_nonzero=sum(af_vals .> 0)
    end

    # Calculate the surface temperature
    # IMPORTANT: Clamp T_sfc to reasonable bounds (210K to T_freeze) INSIDE ice_rhs!
    # This is critical because the ODE solver uses intermediate stages where T_bulk
    # can go outside the post-step clamped range, leading to extremely cold T_sfc
    # values (e.g., T_bulk=200K → T_sfc=128.8K) that can cause issues.
    T_sfc_min = FT(210)  # Minimum physically reasonable ice surface temperature
    T_sfc = @. clamp(ice_surface_temperature(Y.T_bulk, T_base), T_sfc_min, T_freeze)
    
    if do_debug || any(isnan, parent(T_sfc))
        T_sfc_vals = parent(T_sfc)
        n_nan = sum(isnan.(T_sfc_vals))
        @info "ice_rhs! [$(ICE_RHS_DEBUG_COUNTER[])] T_sfc (2*T_bulk - T_base, clamped)" minimum(T_sfc_vals) maximum(T_sfc_vals) mean(T_sfc_vals) n_nan T_base
        if n_nan > 0
            @error "ice_rhs! NaN in T_sfc calculation"
        end
    end

    # Calculate the conductive flux, and set it to zero if the area fraction is zero
    F_conductive = @. k_ice / (h) * (T_base - T_sfc) # fluxes are defined to be positive when upward
    
    if do_debug || any(isnan, parent(F_conductive))
        Fc_vals = parent(F_conductive)
        n_nan = sum(isnan.(Fc_vals))
        @info "ice_rhs! [$(ICE_RHS_DEBUG_COUNTER[])] F_conductive" minimum(Fc_vals) maximum(Fc_vals) n_nan k_ice h
        if n_nan > 0
            @error "ice_rhs! NaN in F_conductive"
        end
    end

    # Compute LW emission term separately for debugging
    LW_emission = @. σ * T_sfc^4
    if do_debug || any(isnan, parent(LW_emission))
        LW_em_vals = parent(LW_emission)
        n_nan = sum(isnan.(LW_em_vals))
        @info "ice_rhs! [$(ICE_RHS_DEBUG_COUNTER[])] LW_emission (σ*T_sfc^4)" minimum(LW_em_vals) maximum(LW_em_vals) n_nan σ
        if n_nan > 0
            @error "ice_rhs! NaN in LW_emission - check if T_sfc has negative values!"
            # Check for negative T_sfc which would be unphysical
            T_sfc_vals = parent(T_sfc)
            n_negative = sum(T_sfc_vals .< 0)
            @error "ice_rhs! T_sfc negative count: $n_negative"
        end
    end

    rhs = @. (
        -p.F_turb_energy +
        (1 - α) * p.SW_d +
        ϵ * (p.LW_d - LW_emission) +
        F_conductive
    ) / (h * ρ * c)

    if do_debug || any(isnan, parent(rhs))
        rhs_vals = parent(rhs)
        n_nan = sum(isnan.(rhs_vals))
        @info "ice_rhs! [$(ICE_RHS_DEBUG_COUNTER[])] rhs (before masking)" minimum(rhs_vals) maximum(rhs_vals) n_nan
        if n_nan > 0
            @error "ice_rhs! NaN in rhs before masking"
        end
    end

    # Zero out tendencies where there is no ice, so that ice temperature remains constant there
    @. rhs = ifelse(p.area_fraction ≈ 0, zero(rhs), rhs)

    if do_debug || any(isnan, parent(rhs))
        rhs_vals = parent(rhs)
        n_nan = sum(isnan.(rhs_vals))
        @info "ice_rhs! [$(ICE_RHS_DEBUG_COUNTER[])] rhs (after masking)" minimum(rhs_vals) maximum(rhs_vals) n_nan
    end

    # Debug: Check the dt value and the limiter term
    # Convert dt to the same float type as the fields (FT) to avoid Float64/Float32 mismatch
    dt_val = FT(float(p.dt))
    limiter = @. (T_freeze - Y.T_bulk) / dt_val
    if do_debug || any(isnan, parent(limiter))
        lim_vals = parent(limiter)
        n_nan = sum(isnan.(lim_vals))
        @info "ice_rhs! [$(ICE_RHS_DEBUG_COUNTER[])] limiter ((T_freeze - T_bulk)/dt)" minimum(lim_vals) maximum(lim_vals) n_nan dt_val T_freeze
        if n_nan > 0
            @error "ice_rhs! NaN in limiter term"
        end
    end

    # If tendencies lead to temperature above freezing, set temperature to freezing
    @. dY.T_bulk = min(rhs, limiter)

    # Final check on output
    if do_debug || any(isnan, parent(dY.T_bulk))
        dY_vals = parent(dY.T_bulk)
        n_nan = sum(isnan.(dY_vals))
        @info "ice_rhs! [$(ICE_RHS_DEBUG_COUNTER[])] OUTPUT dY.T_bulk" minimum(dY_vals) maximum(dY_vals) n_nan
        if n_nan > 0
            @error "ice_rhs! NaN in OUTPUT dY.T_bulk!"
        end
    end
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
        if getproperty(old_cache, p) isa CC.Fields.Field
            ArrayType = ClimaComms.array_type(getproperty(old_cache, p))
            parent(getproperty(old_cache, p)) .=
                ArrayType(parent(getproperty(new_cache, p)))
        end
    end
end
