import SciMLBase
import ClimaCore as CC
import ClimaTimeSteppers as CTS
import ClimaCoupler: Checkpointer, FluxCalculator, Interfacer, Utilities, FieldExchanger
import Insolation.Parameters.InsolationParameters
import Dates
import ClimaAtmos as CA # for albedo calculation

###
### Functions required by ClimaCoupler.jl for a SurfaceModelSimulation
###
"""
    SlabOceanSimulation{P, I}

Equation:

    (h * ρ * c) dTdt = -(F_turb_energy + F_radiative)

"""
struct SlabOceanSimulation{P, I} <: Interfacer.OceanModelSimulation
    params::P
    integrator::I
end

# ocean parameters
Base.@kwdef struct OceanSlabParameters{FT <: AbstractFloat}
    h::FT = 20              # depth of the ocean [m]
    ρ::FT = 1500            # density of the ocean [kg / m3]
    c::FT = 800             # specific heat of the ocean [J / kg / K]
    T_init::FT = 271        # initial temperature of the ocean [K]
    z0m::FT = 5e-4          # roughness length for momentum [m]
    z0b::FT = 5e-4          # roughness length for heat [m]
    α::FT = 0.38            # albedo of the ocean [0, 1]
    evolving_switch::FT = 1 # switch to turn off the evolution of the ocean temperature [0 or 1]
end

"""
    slab_ocean_space_init(space, params)

Initialize the slab ocean prognostic variable (temperature), including an
anomaly in the tropics by default.
"""
function slab_ocean_space_init(space, params)
    FT = CC.Spaces.undertype(space)
    coords = CC.Fields.coordinate_field(space)

    # initial condition
    # FT(271) close to the average of T_1 in atmos
    T_sfc = params.T_init .+ temp_anomaly.(coords)

    # prognostic variable
    Y = CC.Fields.FieldVector(; T_sfc = T_sfc)

    return Y, space
end

"""
    SlabOceanSimulation(::Type{FT}; tspan, dt, saveat, space, area_fraction, stepper = CTS.RK4()) where {FT}

Initializes the `DiffEq` problem, and creates a Simulation-type object containing the necessary information for `step!` in the coupling loop.
"""
function SlabOceanSimulation(
    ::Type{FT},
    start_date;
    tspan,
    dt,
    saveat,
    space,
    area_fraction,
    thermo_params,
    stepper = CTS.RK4(),
    evolving = true,
) where {FT}

    evolving_switch = evolving ? FT(1) : FT(0)
    params = OceanSlabParameters{FT}(evolving_switch = evolving_switch)

    Y, space = slab_ocean_space_init(space, params)
    cache = (
        params = params,
        F_turb_energy = CC.Fields.zeros(space),
        F_radiative = CC.Fields.zeros(space),
        q_sfc = CC.Fields.zeros(space),
        ρ_sfc = CC.Fields.zeros(space),
        area_fraction = area_fraction,
        thermo_params = thermo_params,
        α_direct = CC.Fields.ones(space) .* params.α,
        α_diffuse = CC.Fields.ones(space) .* params.α,
        u_atmos = CC.Fields.zeros(space),
        v_atmos = CC.Fields.zeros(space),
        start_date = start_date,
        # add dss_buffer to cache to avoid runtime dss allocation
        dss_buffer = CC.Spaces.create_dss_buffer(Y),
    )

    ode_algo = CTS.ExplicitAlgorithm(stepper)
    ode_function =
        CTS.ClimaODEFunction(; T_exp! = slab_ocean_rhs!, dss! = (Y, p, t) -> CC.Spaces.weighted_dss!(Y, p.dss_buffer))
    if typeof(dt) isa Number
        dt = Float64(dt)
        tspan = Float64.(tspan)
        saveat = Float64.(saveat)
    end
    problem = SciMLBase.ODEProblem(ode_function, Y, tspan, cache)
    integrator = SciMLBase.init(problem, ode_algo, dt = dt, saveat = saveat, adaptive = false)

    sim = SlabOceanSimulation(params, integrator)

    # DSS state to ensure we have continuous fields
    dss_state!(sim)
    return sim
end

# extensions required by Interfacer
Interfacer.get_field(sim::SlabOceanSimulation, ::Val{:area_fraction}) = sim.integrator.p.area_fraction
Interfacer.get_field(sim::SlabOceanSimulation, ::Val{:roughness_buoyancy}) = sim.integrator.p.params.z0b
Interfacer.get_field(sim::SlabOceanSimulation, ::Val{:roughness_momentum}) = sim.integrator.p.params.z0m
Interfacer.get_field(sim::SlabOceanSimulation, ::Val{:surface_direct_albedo}) = sim.integrator.p.α_direct
Interfacer.get_field(sim::SlabOceanSimulation, ::Val{:surface_diffuse_albedo}) = sim.integrator.p.α_diffuse
Interfacer.get_field(sim::SlabOceanSimulation, ::Val{:surface_temperature}) = sim.integrator.u.T_sfc

"""
    Interfacer.get_field(sim::SlabOceanSimulation, ::Val{:energy})

Extension of Interfacer.get_field to get the energy of the ocean.
It multiplies the the slab temperature by the heat capacity, density, and depth.
"""
Interfacer.get_field(sim::SlabOceanSimulation, ::Val{:energy}) =
    sim.integrator.p.params.ρ .* sim.integrator.p.params.c .* sim.integrator.u.T_sfc .* sim.integrator.p.params.h

function Interfacer.update_field!(sim::SlabOceanSimulation, ::Val{:area_fraction}, field::CC.Fields.Field)
    sim.integrator.p.area_fraction .= field
end
function Interfacer.update_field!(sim::SlabOceanSimulation, ::Val{:air_density}, field)
    Interfacer.remap!(sim.integrator.p.ρ_sfc, field)
end
function Interfacer.update_field!(sim::SlabOceanSimulation, ::Val{:u_atmos}, field::CC.Fields.Field)
    Interfacer.remap!(sim.integrator.p.u_atmos, field)
end
function Interfacer.update_field!(sim::SlabOceanSimulation, ::Val{:v_atmos}, field::CC.Fields.Field)
    Interfacer.remap!(sim.integrator.p.v_atmos, field)
end
function Interfacer.update_field!(sim::SlabOceanSimulation, ::Val{:radiative_energy_flux_sfc}, field)
    Interfacer.remap!(sim.integrator.p.F_radiative, field)
end
function Interfacer.update_field!(sim::SlabOceanSimulation, ::Val{:turbulent_energy_flux}, field)
    Interfacer.remap!(sim.integrator.p.F_turb_energy, field)
end
function Interfacer.update_field!(sim::SlabOceanSimulation, ::Val{:surface_direct_albedo}, field::CC.Fields.Field)
    Interfacer.remap!(sim.integrator.p.α_direct, field)
end
function Interfacer.update_field!(sim::SlabOceanSimulation, ::Val{:surface_diffuse_albedo}, field::CC.Fields.Field)
    Interfacer.remap!(sim.integrator.p.α_diffuse, field)
end

# extensions required by FieldExchanger
Interfacer.step!(sim::SlabOceanSimulation, t) = Interfacer.step!(sim.integrator, t - sim.integrator.t, true)

"""
Extend Interfacer.add_coupler_fields! to add the fields required for SlabOceanSimulation.

The fields added are:
- `:ρ_sfc` (for humidity calculation)
- `:F_radiative` (for radiation input)
- `:u_int` (for water albedo calculation)
- `:v_int` (for water albedo calculation)
"""
function Interfacer.add_coupler_fields!(coupler_field_names, ::SlabOceanSimulation)
    ocean_coupler_fields = [:ρ_sfc, :F_radiative, :u_int, :v_int]
    push!(coupler_field_names, ocean_coupler_fields...)
end

function FluxCalculator.update_turbulent_fluxes!(sim::SlabOceanSimulation, fields::NamedTuple)
    (; F_lh, F_sh) = fields
    @. sim.integrator.p.F_turb_energy = F_lh + F_sh
end

"""
    Checkpointer.get_model_prog_state(sim::SlabOceanSimulation)

Extension of Checkpointer.get_model_prog_state to get the model state.
"""
function Checkpointer.get_model_prog_state(sim::SlabOceanSimulation)
    return sim.integrator.u
end

###
### Slab ocean model-specific functions (not explicitly required by ClimaCoupler.jl)
###
# ode
function slab_ocean_rhs!(dY, Y, cache, t)
    p, F_turb_energy, F_radiative = cache
    rhs = @. -(F_turb_energy + F_radiative) / (p.h * p.ρ * p.c)

    # Note that the area fraction has already been applied to the fluxes,
    #  so we don't need to multiply by it here.
    @. dY.T_sfc = rhs * p.evolving_switch
    @. cache.q_sfc = TD.q_vap_saturation_generic.(cache.thermo_params, Y.T_sfc, cache.ρ_sfc, TD.Liquid())
end

"""
    temp_anomaly(coord)

Calculates a an anomaly to be added to the initial condition for temperature.
This default case includes only an anomaly at the tropics.
"""
function temp_anomaly(coord)
    # include tropics anomaly
    FT = eltype(coord)
    anom = FT(29 * exp(-coord.lat^2 / (2 * 26^2)))
    return anom
end

"""
    dss_state!(sim::SlabOceanSimulation)

Perform DSS on the state of a component simulation, intended to be used
before the initial step of a run. This method acts on slab ocean model sims.
"""
dss_state!(sim::SlabOceanSimulation) = CC.Spaces.weighted_dss!(sim.integrator.u, sim.integrator.p.dss_buffer)

"""
    FieldExchanger.update_sim!(::SlabOceanSimulation, csf, area_fraction)

Updates the air density (needed for the turbulent flux calculation)
and the direct and diffuse albedos of the ocean.
"""
function FieldExchanger.update_sim!(sim::SlabOceanSimulation, csf, area_fraction)
    update_field!(sim, Val(:air_density), csf.ρ_sfc)
    update_field!(sim, Val(:u_atmos), csf.u_int)
    update_field!(sim, Val(:v_atmos), csf.v_int)

    # Update the direct and diffuse albedos with the new atmospheric wind
    set_albedos!(sim, sim.integrator.t)
end

"""
    FieldExchanger.import_atmos_fields!(csf, sim::SlabOceanSimulation, atmos_sim)

Imports quantities from the coupled simulation fields into the ocean simulation.
This is similar to the default method defined in FieldExchanger.jl, but it also
includes atmospheric wind, which is required to compute the ocean albedo.
"""
function FieldExchanger.import_atmos_fields!(csf, sim::SlabOceanSimulation, atmos_sim)
    # surface density - needed for q_sat and requires atmos and sfc states, so it is calculated and saved in the coupler
    Interfacer.remap!(csf.ρ_sfc, FluxCalculator.calculate_surface_air_density(atmos_sim, csf.T_sfc)) # TODO: generalize for PartitionedStateFluxes (#445) (use individual T_sfc)

    # radiative fluxes
    Interfacer.get_field(csf.F_radiative, atmos_sim, Val(:radiative_energy_flux_sfc))

    # precipitation
    Interfacer.get_field(csf.P_liq, atmos_sim, Val(:liquid_precipitation))
    Interfacer.get_field(csf.P_snow, atmos_sim, Val(:snow_precipitation))

    # wind
    Interfacer.get_field(csf.u_int, atmos_sim, Val(:u_int))
    Interfacer.get_field(csf.v_int, atmos_sim, Val(:v_int))
end

"""
    function set_albedos!(sim::SlabOceanSimulation, t)

Set the direct and diffuse albedos of the ocean based on the current date and
the atmospheric wind. The albedos are calculated using the `surface_albedo_direct`
and `surface_albedo_diffuse` functions from the `ClimaAtmos` module, so this
is dependent on running with `ClimaAtmosSimulation` as the atmosphere simulation.
"""
function set_albedos!(sim::SlabOceanSimulation, t)
    u = sim.integrator.u
    p = sim.integrator.p
    FT = eltype(u)

    # Compute the current date
    current_date = t isa ITime ? date(t) : p.start_date + Dates.second(t)

    # TODO: Where does this date0 come from?
    date0 = DateTime("2000-01-01T11:58:56.816")
    insolation_params = InsolationParameters(FT)
    d, δ, η_UTC = FT.(Insolation.helper_instantaneous_zenith_angle(current_date, date0, insolation_params))

    # Get the atmospheric wind vector and the cosine of the zenith angle
    surface_coords = Fields.coordinate_field(CC.Spaces.level(u.T_sfc, CC.Spaces.nlevels(u.T_sfc)))
    (zenith_angle, _, _) = @. instantaneous_zenith_angle(d, δ, η_UTC, surface_coords.long, surface_coords.lat) # the tuple is (zenith angle, azimuthal angle, earth-sun distance)
    wind_atmos = SA.SVector{2, FT}(p.u_atmos, p.v_atmos) # wind vector from components
    λ = FT(0) # spectral wavelength (not used for now)
    max_zenith_angle = FT(π) / 2 - eps(FT)
    cos_zenith = CC.Fields.array2field(cos(min(zenith_angle, max_zenith_angle)), axes(u)) # cosine of the zenith angle, as a ClimaCore field

    # Use the albedo model from ClimaAtmos
    α_model = CA.RegressionFunctionAlbedo{FT}()
    update_field!(sim, Val(:surface_direct_albedo), CA.surface_albedo_direct(α_model).(λ, cos_zenith, wind_atmos))
    update_field!(sim, Val(:surface_diffuse_albedo), CA.surface_albedo_diffuse(α_model).(λ, cos_zenith, wind_atmos))

    return nothing
end
