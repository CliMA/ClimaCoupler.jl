# atmos_init: for ClimaAtmos interface
import StaticArrays
import Statistics
import LinearAlgebra
import ClimaAtmos as CA
import ClimaAtmos: set_surface_albedo!
import ClimaAtmos.Parameters as CAP
import ClimaCore as CC
import ClimaCore.Geometry: ⊗
import SurfaceFluxes as SF
import Thermodynamics as TD
import ClimaCoupler: Checkpointer, FieldExchanger, FluxCalculator, Interfacer, Utilities

include("climaatmos_extra_diags.jl")

###
### Functions required by ClimaCoupler.jl for an AtmosModelSimulation
###
struct ClimaAtmosSimulation{P, Y, D, I} <: Interfacer.AtmosModelSimulation
    params::P
    Y_init::Y
    domain::D
    integrator::I
end
Interfacer.name(::ClimaAtmosSimulation) = "ClimaAtmosSimulation"

function hasradiation(integrator)
    return !isnothing(integrator.p.atmos.radiation_mode)
end

function hasmoisture(integrator)
    return !(integrator.p.atmos.moisture_model isa CA.DryModel)
end

function atmos_init(atmos_config)
    # By passing `parsed_args` to `AtmosConfig`, `parsed_args` overwrites the default atmos config
    FT = atmos_config.parsed_args["FLOAT_TYPE"] == "Float64" ? Float64 : Float32
    simulation = CA.get_simulation(atmos_config)
    (; integrator) = simulation
    Y = integrator.u
    center_space = axes(Y.c.ρe_tot)
    face_space = axes(Y.f.u₃)
    spaces = (; center_space = center_space, face_space = face_space)

    # define shorter references for long variable names to increase readability, and set to zero
    ρ_flux_h_tot = integrator.p.precomputed.sfc_conditions.ρ_flux_h_tot
    ρ_flux_uₕ = integrator.p.precomputed.sfc_conditions.ρ_flux_uₕ
    @. ρ_flux_h_tot = CC.Geometry.Covariant3Vector(FT(0.0))
    ρ_flux_uₕ.components .= Ref(StaticArrays.SMatrix{1, 2}([FT(0), FT(0)]))

    if hasmoisture(integrator)
        ρ_flux_q_tot = integrator.p.precomputed.sfc_conditions.ρ_flux_q_tot
        surface_rain_flux = integrator.p.precipitation.surface_rain_flux
        surface_snow_flux = integrator.p.precipitation.surface_snow_flux
        @. ρ_flux_q_tot = CC.Geometry.Covariant3Vector(FT(0.0))
        surface_rain_flux .= FT(0)
        surface_snow_flux .= FT(0)
    end
    if integrator.p.atmos.precip_model isa CA.Microphysics0Moment
        ᶜS_ρq_tot = integrator.p.precipitation.ᶜS_ρq_tot
        ᶜ3d_rain = integrator.p.precipitation.ᶜ3d_rain
        ᶜ3d_snow = integrator.p.precipitation.ᶜ3d_snow
        ᶜS_ρq_tot .= FT(0)
        ᶜ3d_rain .= FT(0)
        ᶜ3d_snow .= FT(0)
    end
    if hasradiation(integrator)
        ᶠradiation_flux = integrator.p.radiation.ᶠradiation_flux
        @. ᶠradiation_flux = CC.Geometry.WVector(FT(0))
    end

    sim = ClimaAtmosSimulation(integrator.p.params, Y, spaces, integrator)

    # DSS state to ensure we have continuous fields
    dss_state!(sim)
    return sim
end

# Extension of CA.set_surface_albedo! to set the surface albedo to 0.38 at the beginning of the simulation,
# so the initial callback initialization doesn't lead to NaNs in the radiation model.
# Subsequently, the surface albedo will be updated by the coupler, via water_albedo_from_atmosphere!.
function CA.set_surface_albedo!(Y, p, t, ::CA.CouplerAlbedo)
    if t == 0
        FT = eltype(Y)
        # set initial insolation initial conditions
        !(p.atmos.insolation isa CA.IdealizedInsolation) && CA.set_insolation_variables!(Y, p, t, p.atmos.insolation)
        # set surface albedo to 0.38
        @warn "Setting surface albedo to 0.38 at the beginning of the simulation"
        p.radiation.rrtmgp_model.direct_sw_surface_albedo .= FT(0.38)
        p.radiation.rrtmgp_model.diffuse_sw_surface_albedo .= FT(0.38)
    else
        nothing
    end
end

"""
    Checkpointer.get_model_prog_state(sim::ClimaAtmosSimulation)

Extension of Checkpointer.get_model_prog_state to get the model state.
"""
function Checkpointer.get_model_prog_state(sim::ClimaAtmosSimulation)
    return sim.integrator.u
end

"""
Interfacer.get_field(atmos_sim::ClimaAtmosSimulation, ::Val{:radiative_energy_flux_toa})

Extension of Interfacer.get_field to get the net TOA radiation, which is a sum of the
upward and downward longwave and shortwave radiation.
"""
function Interfacer.get_field(atmos_sim::ClimaAtmosSimulation, ::Val{:radiative_energy_flux_toa})
    FT = eltype(atmos_sim.integrator.u)

    if hasradiation(atmos_sim.integrator)
        face_space = axes(atmos_sim.integrator.u.f)
        nz_faces = length(CC.Spaces.vertical_topology(face_space).mesh.faces)

        (; face_lw_flux_dn, face_lw_flux_up, face_sw_flux_dn, face_sw_flux_up) =
            atmos_sim.integrator.p.radiation.rrtmgp_model

        LWd_TOA = CC.Fields.level(CC.Fields.array2field(FT.(face_lw_flux_dn), face_space), nz_faces - CC.Utilities.half)
        LWu_TOA = CC.Fields.level(CC.Fields.array2field(FT.(face_lw_flux_up), face_space), nz_faces - CC.Utilities.half)
        SWd_TOA = CC.Fields.level(CC.Fields.array2field(FT.(face_sw_flux_dn), face_space), nz_faces - CC.Utilities.half)
        SWu_TOA = CC.Fields.level(CC.Fields.array2field(FT.(face_sw_flux_up), face_space), nz_faces - CC.Utilities.half)

        return @. -(LWd_TOA + SWd_TOA - LWu_TOA - SWu_TOA)
    else
        return FT[0]
    end
end

function Interfacer.get_field(atmos_sim::ClimaAtmosSimulation, ::Val{:energy})
    integrator = atmos_sim.integrator
    p = integrator.p


    # return total energy and (if Microphysics0Moment) the energy lost due to precipitation removal
    if p.atmos.precip_model isa CA.Microphysics0Moment
        ᶜts = p.precomputed.ᶜts
        ᶜΦ = p.core.ᶜΦ
        ᶜS_ρq_tot = p.precipitation.ᶜS_ρq_tot
        thermo_params = get_thermo_params(atmos_sim)
        return integrator.u.c.ρe_tot .-
               ᶜS_ρq_tot .* CA.e_tot_0M_precipitation_sources_helper.(Ref(thermo_params), ᶜts, ᶜΦ) .* integrator.dt
    else
        return integrator.u.c.ρe_tot
    end
end

# helpers for get_field extensions, dipatchable on different moisture model options and radiation modes

surface_rain_flux(::CA.DryModel, integrator) = StaticArrays.SVector(eltype(integrator.u)(0))
surface_rain_flux(::Union{CA.EquilMoistModel, CA.NonEquilMoistModel}, integrator) =
    integrator.p.precipitation.surface_rain_flux

surface_snow_flux(::CA.DryModel, integrator) = StaticArrays.SVector(eltype(integrator.u)(0))
surface_snow_flux(::Union{CA.EquilMoistModel, CA.NonEquilMoistModel}, integrator) =
    integrator.p.precipitation.surface_snow_flux

surface_radiation_flux(::Nothing, integrator) = StaticArrays.SVector(eltype(integrator.u)(0))
surface_radiation_flux(::CA.RRTMGPI.AbstractRRTMGPMode, integrator) =
    CC.Fields.level(integrator.p.radiation.ᶠradiation_flux, CC.Utilities.half)

moisture_flux(::CA.DryModel, integrator) = StaticArrays.SVector(eltype(integrator.u)(0))
moisture_flux(::Union{CA.EquilMoistModel, CA.NonEquilMoistModel}, integrator) =
    CC.Geometry.WVector.(integrator.p.precomputed.sfc_conditions.ρ_flux_q_tot)

ρq_tot(::CA.DryModel, integrator) = StaticArrays.SVector(eltype(integrator.u)(0))
ρq_tot(::Union{CA.EquilMoistModel, CA.NonEquilMoistModel}, integrator) = integrator.u.c.ρq_tot

# extensions required by the Interfacer
Interfacer.get_field(sim::ClimaAtmosSimulation, ::Val{:air_density}) =
    TD.air_density.(thermo_params, sim.integrator.p.precomputed.ᶜts)
Interfacer.get_field(sim::ClimaAtmosSimulation, ::Val{:air_temperature}) =
    TD.air_temperature.(thermo_params, sim.integrator.p.precomputed.ᶜts)
Interfacer.get_field(sim::ClimaAtmosSimulation, ::Val{:liquid_precipitation}) =
    surface_rain_flux(sim.integrator.p.atmos.moisture_model, sim.integrator)
Interfacer.get_field(sim::ClimaAtmosSimulation, ::Val{:radiative_energy_flux_sfc}) =
    surface_radiation_flux(sim.integrator.p.atmos.radiation_mode, sim.integrator)
Interfacer.get_field(sim::ClimaAtmosSimulation, ::Val{:snow_precipitation}) =
    surface_snow_flux(sim.integrator.p.atmos.moisture_model, sim.integrator)
Interfacer.get_field(sim::ClimaAtmosSimulation, ::Val{:turbulent_energy_flux}) =
    CC.Geometry.WVector.(sim.integrator.p.precomputed.sfc_conditions.ρ_flux_h_tot)
Interfacer.get_field(sim::ClimaAtmosSimulation, ::Val{:turbulent_moisture_flux}) =
    moisture_flux(sim.integrator.p.atmos.moisture_model, sim.integrator)
Interfacer.get_field(sim::ClimaAtmosSimulation, ::Val{:thermo_state_int}) =
    CC.Spaces.level(sim.integrator.p.precomputed.ᶜts, 1)
Interfacer.get_field(atmos_sim::ClimaAtmosSimulation, ::Val{:water}) =
    ρq_tot(atmos_sim.integrator.p.atmos.moisture_model, atmos_sim.integrator)
function Interfacer.update_field!(sim::ClimaAtmosSimulation, ::Val{:surface_temperature}, csf)
    # note that this field is also being updated internally by the surface thermo state in ClimaAtmos
    # if turbulent fluxes are calculated, to ensure consistency. In case the turbulent fluxes are not
    # calculated, we update the field here.
    sim.integrator.p.radiation.rrtmgp_model.surface_temperature .= CC.Fields.field2array(csf.T_S)
end
# extensions required by FluxCalculator (partitioned fluxes)
Interfacer.get_field(sim::ClimaAtmosSimulation, ::Val{:height_int}) =
    CC.Spaces.level(CC.Fields.coordinate_field(sim.integrator.u.c).z, 1)
Interfacer.get_field(sim::ClimaAtmosSimulation, ::Val{:height_sfc}) =
    CC.Spaces.level(CC.Fields.coordinate_field(sim.integrator.u.f).z, CC.Utilities.half)
function Interfacer.get_field(sim::ClimaAtmosSimulation, ::Val{:uv_int})
    uₕ_int = CC.Geometry.UVVector.(CC.Spaces.level(sim.integrator.u.c.uₕ, 1))
    return @. StaticArrays.SVector(uₕ_int.components.data.:1, uₕ_int.components.data.:2)
end

function Interfacer.update_field!(atmos_sim::ClimaAtmosSimulation, ::Val{:co2}, field)
    if atmos_sim.integrator.p.atmos.radiation_mode isa CA.RRTMGPI.GrayRadiation
        @warn("Gray radiation model initialized, skipping CO2 update", maxlog = 1)
        return
    else
        atmos_sim.integrator.p.radiation.rrtmgp_model.volume_mixing_ratio_co2 .= Statistics.mean(parent(field))
    end
end
# extensions required by the Interfacer
function Interfacer.update_field!(sim::ClimaAtmosSimulation, ::Val{:surface_direct_albedo}, field)
    sim.integrator.p.radiation.rrtmgp_model.direct_sw_surface_albedo .=
        reshape(CC.Fields.field2array(field), 1, length(parent(field)))
end

function Interfacer.update_field!(sim::ClimaAtmosSimulation, ::Val{:surface_diffuse_albedo}, field)
    sim.integrator.p.radiation.rrtmgp_model.diffuse_sw_surface_albedo .=
        reshape(CC.Fields.field2array(field), 1, length(parent(field)))
end

function Interfacer.update_field!(sim::ClimaAtmosSimulation, ::Val{:turbulent_fluxes}, fields)
    (; F_turb_energy, F_turb_moisture, F_turb_ρτxz, F_turb_ρτyz) = fields

    Y = sim.integrator.u
    surface_local_geometry = CC.Fields.level(CC.Fields.local_geometry_field(Y.f), CC.Utilities.half)
    surface_normal = @. CA.C3(CA.unit_basis_vector_data(CA.C3, surface_local_geometry))

    # get template objects for the contravariant components of the momentum fluxes (required by Atmos boundary conditions)
    vec_ct12_ct1 = @. CA.CT12(CA.CT2(CA.unit_basis_vector_data(CA.CT1, surface_local_geometry)), surface_local_geometry)
    vec_ct12_ct2 = @. CA.CT12(CA.CT2(CA.unit_basis_vector_data(CA.CT2, surface_local_geometry)), surface_local_geometry)

    sim.integrator.p.precomputed.sfc_conditions.ρ_flux_uₕ .= (
        surface_normal .⊗
        CA.C12.(
            Utilities.swap_space!(axes(vec_ct12_ct1), F_turb_ρτxz) .* vec_ct12_ct1 .+
            Utilities.swap_space!(axes(vec_ct12_ct2), F_turb_ρτyz) .* vec_ct12_ct2,
            surface_local_geometry,
        )
    )

    parent(sim.integrator.p.precomputed.sfc_conditions.ρ_flux_h_tot) .= parent(F_turb_energy) .* parent(surface_normal) # (shf + lhf)
    parent(sim.integrator.p.precomputed.sfc_conditions.ρ_flux_q_tot) .=
        parent(F_turb_moisture) .* parent(surface_normal) # (evap)

    # TODO: see if Atmos can rever to a simpler solution
end

# extensions required by FieldExchanger
Interfacer.step!(sim::ClimaAtmosSimulation, t) = Interfacer.step!(sim.integrator, t - sim.integrator.t, true)
Interfacer.reinit!(sim::ClimaAtmosSimulation) = Interfacer.reinit!(sim.integrator)

function FieldExchanger.update_sim!(atmos_sim::ClimaAtmosSimulation, csf, turbulent_fluxes)

    u = atmos_sim.integrator.u
    p = atmos_sim.integrator.p
    t = atmos_sim.integrator.t

    !isempty(atmos_sim.integrator.p.radiation) &&
        !(p.atmos.insolation isa CA.IdealizedInsolation) &&
        CA.set_insolation_variables!(u, p, t, p.atmos.insolation)

    if hasradiation(atmos_sim.integrator)
        Interfacer.update_field!(atmos_sim, Val(:surface_direct_albedo), csf.surface_direct_albedo)
        Interfacer.update_field!(atmos_sim, Val(:surface_diffuse_albedo), csf.surface_diffuse_albedo)
    end

    !isempty(atmos_sim.integrator.p.radiation) && Interfacer.update_field!(atmos_sim, Val(:surface_temperature), csf)

    if turbulent_fluxes isa FluxCalculator.PartitionedStateFluxes
        Interfacer.update_field!(atmos_sim, Val(:turbulent_fluxes), csf)
    end
end

"""
    FluxCalculator.calculate_surface_air_density(atmos_sim::ClimaAtmosSimulation, T_S::CC.Fields.Field)

Extension for this function to calculate surface density.
"""
function FluxCalculator.calculate_surface_air_density(atmos_sim::ClimaAtmosSimulation, T_S::CC.Fields.Field)
    thermo_params = get_thermo_params(atmos_sim)
    ts_int = Interfacer.get_field(atmos_sim, Val(:thermo_state_int))
    FluxCalculator.extrapolate_ρ_to_sfc.(Ref(thermo_params), ts_int, Utilities.swap_space!(axes(ts_int.ρ), T_S))
end

# FluxCalculator.get_surface_params required by FluxCalculator (partitioned fluxes)
FluxCalculator.get_surface_params(sim::ClimaAtmosSimulation) = CAP.surface_fluxes_params(sim.integrator.p.params)

###
### ClimaAtmos.jl model-specific functions (not explicitly required by ClimaCoupler.jl)
###
"""
    get_atmos_config_dict(coupler_dict::Dict, job_id::String)

Returns the specified atmospheric configuration (`atmos_config`) overwitten by arguments
in the coupler dictionary (`config_dict`).
The returned `atmos_config` dictionary will then be used to set up the atmosphere simulation.

The `atmos_config_repo` flag allows us to use a configuration specified within
either the ClimaCoupler or ClimaAtmos repos, which is useful for direct
coupled/atmos-only comparisons.

In this function, parameters are overwritten in a specific order, from lowest to highest priority:
    1. Default atmos config
    2. Provided atmos config file (if any)
    3. TOML parameter file
    5. Output directory and timestep are set explicitly based on the coupler config

The TOML parameter file to use is chosen using the following priority:
If a coupler TOML file is provided, it is used. Otherwise we use an atmos TOML
file if it's provided. If neither is provided, we use a default coupler TOML file.
"""
function get_atmos_config_dict(coupler_dict::Dict, job_id::String)
    atmos_config_file = coupler_dict["atmos_config_file"]
    atmos_config_repo = coupler_dict["atmos_config_repo"]
    # override default or specified configs with coupler arguments, and set the correct atmos config_file
    if atmos_config_repo == "ClimaCoupler"
        @assert !isnothing(atmos_config_file) "Must specify `atmos_config_file` within ClimaCoupler."

        @info "Using Atmos configuration from ClimaCoupler in $atmos_config_file"
        atmos_config = merge(
            CA.override_default_config(joinpath(pkgdir(ClimaCoupler), atmos_config_file)),
            coupler_dict,
            Dict("config_file" => atmos_config_file),
        )
    elseif atmos_config_repo == "ClimaAtmos"
        if isnothing(atmos_config_file)
            @info "Using Atmos default configuration"
            atmos_config = merge(CA.default_config_dict(), coupler_dict, Dict("config_file" => atmos_config_file))
        else
            @info "Using Atmos configuration from $atmos_config_file"
            atmos_config = merge(
                CA.override_default_config(joinpath(pkgdir(CA), atmos_config_file)),
                coupler_dict,
                Dict("config_file" => atmos_config_file),
            )
        end
    else
        error("Invalid `atmos_config_repo`; please use \"ClimaCoupler\" or \"ClimaAtmos\"")
    end

    # use coupler toml if atmos is not defined
    atmos_toml_file = atmos_config["toml"]
    coupler_toml_file = coupler_dict["coupler_toml_file"]
    default_toml_file = "toml/default_coarse.toml"

    toml_file = !isempty(atmos_toml_file) ? joinpath(pkgdir(CA), atmos_toml_file[1]) : nothing
    toml_file = !isnothing(coupler_toml_file) ? joinpath(pkgdir(ClimaCoupler), coupler_toml_file) : toml_file
    toml_file = isnothing(toml_file) ? joinpath(pkgdir(ClimaCoupler), default_toml_file) : toml_file

    if !isnothing(toml_file)
        @info "Overwriting Atmos parameters from $toml_file"
        atmos_config = merge(atmos_config, Dict("toml" => [toml_file]))
    end

    # Specify atmos output directory to be inside the coupler output directory
    atmos_output_dir = joinpath(coupler_dict["coupler_output_dir"], job_id, "clima_atmos")
    atmos_config["output_dir"] = atmos_output_dir

    # Access extra atmosphere diagnostics from coupler so we can rename for atmos code
    atmos_config["diagnostics"] = coupler_dict["extra_atmos_diagnostics"]

    # The Atmos `get_simulation` function expects the atmos config to contains its timestep size
    # in the `dt` field. If there is a `dt_atmos` field in coupler_dict, we add it to the atmos config as `dt`
    dt_atmos = haskey(coupler_dict, "dt_atmos") ? coupler_dict["dt_atmos"] : coupler_dict["dt"]
    atmos_config["dt"] = dt_atmos

    # set restart file to the initial file saved in this location if it is not nothing
    # TODO this is hardcoded and should be fixed once we have a better restart system
    if !isnothing(atmos_config["restart_file"])
        atmos_config["restart_file"] = replace(atmos_config["restart_file"], "active" => "0000")
    end
    return atmos_config
end


# flux calculation borrowed from atmos
"""
    CoupledMoninObukhov()
A modified version of a Monin-Obukhov surface for the Coupler, see the link below for more information
https://clima.github.io/SurfaceFluxes.jl/dev/SurfaceFluxes/#Monin-Obukhov-Similarity-Theory-(MOST)
"""
struct CoupledMoninObukhov end
"""
    coupler_surface_setup(::CoupledMoninObukhov, p, csf_sfc = (; T = nothing, z0m = nothing, z0b = nothing, beta = nothing, q_vap = nothing))

Sets up `surface_setup` as a `CC.Fields.Field` of `SurfaceState`s.
"""
function coupler_surface_setup(
    ::CoupledMoninObukhov,
    p,
    T = nothing,
    z0m = nothing,
    z0b = nothing,
    beta = nothing,
    q_vap = nothing,
)

    surface_state(z0m, z0b, T, beta, q_vap) = CA.SurfaceConditions.SurfaceState(;
        parameterization = CA.SurfaceConditions.MoninObukhov(; z0m, z0b),
        T,
        beta,
        q_vap,
    )
    surface_state_field = @. surface_state(z0m, z0b, T, beta, q_vap)
    return surface_state_field
end

"""
    get_new_cache(atmos_sim::ClimaAtmosSimulation, csf)

Returns a new `p` with the updated surface conditions.
"""
function get_new_cache(atmos_sim::ClimaAtmosSimulation, csf)
    if hasmoisture(atmos_sim.integrator)
        csf_sfc = (csf.T_S, csf.z0m_S, csf.z0b_S, csf.beta, csf.q_sfc)
    else
        csf_sfc = (csf.T_S, csf.z0m_S, csf.z0b_S, csf.beta)
    end

    p = atmos_sim.integrator.p

    coupler_sfc_setup = coupler_surface_setup(CoupledMoninObukhov(), p, csf_sfc...)

    p_names = propertynames(p)
    p_values = map(x -> x == :sfc_setup ? coupler_sfc_setup : getproperty(p, x), p_names)

    (; zip(p_names, p_values)...)
end

"""
    FluxCalculator.atmos_turbulent_fluxes_most!(atmos_sim::ClimaAtmosSimulation, csf)

Computes turbulent surface fluxes using ClimaAtmos's `update_surface_conditions!` and
and the Monin Obukhov Similarity Theory. This
requires that we define a new temporary parameter Tuple, `new_p`, and save the new surface state
in it. We do not want `new_p` to live in the atmospheric model permanently, because that would also
trigger flux calculation during Atmos `Interfacer.step!`. We only want to trigger this once per coupling
timestep from ClimaCoupler.

For debigging atmos, we can set the following atmos defaults:
 csf.z0m_S .= 1.0e-5
 csf.z0b_S .= 1.0e-5
 csf.beta .= 1
 csf = merge(csf, (;q_sfc = nothing))
"""
function FluxCalculator.atmos_turbulent_fluxes_most!(atmos_sim::ClimaAtmosSimulation, csf)

    if isnothing(atmos_sim.integrator.p.sfc_setup) # trigger flux calculation if not done in Atmos internally
        new_p = get_new_cache(atmos_sim, csf)
        CA.SurfaceConditions.update_surface_conditions!(atmos_sim.integrator.u, new_p, atmos_sim.integrator.t)
        atmos_sim.integrator.p.precomputed.sfc_conditions .= new_p.precomputed.sfc_conditions
    end
end

"""
    get_thermo_params(sim::ClimaAtmosSimulation)

Returns the thermodynamic parameters from the atmospheric model simulation object.
"""
get_thermo_params(sim::ClimaAtmosSimulation) = CAP.thermodynamics_params(sim.integrator.p.params)

"""
    dss_state!(sim::ClimaAtmosSimulation)

Perform DSS on the state of a component simulation, intended to be used
before the initial step of a run. This method acts on atmosphere simulations.
These sims don't store a dss buffer in their cache, so we must allocate
one here.
"""
function dss_state!(sim::ClimaAtmosSimulation)
    Y = sim.integrator.u
    for key in propertynames(Y)
        field = getproperty(Y, key)
        buffer = CC.Spaces.create_dss_buffer(field)
        CC.Spaces.weighted_dss!(field, buffer)
    end
end

"""
    FluxCalculator.water_albedo_from_atmosphere!(atmos_sim::ClimaAtmosSimulation, direct_albedo::CC.Fields.Field, diffuse_albedo::CC.Fields.Field)

Extension to calculate the water surface albedo from wind speed and insolation. It can be used for prescribed ocean and lakes.
"""
function FluxCalculator.water_albedo_from_atmosphere!(
    atmos_sim::ClimaAtmosSimulation,
    direct_albedo::CC.Fields.Field,
    diffuse_albedo::CC.Fields.Field,
)

    Y = atmos_sim.integrator.u
    p = atmos_sim.integrator.p
    t = atmos_sim.integrator.t

    rrtmgp_model = atmos_sim.integrator.p.radiation.rrtmgp_model
    FT = eltype(Y)
    λ = FT(0) # spectral wavelength (not used for now)

    # update for current zenith angle
    bottom_coords = CC.Fields.coordinate_field(CC.Spaces.level(Y.c, 1))
    μ = CC.Fields.array2field(rrtmgp_model.cos_zenith, axes(bottom_coords))
    FT = eltype(atmos_sim.integrator.u)
    α_model = CA.RegressionFunctionAlbedo{FT}()

    # set the direct and diffuse surface albedos
    direct_albedo .= CA.surface_albedo_direct(α_model).(λ, μ, LinearAlgebra.norm.(CC.Fields.level(Y.c.uₕ, 1)))
    diffuse_albedo .= CA.surface_albedo_diffuse(α_model).(λ, μ, LinearAlgebra.norm.(CC.Fields.level(Y.c.uₕ, 1)))

end
