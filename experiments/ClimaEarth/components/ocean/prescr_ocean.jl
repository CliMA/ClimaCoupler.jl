import ClimaUtilities.TimeVaryingInputs: TimeVaryingInput, evaluate!
import ClimaUtilities.ClimaArtifacts: @clima_artifact
import Interpolations # triggers InterpolationsExt in ClimaUtilities
import Thermodynamics as TD
import ClimaCoupler: Checkpointer, FieldExchanger, Interfacer
import Insolation
import Insolation.Parameters.InsolationParameters
import StaticArrays as SA
import Dates
import ClimaAtmos as CA # for albedo calculation
import LinearAlgebra

"""
    PrescribedOceanSimulation{C}

Sea surface temperature (SST) is prescribed and read in directly from a file
at each timestep.

Ocean roughness follows the https://github.com/NOAA-GFDL/ice_param/blob/main/ocean_rough.F90#L47.

The cache is expected to contain the following variables:
- `T_sfc` (surface temperature [K])
- `z0m` (roughness length for momentum [m])
- `z0b` (roughness length for tracers [m])
- `α_direct` (direct albedo)
- `α_diffuse` (diffuse albedo)
- `u_int`, `v_int` (the surface wind)
- `start_date`
- `area_fraction` (fraction of the grid cell covered by the ocean)
- `phase` (phase of the water used to calculate surface humidity)
- `thermo_params` (thermodynamic parameters)
- `SST_timevaryinginput` (TimeVaryingInput object containing SST data)
- `t` (the current time)
"""
struct PrescribedOceanSimulation{C} <: Interfacer.AbstractSurfaceStub
    cache::C
end

"""
    PrescribedOceanSimulation(
        ::Type{FT},
        space,
        start_date,
        t_start,
        coupled_param_dict,
        thermo_params,
        comms_ctx;
        z0m = FT(5.8e-5),
        z0b = FT(5.8e-5),
        α_direct_val = FT(0.06),
        α_diffuse_val = FT(0.06),
        sst_path::Union{Nothing, String} = nothing,
    )

Initialize the `PrescribedOceanSimulation` object with all required cache fields,
and reading in prescribed SST data.

The SST is read from the file specified by `sst_path`. If `sst_path` is `nothing`,
the model will use the default path from `ClimaArtifacts`.

"""
function PrescribedOceanSimulation(
    ::Type{FT},
    space,
    start_date,
    t_start,
    coupled_param_dict,
    thermo_params,
    comms_ctx;
    z0m = FT(5.8e-5),
    z0b = FT(5.8e-5),
    α_direct_val = FT(0.06),
    α_diffuse_val = FT(0.06),
    sst_path::Union{Nothing, String} = nothing,
) where {FT}
    # Read in initial SST data
    sst_data =
        isnothing(sst_path) ?
        try
            joinpath(
                @clima_artifact("historical_sst_sic", comms_ctx),
                "MODEL.SST.HAD187001-198110.OI198111-202206.nc",
            )
        catch error
            @warn "Using lowres SST. If you want the higher resolution version, you have to obtain it from ClimaArtifacts"
            joinpath(
                @clima_artifact("historical_sst_sic_lowres", comms_ctx),
                "MODEL.SST.HAD187001-198110.OI198111-202206_lowres.nc",
            )
        end : sst_path
    @info "PrescribedOcean: using SST file" sst_data

    C_to_K = coupled_param_dict["temperature_water_freeze"]
    SST_timevaryinginput = TimeVaryingInput(
        sst_data,
        "SST",
        space,
        reference_date = start_date,
        file_reader_kwargs = (; preprocess_func = (data) -> data + C_to_K,), ## convert Celsius to Kelvin
    )

    SST_init = zeros(space)
    evaluate!(SST_init, SST_timevaryinginput, t_start)

    # Create the cache
    cache = (;
        T_sfc = SST_init,
        z0m = z0m,
        z0b = z0b,
        α_direct = ones(space) .* α_direct_val,
        α_diffuse = ones(space) .* α_diffuse_val,
        u_int = zeros(space),
        v_int = zeros(space),
        area_fraction = ones(space),
        phase = TD.Liquid(),
        thermo_params = thermo_params,
        SST_timevaryinginput = SST_timevaryinginput,
        start_date = start_date,
        t = Ref(t_start),
    )
    return PrescribedOceanSimulation(cache)
end

## Extensions of Interfacer and FieldExchanger functions

Interfacer.get_field(sim::PrescribedOceanSimulation, ::Val{:surface_direct_albedo}) =
    sim.cache.α_direct
Interfacer.get_field(sim::PrescribedOceanSimulation, ::Val{:surface_diffuse_albedo}) =
    sim.cache.α_diffuse
Interfacer.get_field(sim::PrescribedOceanSimulation, ::Val{:emissivity}) =
    eltype(sim.cache.T_sfc)(1)
Interfacer.get_field(sim::PrescribedOceanSimulation, ::Val{:surface_temperature}) =
    sim.cache.T_sfc
Interfacer.get_field(sim::PrescribedOceanSimulation, ::Val{:area_fraction}) =
    sim.cache.area_fraction
Interfacer.get_field(sim::PrescribedOceanSimulation, ::Val{:roughness_momentum}) =
    sim.cache.z0m
Interfacer.get_field(sim::PrescribedOceanSimulation, ::Val{:roughness_buoyancy}) =
    sim.cache.z0b
# Specify COARE3 roughness model for PrescribedOceanSimulation
Interfacer.get_field(sim::PrescribedOceanSimulation, ::Val{:roughness_model}) = :coare3

function Interfacer.update_field!(
    sim::PrescribedOceanSimulation,
    ::Val{:u_int},
    field::CC.Fields.Field,
)
    Interfacer.remap!(sim.cache.u_int, field)
end

function Interfacer.update_field!(
    sim::PrescribedOceanSimulation,
    ::Val{:v_int},
    field::CC.Fields.Field,
)
    Interfacer.remap!(sim.cache.v_int, field)
end

"""
    Interfacer.step!(sim::PrescribedOceanSimulation, t)

Update the cached surface temperature field using the prescribed data
at each timestep.
"""
function Interfacer.step!(sim::PrescribedOceanSimulation, t)
    evaluate!(sim.cache.T_sfc, sim.cache.SST_timevaryinginput, t)
    sim.cache.t[] = t
end

function Checkpointer.get_model_cache(sim::PrescribedOceanSimulation)
    return sim.cache
end

"""
Extend Interfacer.add_coupler_fields! to add the fields required for PrescribedOceanSimulation.

The fields added are:
- `:u_int` (for water albedo calculation)
- `:v_int` (for water albedo calculation)
"""
function Interfacer.add_coupler_fields!(coupler_field_names, ::PrescribedOceanSimulation)
    ocean_coupler_fields = [:u_int, :v_int]
    push!(coupler_field_names, ocean_coupler_fields...)
end

function Checkpointer.restore_cache!(sim::PrescribedOceanSimulation, new_cache)
    old_cache = Checkpointer.get_model_cache(sim)
    for p in propertynames(old_cache)
        if getproperty(old_cache, p) isa CC.Fields.Field
            ArrayType = ClimaComms.array_type(getproperty(old_cache, p))
            parent(getproperty(old_cache, p)) .=
                ArrayType(parent(getproperty(new_cache, p)))
        end
    end
end

"""
    FieldExchanger.update_sim!(::PrescribedOceanSimulation, csf)

Update the wind velocity (needed for the turbulent flux calculation) and the
direct and diffuse albedos of the ocean.
"""
function FieldExchanger.update_sim!(sim::PrescribedOceanSimulation, csf)
    Interfacer.update_field!(sim, Val(:u_int), csf.u_int)
    Interfacer.update_field!(sim, Val(:v_int), csf.v_int)

    # Update the direct and diffuse albedos with the new atmospheric wind
    set_albedos!(sim, sim.cache.t[])
end

"""
    FieldExchanger.import_atmos_fields!(csf, sim::PrescribedOceanSimulation, atmos_sim)

Import wind from atmos.
"""
function FieldExchanger.import_atmos_fields!(csf, sim::PrescribedOceanSimulation, atmos_sim)
    Interfacer.get_field!(csf.u_int, atmos_sim, Val(:u_int))
    Interfacer.get_field!(csf.v_int, atmos_sim, Val(:v_int))
end

"""
    function set_albedos!(sim::PrescribedOceanSimulation, t)

Set the direct and diffuse albedos of the ocean based on the current date and
the atmospheric wind. The albedos are calculated using the `surface_albedo_direct`
and `surface_albedo_diffuse` functions from the `ClimaAtmos` module, so this
is dependent on running with `ClimaAtmosSimulation` as the atmosphere simulation.
"""
function set_albedos!(sim::PrescribedOceanSimulation, t)
    p = sim.cache
    FT = CC.Spaces.undertype(axes(sim.cache.T_sfc))

    # Compute the current date
    current_date =
        t isa ClimaUtilities.TimeManager.ITime ? date(t) : p.start_date + Dates.Second(t)

    insolation_params = InsolationParameters(FT)
    d, δ, η_UTC =
        FT.(Insolation.helper_instantaneous_zenith_angle(current_date, insolation_params))

    # Get the atmospheric wind vector and the cosine of the zenith angle
    surface_coords = CC.Fields.coordinate_field(axes(sim.cache.T_sfc))
    insolation_tuple =
        Insolation.instantaneous_zenith_angle.(
            d,
            δ,
            η_UTC,
            surface_coords.long,
            surface_coords.lat,
        ) # the tuple is (zenith angle, azimuthal angle, earth-sun distance)
    zenith_angle = insolation_tuple.:1
    wind_atmos = LinearAlgebra.norm.(CC.Geometry.Covariant12Vector.(p.u_int, p.v_int)) # wind vector from components
    λ = FT(0) # spectral wavelength (not used for now)
    max_zenith_angle = FT(π) / 2 - eps(FT)
    cos_zenith = @. cos(min(zenith_angle, max_zenith_angle)) # cosine of the zenith angle

    # Use the albedo model from ClimaAtmos
    α_model = CA.RegressionFunctionAlbedo{FT}()
    Interfacer.update_field!(
        sim,
        Val(:surface_direct_albedo),
        CA.surface_albedo_direct(α_model).(λ, cos_zenith, wind_atmos),
    )
    Interfacer.update_field!(
        sim,
        Val(:surface_diffuse_albedo),
        CA.surface_albedo_diffuse(α_model).(λ, cos_zenith, wind_atmos),
    )

    return nothing
end
