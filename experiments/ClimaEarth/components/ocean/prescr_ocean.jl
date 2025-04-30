import ClimaUtilities.TimeVaryingInputs: TimeVaryingInput, evaluate!
import ClimaUtilities.ClimaArtifacts: @clima_artifact
import Interpolations # triggers InterpolationsExt in ClimaUtilities
import Thermodynamics as TD
import ClimaCoupler: Checkpointer, FieldExchanger, Interfacer

"""
    PrescribedOceanSimulation{C}

Sea surface temperature (SST) is prescribed and read in directly from a file
at each timestep.

Ocean roughness follows the https://github.com/NOAA-GFDL/ice_param/blob/main/ocean_rough.F90#L47.

The cache is expected to contain the following variables:
- `T_sfc` (surface temperature [K])
- `ρ_sfc` (surface air density [kg / m3])
- `z0m` (roughness length for momentum [m])
- `z0b` (roughness length for tracers [m])
- `beta` (evaporation scaling factor)
- `α_direct` (direct albedo)
- `α_diffuse` (diffuse albedo)
- `area_fraction` (fraction of the grid cell covered by the ocean)
- `phase` (phase of the water used to calculate surface humidity)
- `thermo_params` (thermodynamic parameters)
- `SST_timevaryinginput` (TimeVaryingInput object containing SST data)
"""
struct PrescribedOceanSimulation{C} <: Interfacer.AbstractSurfaceStub
    cache::C
end

"""
    PrescribedOceanSimulation(
        ::Type{FT},
        space,
        start_date,
        area_fraction,
        thermo_params,
        comms_ctx;
        z0m = FT(5.8e-5),
        z0b = FT(5.8e-5),
        beta = FT(1),
        α_direct_val = FT(0.06),
        α_diffuse_val = FT(0.06),
    )

Initialize the `PrescribedOceanSimulation` object with all required cache fields,
and reading in prescribed SST data.
"""
function PrescribedOceanSimulation(
    ::Type{FT},
    space,
    start_date,
    t_start,
    area_fraction,
    thermo_params,
    comms_ctx;
    z0m = FT(5.8e-5),
    z0b = FT(5.8e-5),
    beta = FT(1),
    α_direct_val = FT(0.06),
    α_diffuse_val = FT(0.06),
) where {FT}
    # Read in initial SST data
    sst_data = try
        joinpath(@clima_artifact("historical_sst_sic", comms_ctx), "MODEL.SST.HAD187001-198110.OI198111-202206.nc")
    catch error
        @warn "Using lowres SST. If you want the higher resolution version, you have to obtain it from ClimaArtifacts"
        joinpath(
            @clima_artifact("historical_sst_sic_lowres", comms_ctx),
            "MODEL.SST.HAD187001-198110.OI198111-202206_lowres.nc",
        )
    end

    SST_timevaryinginput = TimeVaryingInput(
        sst_data,
        "SST",
        space,
        reference_date = start_date,
        file_reader_kwargs = (; preprocess_func = (data) -> data + FT(273.15),), ## convert to Kelvin
    )

    SST_init = zeros(space)
    evaluate!(SST_init, SST_timevaryinginput, t_start)

    # Create the cache
    cache = (;
        T_sfc = SST_init,
        ρ_sfc = zeros(space),
        z0m = z0m,
        z0b = z0b,
        beta = beta,
        α_direct = ones(space) .* α_direct_val,
        α_diffuse = ones(space) .* α_diffuse_val,
        area_fraction = area_fraction,
        phase = TD.Liquid(),
        thermo_params = thermo_params,
        SST_timevaryinginput = SST_timevaryinginput,
    )
    return PrescribedOceanSimulation(cache)
end

## Extensions of Interfacer and FieldExchanger functions

"""
    Interfacer.step!(sim::PrescribedOceanSimulation, t)

Update the cached surface temperature field using the prescribed data
at each timestep.
"""
function Interfacer.step!(sim::PrescribedOceanSimulation, t)
    evaluate!(sim.cache.T_sfc, sim.cache.SST_timevaryinginput, t)
end

function Checkpointer.get_model_cache(sim::PrescribedOceanSimulation)
    return sim.cache
end

function Checkpointer.restore_cache!(sim::PrescribedOceanSimulation, new_cache)
    old_cache = Checkpointer.get_model_cache(sim)
    for p in propertynames(old_cache)
        if getproperty(old_cache, p) isa Field
            ArrayType = ClimaComms.array_type(getproperty(old_cache, p))
            parent(getproperty(old_cache, p)) .= ArrayType(parent(getproperty(new_cache, p)))
        end
    end
end
