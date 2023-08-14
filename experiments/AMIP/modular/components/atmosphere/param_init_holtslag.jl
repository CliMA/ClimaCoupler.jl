import ClimaAtmos.TurbulenceConvection as TC
import ClimaAtmos.TurbulenceConvection.Parameters as TCP
import CLIMAParameters as CP
import RRTMGP.Parameters as RP
import SurfaceFluxes as SF
import SurfaceFluxes.UniversalFunctions as UF
import ClimaCore
import ClimaCore as CC
import Insolation.Parameters as IP
import Thermodynamics as TD
import CloudMicrophysics as CM
import ClimaAtmos as CA
import ClimaAtmos: create_climaatmos_parameter_set

"""
    function create_climaatmos_parameter_set(
        toml_dict::CP.AbstractTOMLDict,
        parsed_args,
        overrides::NamedTuple = NamedTuple(),
    )

ClimaAtmos function that creates a parameter set for ClimaAtmos from a TOML dictionary and parsed arguments.
This overwriting version allows additional specification of parameters.
"""
function create_climaatmos_parameter_set(
    toml_dict::CP.AbstractTOMLDict,
    parsed_args,
    overrides::NamedTuple = NamedTuple(),
)
    FT = CP.float_type(toml_dict)
    FTD = FT # can change to Dual for testing duals

    aliases = string.(fieldnames(TD.Parameters.ThermodynamicsParameters))
    pairs = CP.get_parameter_values!(toml_dict, aliases, "Thermodynamics")
    pairs = CA.override_climaatmos_defaults((; pairs...), overrides)
    thermo_params = TD.Parameters.ThermodynamicsParameters{FTD}(; pairs...)
    TP = typeof(thermo_params)

    aliases = string.(fieldnames(CM.Parameters.CloudMicrophysicsParameters))
    aliases = setdiff(aliases, ["thermo_params"])
    pairs = CP.get_parameter_values!(toml_dict, aliases, "CloudMicrophysics")
    pairs = CA.override_climaatmos_defaults((; pairs...), overrides)
    microphys_params = CM.Parameters.CloudMicrophysicsParameters{FTD, TP}(;
        pairs...,
        thermo_params,
    )
    MP = typeof(microphys_params)

    aliases = [
        "Pr_0_Holtslag",
        "a_m_Holtslag",
        "a_h_Holtslag",
        "b_m_Holtslag",
        "b_h_Holtslag",
        "c_m_Holtslag",
        "c_h_Holtslag",
        "d_m_Holtslag",
        "d_h_Holtslag",
        "ζ_a_Holtslag",
        "γ_Holtslag",
    ]
    pairs = CP.get_parameter_values!(toml_dict, aliases, "UniversalFunctions")
    pairs = (; pairs...) # convert to NamedTuple
    pairs = (;
        Pr_0 = pairs.Pr_0_Holtslag,
        a_m = pairs.a_m_Holtslag,
        a_h = pairs.a_h_Holtslag,
        b_m = pairs.b_m_Holtslag,
        b_h = pairs.b_h_Holtslag,
        c_m = pairs.c_m_Holtslag,
        c_h = pairs.c_h_Holtslag,
        d_m = pairs.d_m_Holtslag,
        d_h = pairs.d_h_Holtslag,
        ζ_a = pairs.ζ_a_Holtslag,
        γ = pairs.γ_Holtslag,
    )
    pairs = CA.override_climaatmos_defaults((; pairs...), overrides)
    ufp = UF.HoltslagParams{FTD}(; pairs...)
    UFP = typeof(ufp)

    pairs = CP.get_parameter_values!(
        toml_dict,
        ["von_karman_const"],
        "SurfaceFluxesParameters",
    )
    pairs = CA.override_climaatmos_defaults((; pairs...), overrides)
    surf_flux_params = SF.Parameters.SurfaceFluxesParameters{FTD, UFP, TP}(;
        pairs...,
        ufp,
        thermo_params,
    )
    SFP = typeof(surf_flux_params)

    aliases = string.(fieldnames(TCP.TurbulenceConvectionParameters))
    pairs = CP.get_parameter_values!(toml_dict, aliases, "EDMF")
    pairs = CA.override_climaatmos_defaults((; pairs...), overrides)
    tc_params = TCP.TurbulenceConvectionParameters{FTD, MP, SFP}(;
        pairs...,
        microphys_params,
        surf_flux_params,
    )

    aliases = string.(fieldnames(RP.RRTMGPParameters))
    pairs = CP.get_parameter_values!(toml_dict, aliases, "RRTMGP")
    params = CA.override_climaatmos_defaults((; pairs...), overrides) # overrides
    rrtmgp_params = RP.RRTMGPParameters{FTD}(; params...)

    aliases = string.(fieldnames(IP.InsolationParameters))
    pairs = CP.get_parameter_values!(toml_dict, aliases, "Insolation")
    params = CA.override_climaatmos_defaults((; pairs...), overrides) # overrides
    insolation_params = IP.InsolationParameters{FTD}(; params...)

    pairs = CP.get_parameter_values!(
        toml_dict,
        ["Omega", "planet_radius", "astro_unit", "ΔT_y_dry", "ΔT_y_wet"],
        "ClimaAtmos",
    )
    pairs = (; pairs...) # convert to NamedTuple
    pairs = CA.override_climaatmos_defaults((; pairs...), overrides)

    param_set = CAP.ClimaAtmosParameters(;
        ug = FTD(1.0), # for Ekman problem
        vg = FTD(0.0), # for Ekman problem
        f = FTD(5e-5), # for Ekman problem
        Cd = FTD(0.01 / (2e2 / 30)), # for Ekman problem
        Omega = FTD(pairs.Omega),
        planet_radius = FTD(pairs.planet_radius),
        astro_unit = FTD(pairs.astro_unit),
        f_plane_coriolis_frequency = FTD(0),
        thermodynamics_params = thermo_params,
        microphysics_params = microphys_params,
        insolation_params = insolation_params,
        rrtmgp_params = rrtmgp_params,
        surfacefluxes_params = surf_flux_params,
        turbconv_params = tc_params,
        entr_coeff = FTD(parsed_args["entr_coeff"]),
        detr_coeff = FTD(parsed_args["detr_coeff"]),
        ΔT_y_dry = FTD(pairs.ΔT_y_dry),
        ΔT_y_wet = FTD(pairs.ΔT_y_wet),
    )
    # logfilepath = joinpath(@__DIR__, "logfilepath_$FT.toml")
    # CP.log_parameter_information(toml_dict, logfilepath)
    return param_set
end