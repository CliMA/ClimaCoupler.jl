import ClimaAtmos.TurbulenceConvection as TC
import ClimaAtmos.TurbulenceConvection.Parameters as TCP
import CLIMAParameters as CP
import RRTMGP.Parameters as RP
import SurfaceFluxes as SF
import SurfaceFluxes.UniversalFunctions as UF
import Insolation.Parameters as IP
import Thermodynamics as TD
import CloudMicrophysics as CM

import ClimaAtmos: create_parameter_set

# Grachev
function get_uf_params(::UF.GrachevType, toml_dict::CP.AbstractTOMLDict)
    aliases = [
        "Pr_0_Grachev",
        "a_m_Grachev",
        "a_h_Grachev",
        "b_m_Grachev",
        "b_h_Grachev",
        "c_h_Grachev",
        "ζ_a_Grachev",
        "γ_Grachev",
    ]
    pairs = CP.get_parameter_values!(toml_dict, aliases, "UniversalFunctions")
    pairs = (; pairs...) # convert to NamedTuple
    pairs = (;
        Pr_0 = pairs.Pr_0_Grachev,
        a_m = pairs.a_m_Grachev,
        a_h = pairs.a_h_Grachev,
        b_m = pairs.b_m_Grachev,
        b_h = pairs.b_h_Grachev,
        c_h = pairs.c_h_Grachev,
        ζ_a = pairs.ζ_a_Grachev,
        γ = pairs.γ_Grachev,
    )

    UF.GrachevParams{FT}(; pairs...)
end

# Beljaars
function get_uf_params(::UF.BeljaarsType, toml_dict::CP.AbstractTOMLDict)
    aliases = [
        "Pr_0_Beljaars",
        "a_m_Beljaars",
        "a_h_Beljaars",
        "b_m_Beljaars",
        "b_h_Beljaars",
        "c_h_Beljaars",
        "c_m_Beljaars",
        "d_h_Beljaars",
        "d_m_Beljaars",
        "ζ_a_Beljaars",
        "γ_Beljaars",
    ]
    pairs = CP.get_parameter_values!(toml_dict, aliases, "UniversalFunctions")
    pairs = (; pairs...) # convert to NamedTuple
    pairs = (;
        Pr_0 = pairs.Pr_0_Beljaars,
        a_m = pairs.a_m_Beljaars,
        a_h = pairs.a_h_Beljaars,
        b_m = pairs.b_m_Beljaars,
        b_h = pairs.b_h_Beljaars,
        c_m = pairs.c_m_Beljaars,
        c_h = pairs.c_h_Beljaars,
        d_m = pairs.d_m_Beljaars,
        d_h = pairs.d_h_Beljaars,
        ζ_a = pairs.ζ_a_Beljaars,
        γ = pairs.γ_Beljaars,
    )

    UF.BeljaarsParams{FT}(; pairs...)
end

# Cheng
function get_uf_params(::UF.ChengType, toml_dict::CP.AbstractTOMLDict)
    aliases = [
        "Pr_0_Cheng",
        "a_m_Cheng",
        "a_h_Cheng",
        "b_m_Cheng",
        "b_h_Cheng",
        "ζ_a_Cheng",
        "γ_Cheng",
    ]
    pairs = CP.get_parameter_values!(toml_dict, aliases, "UniversalFunctions")
    pairs = (; pairs...) # convert to NamedTuple
    pairs = (;
        Pr_0 = pairs.Pr_0_Cheng,
        a_m = pairs.a_m_Cheng,
        a_h = pairs.a_h_Cheng,
        b_m = pairs.b_m_Cheng,
        b_h = pairs.b_h_Cheng,
        ζ_a = pairs.ζ_a_Cheng,
        γ = pairs.γ_Cheng,
    )
    UF.ChengParams{FT}(; pairs...)
end

# Holtslag
function get_uf_params(::UF.HoltslagType, toml_dict::CP.AbstractTOMLDict)
    aliases = [
        "Pr_0_Holtslag",
        "a_m_Holtslag",
        "a_h_Holtslag",
        "b_m_Holtslag",
        "b_h_Holtslag",
        "c_h_Holtslag",
        "c_m_Holtslag",
        "d_h_Holtslag",
        "d_m_Holtslag",
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

    UF.HoltslagParams{FT}(; pairs...)
end

# Businger
function get_uf_params(::UF.BusingerType, toml_dict::CP.AbstractTOMLDict)
    aliases = [
        "Pr_0_Businger",
        "a_m_Businger",
        "a_h_Businger",
        "ζ_a_Businger",
        "γ_Businger",
    ]
    pairs = CP.get_parameter_values!(toml_dict, aliases, "UniversalFunctions")
    pairs = (; pairs...) # convert to NamedTuple
    pairs = (;
        Pr_0 = pairs.Pr_0_Businger,
        a_m = pairs.a_m_Businger,
        a_h = pairs.a_h_Businger,
        ζ_a = pairs.ζ_a_Businger,
        γ = pairs.γ_Businger,
    )
    UF.BusingerParams{FT}(; pairs...)
end

# Gryanik
function get_uf_params(::UF.GryanikType, toml_dict::CP.AbstractTOMLDict)
    aliases = [
        "Pr_0_Gryanik",
        "a_m_Gryanik",
        "a_h_Gryanik",
        "b_m_Gryanik",
        "b_h_Gryanik",
        "ζ_a_Gryanik",
        "γ_Gryanik",
    ]
    pairs = CP.get_parameter_values!(toml_dict, aliases, "UniversalFunctions")
    pairs = (; pairs...) # convert to NamedTuple
    pairs = (;
        Pr_0 = pairs.Pr_0_Gryanik,
        a_m = pairs.a_m_Gryanik,
        a_h = pairs.a_h_Gryanik,
        b_m = pairs.a_m_Gryanik,
        b_h = pairs.a_h_Gryanik,
        ζ_a = pairs.ζ_a_Gryanik,
        γ = pairs.γ_Gryanik,
    )
    UF.GryanikParams{FT}(; pairs...)
end

function create_parameter_set(config::CA.AtmosConfig)
    # Helper function that creates a parameter struct. If a struct has nested
    # parameter structs, they must be passed to subparam_structs as a NamedTuple.
    function create_parameter_struct(param_struct; subparam_structs = (;))
        aliases = string.(fieldnames(param_struct))
        aliases = setdiff(aliases, string.(propertynames(subparam_structs)))
        pairs = CP.get_parameter_values!(toml_dict, aliases)
        # Workaround for setting τ_precip = dt
        if parsed_args["override_τ_precip"] &&
           param_struct == CM.Parameters.CloudMicrophysicsParameters
            pairs = (;
                pairs...,
                τ_precip = FT(CA.time_to_seconds(parsed_args["dt"])),
            )
        end
        return param_struct{FT, typeof.(values(subparam_structs))...}(;
            pairs...,
            subparam_structs...,
        )
    end

    (; toml_dict, parsed_args) = config
    FT = CP.float_type(toml_dict)

    thermo_params =
        create_parameter_struct(TD.Parameters.ThermodynamicsParameters)
    rrtmgp_params = create_parameter_struct(RP.RRTMGPParameters)
    insolation_params = create_parameter_struct(IP.InsolationParameters)
    microphys_params = create_parameter_struct(
        CM.Parameters.CloudMicrophysicsParameters;
        subparam_structs = (; thermo_params),
    )

    if config_dict["universal_function"] == "Businger" # TODO: config_dict must be global here (port to Atmos / SF)
        universal_function = UF.Businger()
    elseif config_dict["universal_function"] == "Gryanik"
        universal_function = UF.Gryanik()
    elseif config_dict["universal_function"] == "Grachev"
        universal_function = UF.Grachev()
    elseif config_dict["universal_function"] == "Beljaars"
        universal_function = UF.Beljaars()
    elseif config_dict["universal_function"] == "Cheng"
        universal_function = UF.Cheng()
    elseif config_dict["universal_function"] == "Holtslag"
        universal_function = UF.Holtslag()
    else
        error("Unknown universal function: $(config_dict["universal_function"])")
    end

    ufp = get_uf_params(universal_function, toml_dict)

    surf_flux_params = create_parameter_struct(
        SF.Parameters.SurfaceFluxesParameters;
        subparam_structs = (; ufp, thermo_params),
    )
    turbconv_params = create_parameter_struct(
        TCP.TurbulenceConvectionParameters;
        subparam_structs = (; microphys_params, surf_flux_params),
    )
    return create_parameter_struct(
        CAP.ClimaAtmosParameters;
        subparam_structs = (;
            thermodynamics_params = thermo_params,
            rrtmgp_params,
            insolation_params,
            microphysics_params = microphys_params,
            surface_fluxes_params = surf_flux_params,
            turbconv_params,
        ),
    )
end
