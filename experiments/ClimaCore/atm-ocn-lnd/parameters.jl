import CLIMAParameters as CP
import Thermodynamics.Parameters as TDP
import Insolation.Parameters as IP
import SurfaceFluxes.Parameters as SFP
import SurfaceFluxes.UniversalFunctions as UF


abstract type AbstractCouplerParameters end
const ACP = AbstractCouplerParameters

Base.@kwdef struct CouplerParameters{TP, SFP, IP} <: ACP
    thermo_params::TP
    surf_flux_params::SFP
    insol_params::IP
end

Base.eltype(::CouplerParameters{FT}) where {FT} = FT

function create_parameters(FT)
    toml_dict = CP.create_toml_dict(FT; dict_type = "alias")

    # Thermodynamics
    aliases = string.(fieldnames(TDP.ThermodynamicsParameters))
    param_pairs = CP.get_parameter_values!(toml_dict, aliases, "Thermodynamics")
    thermo_params = TDP.ThermodynamicsParameters{FT}(; param_pairs...)
    TP = typeof(thermo_params)

    # Insolation
    insol_aliases = string.(fieldnames(IP.InsolationParameters))
    insol_param_pairs =
        CP.get_parameter_values!(toml_dict, insol_aliases, "Insolation")
    insol_params = IP.InsolationParameters{FT}(; insol_param_pairs...)
    IPT = typeof(insol_params)

    # SurfaceFluxes
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
    ufp = UF.BusingerParams{FT}(; pairs...)
    UFP = typeof(ufp)

    pairs = CP.get_parameter_values!(
        toml_dict,
        ["von_karman_const"],
        "SurfaceFluxesParameters",
    )
    surf_flux_params =
        SFP.SurfaceFluxesParameters{FT, UFP, TP}(; pairs..., ufp, thermo_params)
    SFPS = typeof(surf_flux_params)

    # combine into CouplerParameters
    param_set = CouplerParameters{TP, SFPS, IPT}(;
        thermo_params,
        surf_flux_params,
        insol_params,
    )

    # add ClimaAtmos, LSM, ClimaOcean params
    # logfilepath = joinpath(@__DIR__, "logfilepath_$FT.toml")
    # CP.log_parameter_information(toml_dict, logfilepath)
    return param_set
end



# # wrapper methods:
# ρ_cloud_liq(ps::ACP) = ps.ρ_cloud_liq
# ρ_cloud_ice(ps::ACP) = ps.ρ_cloud_ice
# T_freeze(ps::ACP) = ps.T_freeze

# # Dependency parameter wrappers
# thermodynamic_parameters(ps::ACP) = ps.thermo_params
# surface_fluxes_parameters(ps::ACP) = ps.surf_flux_params
# insolation_parameters(ps::ACP) = ps.insol_params

# TODO: add other component parameters