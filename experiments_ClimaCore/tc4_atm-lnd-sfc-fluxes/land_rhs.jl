"""
    SoilModelCoupler{FT, domain, em <: AbstractSoilModel, hm <: AbstractSoilModel, bc, A,B}
The model type for the soil model.
# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct SoilModel{
    FT,
    dm <: AbstractVerticalDomain{FT},
    em <: AbstractSoilComponentModel,
    hm <: AbstractSoilComponentModel,
    bc,
    A,
    B,
} <: AbstractModel
    domain::dm
    "Soil energy model - prescribed or dynamics"
    energy_model::em
    "Soil hydrology model - prescribed or dynamic"
    hydrology_model::hm
    "Boundary conditions tuple"
    boundary_conditions::bc
    "Soil parameters"
    soil_param_set::A
    "Earth parameter set"
    earth_param_set::B
    "name"
    name::Symbol = :soil
    "variables"
    variables::Tuple = (:ϑ_l, :θ_i, :ρe_int)
end

"""
    make_rhs!(energy::SoilEnergyModel, hydrology::PrescribedHydrologyModel, model::SoilModel)
"""

function make_rhs!(
    energy::SoilEnergyModel,
    hydrology::PrescribedHydrologyModel,
    model::SoilModelCoupler,
)
    function rhs!(dY, Y, _, t)
        dϑ_l = dY.ϑ_l
        dθ_i = dY.θ_i
        dρe_int = dY.ρe_int
        ϑ_l = Y.ϑ_l
        θ_i = Y.θ_i
        ρe_int = Y.ρe_int

        cspace = axes(ϑ_l)
        zc = Fields.coordinate_field(cspace)
        FT = eltype(zc)

        # boundary conditions and parameters
        bc = model.boundary_conditions
        top_heat_flux, btm_heat_flux = compute_vertical_flux(bc.top.energy),compute_vertical_flux(bc.bottom.energy)
        sp = model.soil_param_set
        param_set = model.earth_param_set
        @unpack ν, ρc_ds, κ_sat_unfrozen, κ_sat_frozen = sp

        # update water content based on prescribed profiles, set RHS to zero.
        ϑ_l = hydrology.ϑ_l_profile.(zc, t)
        θ_i = hydrology.θ_i_profile.(zc, t)
        dθ_i = Fields.zeros(FT, cspace)
        dρe_int = Fields.zeros(FT, space)

        # Compute center values of everything
        ρc_s = volumetric_heat_capacity.(θ_l, θ_i, ρc_ds, Ref(param_set))
        T = temperature_from_ρe_int.(ρe_int, θ_i, ρc_s, Ref(param_set))
        κ_dry = k_dry(param_set, sp)
        S_r = relative_saturation.(θ_l, θ_i, ν)
        kersten = kersten_number.(θ_i, S_r, Ref(sp))
        κ_sat =
            saturated_thermal_conductivity.(
                θ_l,
                θ_i,
                κ_sat_unfrozen,
                κ_sat_frozen,
            )
        κ = thermal_conductivity.(κ_dry, kersten, κ_sat)

        # rhs operators
        interpc2f = Operators.InterpolateC2F()
        gradc2f_heat = Operators.GradientC2F()
        divf2c_heat = Operators.DivergenceF2C(
            top = Operators.SetValue(top_heat_flux),
            bottom = Operators.SetValue(btm_heat_flux),
        )
        @. dρe_int = -divf2c_heat(-interpc2f(κ) * gradc2f_heat(T))
        return dY
    end
    return rhs!
end
