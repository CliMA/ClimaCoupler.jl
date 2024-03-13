# these extensions add extra diagnostics to the atmos model output
import ClimaAtmos.Diagnostics as CAD

import ClimaAtmos.Diagnostics: add_diagnostic_variable!

"""
    add_diagnostic_variable!(short_name::String, long_name::String, standard_name::String, units::String, comments::String, compute!::Function)

This is an extension of the `add_diagnostic_variable!` function from ClimaAtmos.Diagnostics, and it adds the specified variable
to the ClimaAtmos.Diagnostics.ALL_DIAGNOSTICS dictionary of possible diagnostic variables.
The `compute!` function is called at every atmos time step to compute the diagnostic variable.

To output these variables, short_name needs to be specified under diagnostics in the required yml file.
"""
### TODO This has been fixed in a recent update

#add_diagnostic_variable!(
#    short_name = "mse",
#    long_name = "Moist static energy",
#    standard_name = "moist_static_energy",
#    units = "J/kg",
#    comments = "Moist static energy",
#    compute! = (out, state, cache, time) -> begin
#        (; params) = cache
#        (; ᶜts) = cache.precomputed
#        c_space = axes(state.c)
#        thermo_params = CAP.thermodynamics_params(params)
#        e_pot = CAP.grav(params) .* Fields.coordinate_field(c_space).z
#        if isnothing(out)
#            return TD.moist_static_energy.(thermo_params, ᶜts, e_pot)
#        else
#            out .= TD.moist_static_energy.(thermo_params, ᶜts, e_pot)
#        end
#    end,
#)
#
#add_diagnostic_variable!(
#    short_name = "lr",
#    long_name = "Lapse rate",
#    standard_name = "lapse_rate",
#    units = "K/m",
#    comments = "Lapse rate",
#    compute! = (out, state, cache, time) -> begin
#        (; params) = cache
#        (; ᶜts) = cache.precomputed
#        thermo_params = CAP.thermodynamics_params(params)
#        ᶜT = @. TD.air_temperature(thermo_params, ᶜts)
#        if isnothing(out)
#            return ClimaCore.Geometry.WVector.(CAD.ᶜgradᵥ.(CAD.ᶠinterp.(ᶜT))).components.data.:1
#        else
#            out .= ClimaCore.Geometry.WVector.(CAD.ᶜgradᵥ.(CAD.ᶠinterp.(ᶜT))).components.data.:1
#        end
#
#    end,
#)
#
#add_diagnostic_variable!(
#    short_name = "ediff",
#    long_name = "Eddy diffusivity",
#    standard_name = "eddy_diffusivity",
#    units = "m2/s",
#    comments = "Eddy diffusivity consistent with the VerticalDiffusion scheme in ClimaAtmos",
#    compute! = (out, state, cache, time) -> begin
#        (; ᶜp) = cache.precomputed
#        (; C_E) = cache.atmos.vert_diff
#        interior_uₕ = Fields.level(state.c.uₕ, 1)
#        ᶠp = ᶠK_E = cache.scratch.ᶠtemp_scalar
#        Fields.bycolumn(axes(ᶜp)) do colidx
#            @. ᶠp[colidx] = CAD.ᶠinterp(ᶜp[colidx])
#            ᶜΔz_surface = Fields.Δz_field(interior_uₕ)
#            @. ᶠK_E[colidx] = CA.eddy_diffusivity_coefficient(
#                C_E,
#                CA.norm(interior_uₕ[colidx]),
#                ᶜΔz_surface[colidx] / 2,
#                ᶠp[colidx],
#            )
#        end
#        if isnothing(out)
#            return CAD.ᶜinterp.(ᶠK_E)
#        else
#            out .= CAD.ᶜinterp.(ᶠK_E)
#        end
#
#    end,
#)
