using ClimaCore.Geometry: ⊗
using ClimaCore.Utilities: half, PlusHalf
import ClimaAtmos: get_surface_fluxes!

"""
    set_ρ_sfc!(ρ_sfc, T_S, integrator)

sets the value of the ρ_sfc field based on the temperature of the surface,
the temperature of the atmosphere at the lowest level, and the heigh
of the lowest level.
"""
function set_ρ_sfc!(ρ_sfc, T_S, integrator)
    ts = integrator.p.ᶜts
    thermo_params = CAP.thermodynamics_params(integrator.p.params)
    ts_int = Spaces.level(ts, 1)
    parent(ρ_sfc) .= parent(ρ_sfc_at_point.(thermo_params, ts_int, swap_space!(zeros(axes(ts_int)), T_S)))
end

"""
    ρ_sfc_at_point(params, ts_int, T_sfc)

Computes the surface density at a point given the atmospheric state
at the lowest level, the surface temperature, and the assumption of
an ideal gas and hydrostatic balance.

Required because the surface models do not compute air density as a
variable.
"""
function ρ_sfc_at_point(params, ts_int, T_sfc)
    T_int = TD.air_temperature(params, ts_int)
    Rm_int = TD.gas_constant_air(params, ts_int)
    ρ_air = TD.air_density(params, ts_int)
    ρ_sfc = ρ_air * (T_sfc / T_int)^(TD.cv_m(params, ts_int) / Rm_int)  # use ideal gas law and hydrostatic balance to extrapolate for surface density
    return ρ_sfc
end

"""
    calculate_surface_fluxes_atmos_grid!(integrator)

Calculates surface fluxes using adapter function `get_surface_fluxes!`
from ClimaAtmos that calls `SurfaceFluxes.jl`. The coupler updates in
atmos model cache fluxes at each coupling timestep.

- TODO: generalize interface for regridding and take land state out of atmos's integrator.p
"""

function calculate_surface_fluxes_atmos_grid!(integrator, info_sfc)
    Y = integrator.u
    p = integrator.p
    t = integrator.t
    ice_fraction = info_sfc.ice_fraction

    Fields.bycolumn(axes(Y.c.uₕ)) do colidx
        get_surface_fluxes!(Y, p, colidx)
        # corrections (accounting for inhomogeneous surfaces)
        @. p.dif_flux_energy_bc[colidx] = # todo: get rid - shouldn't make any difference anyway
            Geometry.WVector(correct_e_over_ice(p.surface_conditions[colidx], ice_fraction[colidx]))
        @. p.dif_flux_ρq_tot_bc[colidx] =
            Geometry.WVector(correct_q_over_ice(p.surface_conditions[colidx], ice_fraction[colidx]))

    end
end
correct_e_over_ice(surface_conditions, ice_fraction) =
    .-surface_conditions.shf .- surface_conditions.lhf .* (FT(1) .- ice_fraction)
correct_q_over_ice(surface_conditions, ice_fraction) = .-surface_conditions.E .* (FT(1) .- ice_fraction)
