using ClimaCore.Geometry: ⊗
using ClimaCore.Utilities: half, PlusHalf
import ClimaAtmos: get_surface_fluxes!


"""
    calculate_surface_fluxes_in_atmos!(i, surface_fractions)

Calculates surface fluxes using adapter function `get_surface_fluxes!`
from ClimaAtmos that calls `SurfaceFluxes.jl`. The coupler updates in
atmos model cache fluxes at each coupling timestep.

- TODO: generalize interface for regridding and take land state out of atmos's integrator.p
"""

function calculate_surface_fluxes_in_atmos!(i, surface_fractions)

    thermo_params = CAP.thermodynamics_params(i.p.params)

    Fields.bycolumn(axes(i.p.ts_sfc)) do colidx
        ClimaAtmos.set_surface_thermo_state!(
            ClimaAtmos.Decoupled(),
            i.p.surface_scheme.sfc_thermo_state_type,
            i.p.ts_sfc[colidx],
            i.p.T_sfc[colidx],
            Spaces.level(i.p.ᶜts[colidx], 1),
            thermo_params,
            i.t,
        )

        get_surface_fluxes!(
            i.u,
            i.p,
            i.t,
            colidx,
            i.p.atmos.vert_diff,
        )

        # TODO remove correction below upon ClimaAtmos 0.10 release update, add beta where needed
        # TODO check that vapour fluxes (subplimation/condensation) are ok over ice
    end

end
