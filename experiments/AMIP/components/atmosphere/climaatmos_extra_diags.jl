# these extensions add extra diagnostics to the atmos model output
import ClimaAtmos.Diagnostics as CAD

"""
    add_diagnostic_variable!(short_name::String, long_name::String, standard_name::String, units::String, comments::String, compute!::Function)

This is an extension of the `add_diagnostic_variable!` function from ClimaAtmos.Diagnostics, and it adds the specified variable
to the ClimaAtmos.Diagnostics.ALL_DIAGNOSTICS dictionary of possible diagnostic variables.
The `compute!` function is called at every atmos time step to compute the diagnostic variable.

To output these variables, short_name needs to be specified under diagnostics in the required yml file.
"""

CAD.add_diagnostic_variable!(
    short_name = "mse",
    long_name = "Moist static energy",
    standard_name = "moist_static_energy",
    units = "J/kg",
    comments = "Moist static energy",
    compute! = (out, state, cache, time) -> begin
        (; params) = cache
        (; ᶜts) = cache.precomputed
        c_space = axes(state.c)
        thermo_params = CAP.thermodynamics_params(params)
        e_pot = CAP.grav(params) .* CC.Fields.coordinate_field(c_space).z
        if isnothing(out)
            return TD.moist_static_energy.(thermo_params, ᶜts, e_pot)
        else
            out .= TD.moist_static_energy.(thermo_params, ᶜts, e_pot)
        end
    end,
)

CAD.add_diagnostic_variable!(
    short_name = "lr",
    long_name = "Lapse rate",
    standard_name = "lapse_rate",
    units = "K/m",
    comments = "Lapse rate",
    compute! = (out, state, cache, time) -> begin
        (; params) = cache
        (; ᶜts) = cache.precomputed
        thermo_params = CAP.thermodynamics_params(params)
        ᶜT = @. TD.air_temperature(thermo_params, ᶜts)
        if isnothing(out)
            return CC.Geometry.WVector.(CAD.ᶜgradᵥ.(CAD.ᶠinterp.(ᶜT))).components.data.:1
        else
            out .= CC.Geometry.WVector.(CAD.ᶜgradᵥ.(CAD.ᶠinterp.(ᶜT))).components.data.:1
        end

    end,
)

CAD.add_diagnostic_variable!(
    short_name = "ediff",
    long_name = "Eddy diffusivity",
    standard_name = "eddy_diffusivity",
    units = "m2/s",
    comments = "Eddy diffusivity consistent with the VerticalDiffusion scheme in ClimaAtmos",
    compute! = (out, state, cache, time) -> begin
        (; ᶜp) = cache.precomputed
        (; C_E) = cache.atmos.vert_diff
        interior_uₕ = CC.Fields.level(state.c.uₕ, 1)
        ᶠp = ᶠK_E = cache.scratch.ᶠtemp_scalar
        CC.Fields.bycolumn(axes(ᶜp)) do colidx
            @. ᶠp[colidx] = CAD.ᶠinterp(ᶜp[colidx])
            ᶜΔz_surface = CC.Fields.Δz_field(interior_uₕ)
            @. ᶠK_E[colidx] = CA.eddy_diffusivity_coefficient(
                C_E,
                CA.norm(interior_uₕ[colidx]),
                ᶜΔz_surface[colidx] / 2,
                ᶠp[colidx],
            )
        end
        if isnothing(out)
            return CAD.ᶜinterp.(ᶠK_E)
        else
            out .= CAD.ᶜinterp.(ᶠK_E)
        end

    end,
)

###
# Add the mass streamfunction diagnostic variable
###

# TODO: this should be integrated from top-down!
CAD.add_diagnostic_variable!(
    short_name = "mass_streamfunction",
    long_name = "Meridional Mass Streamfunction (Hartmann Eq 6.9)",
    standard_name = "meridional_mass_streamfunction",
    units = "kg m^-1 s^-1",
    compute! = (out, state, cache, time) -> begin
        FT = eltype(state)
        lat = CC.Fields.coordinate_field(axes(state.f)).lat
        ρv = ClimaCore.Geometry.UVVector.(state.c.uₕ).components.data.:2 .* state.c.ρ
        ᶠ∫_0_z_ρv = cache.scratch.ᶠtemp_scalar
        ClimaCore.Fields.bycolumn(axes(ρv)) do colidx
            ClimaCore.Operators.column_integral_indefinite!(ᶠ∫_0_z_ρv[colidx], ρv[colidx])
        end
        twoπa = FT(2π * 6371e3)
        if isnothing(out)
            return ᶠ∫_0_z_ρv .* twoπa .* cosd.(lat)
        else
            out .= ᶠ∫_0_z_ρv .* twoπa .* cosd.(lat)
        end
    end,
)

CAD.add_diagnostic_variable!(
    short_name = "stab",
    long_name = "Static Stability: N^2 = -g/theta dtheta/dz",
    standard_name = "static_stability",
    units = "s^-2",
    compute! = (out, state, cache, time) -> begin
        thermo_params = CAP.thermodynamics_params(cache.params)
        ᶜts = cache.precomputed.ᶜts
        θ_virt = TD.virtual_pottemp.(thermo_params, ᶜts)
        dθ_virtdz = @. ᶜgradᵥ_.(ᶠinterp_(θ_virt))
        temp = cache.scratch.ᶜtemp_scalar
        parent(temp) .= parent(CAP.grav(cache.params) ./ θ_virt .* dθ_virtdz)
        if isnothing(out)
            return copy(temp)
        else
            out .= temp
        end
    end,
)

CAD.add_diagnostic_variable!(
    short_name = "vT",
    long_name = "Product of meridional wind and temperature",
    standard_name = "vT",
    units = "K m/s",
    compute! = (out, state, cache, time) -> begin
        thermo_params = CAP.thermodynamics_params(cache.params)
        T = TD.air_temperature.(thermo_params, cache.precomputed.ᶜts)
        v = ClimaCore.Geometry.UVVector.(state.c.uₕ).components.data.:2
        if isnothing(out)
            return  v .* T
        else
            out .= v .* T
        end
    end,
)

const ᶜgradᵥ_ = ClimaCore.Operators.GradientF2C()
const ᶠinterp_ = ClimaCore.Operators.InterpolateC2F(
    bottom = ClimaCore.Operators.Extrapolate(),
    top = ClimaCore.Operators.Extrapolate(),
)

CAD.add_diagnostic_variable!(
    short_name = "egr",
    long_name = "max. Eady growth rate (0.31 f/N du/dz)",
    standard_name = "Eady_growth_rate",
    units = "s^-1",
    compute! = (out, state, cache, time) -> begin
        FT = eltype(state)
        thermo_params = CAP.thermodynamics_params(cache.params)
        ᶜts = cache.precomputed.ᶜts
        θ_virt = TD.virtual_pottemp.(thermo_params, ᶜts)
        dθ_virtdz = @. ᶜgradᵥ_.(ᶠinterp_(θ_virt))
        N = cache.scratch.ᶜtemp_scalar
        Nsq = CAP.grav(cache.params) ./ θ_virt .* ClimaCore.Geometry.WVector.(dθ_virtdz)
        mask = parent(Nsq) .> 0
        parent(N) .= sqrt.((parent(Nsq) .* mask ))

        Ω = CAP.Omega(cache.params)
        φ = ClimaCore.Fields.coordinate_field(N).lat
        f = 2Ω .* sin.(φ)
        u = ClimaCore.Geometry.UVVector.(state.c.uₕ).components.data.:1
        grad_u = @. ᶜgradᵥ_.(ᶠinterp_(u))
        du_dz = ClimaCore.Geometry.WVector.(grad_u)
        temp = cache.scratch.ᶜtemp_scalar
        parent(temp) .= parent(FT(0.31) .* f ./ N .* du_dz)
        if isnothing(out)
            return copy(temp)
        else
            out .= temp
        end
    end,
)

# climate diagnostics: ta, ua, (hus), *mass_streamfunction, *stab, mse

# storm track diagnostics: [vT], [v][T], egr = f/N dTdy

# vars = [mse, lr, ediff, mass_streamfunction, stab, vT, egr]