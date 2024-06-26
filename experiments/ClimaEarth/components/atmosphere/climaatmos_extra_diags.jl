# these extensions add extra diagnostics to the atmos model output
import ClimaAtmos.Diagnostics as CAD

const ᶜgradᵥ_ = CC.Operators.GradientF2C()
const ᶠinterp_ = CC.Operators.InterpolateC2F(bottom = CC.Operators.Extrapolate(), top = CC.Operators.Extrapolate())

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
    units = "J kg^-1",
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
    units = "K m^-1",
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
    units = "m^2 s^-1",
    comments = "Eddy diffusivity consistent with the VerticalDiffusion scheme in ClimaAtmos.",
    compute! = (out, state, cache, time) -> begin
        (; ᶜp) = cache.precomputed
        (; C_E) = cache.atmos.vert_diff
        interior_uₕ = CC.Fields.level(state.c.uₕ, 1)
        ᶠp = ᶠK_E = cache.scratch.ᶠtemp_scalar
        @. ᶠp = CAD.ᶠinterp(ᶜp)
        ᶜΔz_surface = CC.Fields.Δz_field(interior_uₕ)
        @. ᶠK_E = CA.eddy_diffusivity_coefficient(C_E, CA.norm(interior_uₕ), ᶜΔz_surface / 2, ᶠp)
        if isnothing(out)
            return CAD.ᶜinterp.(ᶠK_E)
        else
            out .= CAD.ᶜinterp.(ᶠK_E)
        end

    end,
)

CAD.add_diagnostic_variable!(
    short_name = "mass_strf",
    long_name = "Meridional Mass Streamfunction",
    standard_name = "meridional_mass_streamfunction",
    units = "kg m^-1 s^-1",
    comments = "Calculated as the vertical integral of the meridional mass flux: 2πa cos(φ) ∫_0^z ρv dz, with a the Earth's radius, φ the latitude, and ρv the meridional mass flux.",
    compute! = (out, state, cache, time) -> begin
        FT = eltype(state)
        φ = CC.Fields.coordinate_field(axes(state.f)).lat
        ρv = CC.Geometry.UVVector.(state.c.uₕ).components.data.:2 .* state.c.ρ
        ᶠ∫_0_z_ρv = cache.scratch.ᶠtemp_scalar
        CC.Operators.column_integral_indefinite!(ᶠ∫_0_z_ρv, ρv)
        twoπa = FT(2π * 6371e3)
        if isnothing(out)
            return ᶠ∫_0_z_ρv .* twoπa .* cosd.(φ)
        else
            out .= ᶠ∫_0_z_ρv .* twoπa .* cosd.(φ)
        end
    end,
)

CAD.add_diagnostic_variable!(
    short_name = "stab",
    long_name = "Static Stability",
    standard_name = "static_stability",
    units = "s^-2",
    comments = "Static stability, N^2 = g/θᵥ ∂θᵥ∂z, where g is the gravity, θ is the potential temperature, and dθ/dz is the vertical gradient of the potential temperature. It is a measure of the vertical thermal stratification of the atmosphere.",
    compute! = (out, state, cache, time) -> begin
        if isnothing(out)
            return static_stability(cache)
        else
            out .= static_stability(cache)
        end
    end,
)

CAD.add_diagnostic_variable!(
    short_name = "vt",
    long_name = "Meridional temperature flux",
    standard_name = "vt",
    units = "K m s^-1",
    comments = "Product of the meridional wind and the air temperature.",
    compute! = (out, state, cache, time) -> begin
        thermo_params = CAP.thermodynamics_params(cache.params)
        T = TD.air_temperature.(thermo_params, cache.precomputed.ᶜts)
        v = CC.Geometry.UVVector.(state.c.uₕ).components.data.:2
        if isnothing(out)
            return v .* T
        else
            out .= v .* T
        end
    end,
)

CAD.add_diagnostic_variable!(
    short_name = "egr",
    long_name = "(zonal) max. Eady growth rate = 0.31 f/N du/dz",
    standard_name = "max_eady_growth_rate",
    units = "s^-1",
    comments = "The maximum Eady growth rate, which is the product of the Coriolis parameter f = 2Ω sin(φ), the Brunt-Väisälä frequency N, and the vertical gradient of the zonal wind du/dz. It is a theoretical measure of the potential for baroclinic instability.",
    compute! = (out, state, cache, time) -> begin
        FT = eltype(state)
        Ω = CAP.Omega(cache.params)
        φ = CC.Fields.coordinate_field(state.c).lat
        # calculate du/dz
        u = CC.Geometry.UVVector.(state.c.uₕ).components.data.:1
        du_dz = cache.scratch.ᶜtemp_scalar
        @. du_dz = CC.Geometry.WVector.(ᶜgradᵥ_.(ᶠinterp_(u))).components.data.:1
        # calculate Brunt-Väisälä frequency N from static stability, mask out imaginary values
        N² = cache.scratch.ᶜtemp_scalar_2
        N² .= static_stability(cache)
        N = cache.scratch.ᶜtemp_scalar_3
        @. N = sqrt(max(N², 0))
        # calculate max. Eady growth rate
        if isnothing(out)
            return FT(0.31) .* 2Ω .* sind.(φ) ./ N .* du_dz
        else
            out .= FT(0.31) .* 2Ω .* sind.(φ) ./ N .* du_dz
        end
    end,
)

"""
    static_stability(cache)

This function computes the static stability, or the square of the Brunt-Väisälä
frequency N² = g/θᵥ ∂θᵥ∂z, where g is the gravity, θᵥ is the virtual potential temperature and ∂θᵥ∂z is its vertical gradient.
"""
function static_stability(cache)
    thermo_params = CAP.thermodynamics_params(cache.params)
    ᶜts = cache.precomputed.ᶜts
    θᵥ = TD.virtual_pottemp.(thermo_params, ᶜts)
    ∂θᵥ∂z = cache.scratch.ᶜtemp_scalar
    @. ∂θᵥ∂z = CC.Geometry.WVector.(ᶜgradᵥ_.(ᶠinterp_(θᵥ))).components.data.:1
    return CAP.grav(cache.params) ./ θᵥ .* ∂θᵥ∂z
end
