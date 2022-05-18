# Atmos boundary condition for a coupled simulation, which calculates and accumulates the boundary flux
using ClimaCore.Utilities: PlusHalf

struct CoupledFlux <: BCtag end
function bc_divF2C_bottom!(::CoupledFlux, dY, Y, p, t)
    # flux calculation
    coords = Fields.coordinate_field(p.domain.hv_center_space)

    z = coords.z
    Yc = Y.Yc

    uₕ = Yc.ρuₕ ./ Yc.ρ
    ρw = Y.ρw
    If2c = Operators.InterpolateF2C()
    Ic2f = Operators.InterpolateC2F(bottom = Operators.Extrapolate(), top = Operators.Extrapolate())
    w = If2c.(ρw) ./ Yc.ρ
    cuv = @. Geometry.UWVector(uₕ)
    windspeed = @. norm(cuv)
    windspeed_f = @. Ic2f(windspeed)
    windspeed_boundary = Fields.level(windspeed_f, PlusHalf(1))
    θ = Yc.ρθ ./ Yc.ρ
    θ_boundary = Fields.level(Ic2f.(θ), PlusHalf(1))
    ρ_f = @. Ic2f(Yc.ρ)
    ρ_boundary = Fields.level(ρ_f, PlusHalf(1))

    # build atmos face fields on surface  boundary space to enable broadcasting
    windspeed_boundary = Fields.Field(Fields.field_values(windspeed_boundary), axes(p.T_sfc))
    θ_boundary = Fields.Field(Fields.field_values(θ_boundary), axes(p.T_sfc))
    ρ_boundary = Fields.Field(Fields.field_values(ρ_boundary), axes(p.T_sfc))

    λ = @. p.cpl_p.C_p * p.cpl_p.C_H * ρ_boundary * windspeed_boundary
    dθ = @. θ_boundary - p.T_sfc
    heat_flux = @. -λ * dθ
    @. dY.F_sfc += heat_flux # accumulation

    return Operators.SetValue(Geometry.WVector.(heat_flux))
end
