#=
# Per-polygon turbulent fluxes on the exchange grid

One `SurfaceFluxes.surface_fluxes` evaluation per exchange-grid polygon, from
atmospheric state gathered off the SE nodes and surface state read from the
owning FV cell, aggregated conservatively to both sides:

- to the FV grid with per-polygon sea-ice-concentration weighting;
- to the SE boundary space via the L2 scatter, `weighted_dss!`, then division
  by the DSS'd nodal wet coverage — dividing only after DSS avoids averaging
  real coastal fluxes with no-data zeros from dry neighbor elements.

Momentum is handled in the global UV (east/north) basis — whose components
are DSS-safe scalars — and converted to the local CT1/CT2 components expected
in `csf.F_turb_ρτxz/yz` only at the very end. All per-step operations run on
the compute device.
=#

import StaticArrays
import NVTX

"""
    ExchangeFluxState{FT, VF}

Device-resident per-polygon scratch for the exchange-grid flux computation of
one surface model (one instance for the ocean, one for sea ice): the gathered
atmospheric and surface state (momentum in the UV basis; `T_sfc` in [K]), the
flux outputs (`F_*`), running time accumulators for the slow-surface path
(`acc_*`, with `n_acc` counting contributions), and two generic scratch
vectors.
"""
struct ExchangeFluxState{FT, VF <: AbstractVector{FT}}
    T_atmos::VF
    q_tot::VF
    q_liq::VF
    q_ice::VF
    ρ_atmos::VF
    u_atmos::VF
    v_atmos::VF
    height_int::VF
    height_sfc::VF
    T_sfc::VF
    sic::VF
    F_sh::VF
    F_lh::VF
    F_moisture::VF
    F_τu::VF
    F_τv::VF
    acc_F_sh::VF
    acc_F_lh::VF
    acc_F_moisture::VF
    acc_F_τu::VF
    acc_F_τv::VF
    scratch1::VF
    scratch2::VF
    n_acc::Base.RefValue{Int}
end

Adapt.@adapt_structure ExchangeFluxState

function ExchangeFluxState{FT}(arch, n_poly::Int) where {FT}
    state = ExchangeFluxState{FT, Vector{FT}}(
        ntuple(_ -> zeros(FT, n_poly), 23)...,
        Ref(0),
    )
    return on_device(arch, state)
end

"""
    IceExchangeState{FT, VF}

Per-polygon scratch for the sea-ice exchange-grid flux computation: the
common [`ExchangeFluxState`](@ref) plus the ice-specific inputs of the
skin-temperature flux balance (`R` — conductive resistance, `T_i` — ice
internal/interface temperature [K], downwelling `SW_d`/`LW_d`) and the
diagnosed surface temperature output `T_sfc_new` [K].
"""
struct IceExchangeState{FT, VF <: AbstractVector{FT}}
    fluxes::ExchangeFluxState{FT, VF}
    R::VF
    T_i::VF
    SW_d::VF
    LW_d::VF
    T_sfc_new::VF
end

Adapt.@adapt_structure IceExchangeState

function IceExchangeState{FT}(arch, n_poly::Int) where {FT}
    fluxes = ExchangeFluxState{FT, Vector{FT}}(
        ntuple(_ -> zeros(FT, n_poly), 23)...,
        Ref(0),
    )
    state = IceExchangeState(fluxes, ntuple(_ -> zeros(FT, n_poly), 5)...)
    return on_device(arch, state)
end

"""
    average_and_reset_exchange_accumulators!(fs::ExchangeFluxState) -> Bool

Write the time-averaged accumulated fluxes into the flux-state outputs
(`F_*`), zero the accumulators and the contribution counter. Returns `false`
without touching anything when no contributions have been accumulated.
"""
function average_and_reset_exchange_accumulators!(fs::ExchangeFluxState)
    n = fs.n_acc[]
    iszero(n) && return false
    @. fs.F_sh = fs.acc_F_sh / n
    @. fs.F_lh = fs.acc_F_lh / n
    @. fs.F_moisture = fs.acc_F_moisture / n
    @. fs.F_τu = fs.acc_F_τu / n
    @. fs.F_τv = fs.acc_F_τv / n
    fill!(fs.acc_F_sh, 0)
    fill!(fs.acc_F_lh, 0)
    fill!(fs.acc_F_moisture, 0)
    fill!(fs.acc_F_τu, 0)
    fill!(fs.acc_F_τv, 0)
    fs.n_acc[] = 0
    return true
end

"""
    gather_atmos_state_to_polys!(fs::ExchangeFluxState, eg::ExchangeGrid, csf,
                                 temp_uv_vec)

Gather the atmospheric near-surface state from the coupler fields (SE nodal)
onto the exchange-grid polygons. `csf.u_int`/`csf.v_int` are converted from
the local CT1/CT2 basis to the global UV basis before gathering, so
per-polygon winds and stresses are basis-consistent across the nodes a
polygon touches.
"""
NVTX.@annotate function gather_atmos_state_to_polys!(
    fs::ExchangeFluxState,
    eg::ExchangeGrid,
    csf,
    temp_uv_vec,
)
    CRExt = get_ConservativeRegriddingCCExt()
    contravariant_to_cartesian!(temp_uv_vec, csf.u_int, csf.v_int)
    u_uv = temp_uv_vec.components.data.:1
    v_uv = temp_uv_vec.components.data.:2
    gather_nodes_to_polys!(fs.u_atmos, eg, CRExt.se_field_to_vec(u_uv))
    gather_nodes_to_polys!(fs.v_atmos, eg, CRExt.se_field_to_vec(v_uv))
    gather_nodes_to_polys!(fs.T_atmos, eg, CRExt.se_field_to_vec(csf.T_atmos))
    gather_nodes_to_polys!(fs.q_tot, eg, CRExt.se_field_to_vec(csf.q_tot_atmos))
    gather_nodes_to_polys!(fs.q_liq, eg, CRExt.se_field_to_vec(csf.q_liq_atmos))
    gather_nodes_to_polys!(fs.q_ice, eg, CRExt.se_field_to_vec(csf.q_ice_atmos))
    gather_nodes_to_polys!(fs.ρ_atmos, eg, CRExt.se_field_to_vec(csf.ρ_atmos))
    gather_nodes_to_polys!(fs.height_int, eg, CRExt.se_field_to_vec(csf.height_int))
    gather_nodes_to_polys!(fs.height_sfc, eg, CRExt.se_field_to_vec(csf.height_sfc))
    return nothing
end

# SurfaceFluxes evaluation for one polygon. Mirrors the nodal computation in
# `FluxCalculator.compute_surface_fluxes!` (surface density and saturation
# humidity from the atmospheric state, zero surface velocity, zero
# displacement height), but with scalar inputs so it can run inside a kernel.
@inline function _polygon_surface_fluxes(
    surface_fluxes_params,
    thermo_params,
    config,
    T_atmos,
    q_tot,
    q_liq,
    q_ice,
    ρ_atmos,
    u_atmos,
    v_atmos,
    height_int,
    height_sfc,
    T_sfc,
)
    FT = typeof(T_atmos)
    ρ_sfc = SF.surface_density(
        surface_fluxes_params,
        T_atmos,
        ρ_atmos,
        T_sfc,
        height_int - height_sfc,
        q_tot,
        FT(0), # q_liq
        FT(0), # q_ice
    )
    q_sfc = TD.q_vap_saturation(thermo_params, T_sfc, ρ_sfc, FT(0), FT(0))
    uv_int = StaticArrays.SVector(u_atmos, v_atmos)
    return FluxCalculator.get_surface_fluxes(
        surface_fluxes_params,
        uv_int,
        T_atmos,
        q_tot,
        q_liq,
        q_ice,
        ρ_atmos,
        height_int,
        zero(uv_int), # u_sfc
        T_sfc,
        q_sfc,
        height_sfc,
        FT(0), # d
        config,
    )
end

# The flux-state vectors as a NamedTuple, for passing to GPU kernels (the
# struct itself holds a host `RefValue` and cannot be a kernel argument).
@inline _kernel_state(fs::ExchangeFluxState) = (;
    fs.T_atmos,
    fs.q_tot,
    fs.q_liq,
    fs.q_ice,
    fs.ρ_atmos,
    fs.u_atmos,
    fs.v_atmos,
    fs.height_int,
    fs.height_sfc,
    fs.T_sfc,
    fs.sic,
    fs.F_sh,
    fs.F_lh,
    fs.F_moisture,
    fs.F_τu,
    fs.F_τv,
    fs.acc_F_sh,
    fs.acc_F_lh,
    fs.acc_F_moisture,
    fs.acc_F_τu,
    fs.acc_F_τv,
)

@inline function _store_ocean_polygon_fluxes!(fs, k, out)
    @inbounds begin
        fs.F_sh[k] = out.F_sh
        fs.F_lh[k] = out.F_lh
        fs.F_moisture[k] = out.F_turb_moisture
        fs.F_τu[k] = out.F_turb_ρτxz
        fs.F_τv[k] = out.F_turb_ρτyz
        fs.acc_F_sh[k] += out.F_sh
        fs.acc_F_lh[k] += out.F_lh
        fs.acc_F_moisture[k] += out.F_turb_moisture
        fs.acc_F_τu[k] += out.F_turb_ρτxz
        fs.acc_F_τv[k] += out.F_turb_ρτyz
    end
    return nothing
end

@kernel function _ocean_polygon_fluxes_kernel!(
    fs,
    surface_fluxes_params,
    thermo_params,
    config,
)
    k = @index(Global)
    @inbounds out = _polygon_surface_fluxes(
        surface_fluxes_params,
        thermo_params,
        config,
        fs.T_atmos[k],
        fs.q_tot[k],
        fs.q_liq[k],
        fs.q_ice[k],
        fs.ρ_atmos[k],
        fs.u_atmos[k],
        fs.v_atmos[k],
        fs.height_int[k],
        fs.height_sfc[k],
        fs.T_sfc[k],
    )
    _store_ocean_polygon_fluxes!(fs, k, out)
end

"""
    compute_ocean_polygon_fluxes!(fs::ExchangeFluxState, surface_fluxes_params,
                                  thermo_params, config)

Run the SurfaceFluxes evaluation for every exchange-grid polygon, storing the
flux outputs and adding them to the running time accumulators (fused; call
sites increment `fs.n_acc`).
"""
NVTX.@annotate function compute_ocean_polygon_fluxes!(
    fs::ExchangeFluxState,
    surface_fluxes_params,
    thermo_params,
    config,
)
    if _is_cpu(fs.F_sh)
        for k in eachindex(fs.F_sh)
            @inbounds out = _polygon_surface_fluxes(
                surface_fluxes_params,
                thermo_params,
                config,
                fs.T_atmos[k],
                fs.q_tot[k],
                fs.q_liq[k],
                fs.q_ice[k],
                fs.ρ_atmos[k],
                fs.u_atmos[k],
                fs.v_atmos[k],
                fs.height_int[k],
                fs.height_sfc[k],
                fs.T_sfc[k],
            )
            _store_ocean_polygon_fluxes!(fs, k, out)
        end
    else
        backend = KernelAbstractions.get_backend(fs.F_sh)
        _ocean_polygon_fluxes_kernel!(backend, _KA_WORKGROUP)(
            _kernel_state(fs),
            surface_fluxes_params,
            thermo_params,
            config;
            ndrange = length(fs.F_sh),
        )
    end
    return nothing
end

@inline _kernel_state(is::IceExchangeState) =
    merge(_kernel_state(is.fluxes), (; is.R, is.T_i, is.SW_d, is.LW_d, is.T_sfc_new))

# SurfaceFluxes evaluation with skin-temperature diagnosis for one sea-ice
# polygon. Polygons without ice short-circuit to zero flux. Mirrors the nodal
# computation in the boundary-space `compute_surface_fluxes!` for
# `ClimaSeaIceSimulation`.
@inline function _polygon_ice_fluxes(
    surface_fluxes_params,
    thermo_params,
    config,
    σ,
    ϵ,
    α_albedo,
    T_melt,
    T_atmos,
    q_tot,
    q_liq,
    q_ice,
    ρ_atmos,
    u_atmos,
    v_atmos,
    height_int,
    height_sfc,
    T_sfc_guess,
    sic,
    R,
    T_i,
    SW_d,
    LW_d,
)
    FT = typeof(T_atmos)
    if sic <= 0
        z = FT(0)
        return (;
            F_turb_ρτxz = z,
            F_turb_ρτyz = z,
            F_sh = z,
            F_lh = z,
            F_turb_moisture = z,
            T_sfc_new = T_sfc_guess,
        )
    end
    update_T_sfc_cb = update_T_sfc(R, T_i, σ, ϵ, SW_d, LW_d, α_albedo, T_melt)
    ρ_sfc = SF.surface_density(
        surface_fluxes_params,
        T_atmos,
        ρ_atmos,
        T_sfc_guess,
        height_int - height_sfc,
        q_tot,
        FT(0),
        FT(0),
    )
    q_sfc = TD.q_vap_saturation(thermo_params, T_sfc_guess, ρ_sfc, FT(0), FT(0))
    uv_int = StaticArrays.SVector(u_atmos, v_atmos)
    return FluxCalculator.get_surface_fluxes(
        surface_fluxes_params,
        uv_int,
        T_atmos,
        q_tot,
        q_liq,
        q_ice,
        ρ_atmos,
        height_int,
        zero(uv_int),
        T_sfc_guess,
        q_sfc,
        height_sfc,
        FT(0), # d
        config,
        update_T_sfc_cb,
    )
end

@inline function _store_ice_polygon_fluxes!(vecs, k, out)
    _store_ocean_polygon_fluxes!(vecs, k, out)
    @inbounds vecs.T_sfc_new[k] = out.T_sfc_new
    return nothing
end

@inline function _ice_polygon_fluxes_at(
    vecs,
    k,
    surface_fluxes_params,
    thermo_params,
    config,
    σ,
    ϵ,
    α_albedo,
    T_melt,
)
    @inbounds out = _polygon_ice_fluxes(
        surface_fluxes_params,
        thermo_params,
        config,
        σ,
        ϵ,
        α_albedo,
        T_melt,
        vecs.T_atmos[k],
        vecs.q_tot[k],
        vecs.q_liq[k],
        vecs.q_ice[k],
        vecs.ρ_atmos[k],
        vecs.u_atmos[k],
        vecs.v_atmos[k],
        vecs.height_int[k],
        vecs.height_sfc[k],
        vecs.T_sfc[k],
        vecs.sic[k],
        vecs.R[k],
        vecs.T_i[k],
        vecs.SW_d[k],
        vecs.LW_d[k],
    )
    _store_ice_polygon_fluxes!(vecs, k, out)
    return nothing
end

@kernel function _ice_polygon_fluxes_kernel!(
    vecs,
    surface_fluxes_params,
    thermo_params,
    config,
    σ,
    ϵ,
    α_albedo,
    T_melt,
)
    k = @index(Global)
    _ice_polygon_fluxes_at(
        vecs,
        k,
        surface_fluxes_params,
        thermo_params,
        config,
        σ,
        ϵ,
        α_albedo,
        T_melt,
    )
end

"""
    compute_ice_polygon_fluxes!(is::IceExchangeState, surface_fluxes_params,
                                thermo_params, config, σ, ϵ, α_albedo, T_melt)

Run the SurfaceFluxes evaluation with skin-temperature diagnosis for every
exchange-grid polygon with ice (`sic > 0`), storing the flux outputs, the
diagnosed `T_sfc_new`, and adding the fluxes to the running accumulators.
"""
NVTX.@annotate function compute_ice_polygon_fluxes!(
    is::IceExchangeState,
    surface_fluxes_params,
    thermo_params,
    config,
    σ,
    ϵ,
    α_albedo,
    T_melt,
)
    vecs = _kernel_state(is)
    if _is_cpu(is.fluxes.F_sh)
        for k in eachindex(is.fluxes.F_sh)
            _ice_polygon_fluxes_at(
                vecs,
                k,
                surface_fluxes_params,
                thermo_params,
                config,
                σ,
                ϵ,
                α_albedo,
                T_melt,
            )
        end
    else
        backend = KernelAbstractions.get_backend(is.fluxes.F_sh)
        _ice_polygon_fluxes_kernel!(backend, _KA_WORKGROUP)(
            vecs,
            surface_fluxes_params,
            thermo_params,
            config,
            σ,
            ϵ,
            α_albedo,
            T_melt;
            ndrange = length(is.fluxes.F_sh),
        )
    end
    return nothing
end

"""
    cartesian_to_contravariant!(ρτxz, ρτyz, uv_field)

Convert a `UVVector` field into the local CT1/CT2 unit-basis components
expected by the coupler flux fields (`csf.F_turb_ρτxz/yz`), via
`CA.projected_vector_data` — the same projection that defines `csf.u_int` in
the ClimaAtmos extension. Exact inverse of `contravariant_to_cartesian!`.
"""
function cartesian_to_contravariant!(ρτxz, ρτyz, uv_field)
    local_geometry = CC.Fields.local_geometry_field(ρτxz)
    @. ρτxz = CA.projected_vector_data(CT1, uv_field, local_geometry)
    @. ρτyz = CA.projected_vector_data(CT2, uv_field, local_geometry)
    return nothing
end

"""
    scatter_poly_fluxes_to_boundary!(remapping, eg::ExchangeGrid,
                                     fs::ExchangeFluxState, weight;
                                     cov_cutoff = 1e-3)

Aggregate per-polygon fluxes onto the SE boundary space as a `weight`-weighted
average, filling the `remapping.flux_scratch` fields for
`FluxCalculator.update_flux_fields!`. `weight` selects the sub-surface the
fluxes apply to — `1 - sic` for open ocean, `sic` for sea ice — so the nodal
result is a per-unit-*weighted*-area flux, consistent with the area fraction
the coupler multiplies it by: L2-scatter `weight * F` (momentum in UV) and
`weight` itself, `weighted_dss!` each scalar, divide by the DSS'd weighted
coverage, then convert momentum UV → CT1/CT2. Nodes with relative coverage
below `cov_cutoff` get zero flux (they are essentially not covered by wet
ocean; their area fraction vanishes there too, so they never contribute to
the coupler sums). Uses `fs.scratch1` internally; `weight` must not alias it.
"""
NVTX.@annotate function scatter_poly_fluxes_to_boundary!(
    remapping,
    eg::ExchangeGrid,
    fs::ExchangeFluxState,
    weight;
    cov_cutoff = 1e-3,
)
    CRExt = get_ConservativeRegriddingCCExt()
    fx = remapping.flux_scratch
    cov = remapping.weight_cov_scratch
    scatter_polys_to_nodes!(CRExt.se_field_to_vec(cov), eg, weight)

    @. fs.scratch1 = weight * fs.F_sh
    scatter_polys_to_nodes!(CRExt.se_field_to_vec(fx.F_sh), eg, fs.scratch1)
    @. fs.scratch1 = weight * fs.F_lh
    scatter_polys_to_nodes!(CRExt.se_field_to_vec(fx.F_lh), eg, fs.scratch1)
    @. fs.scratch1 = weight * fs.F_moisture
    scatter_polys_to_nodes!(CRExt.se_field_to_vec(fx.F_turb_moisture), eg, fs.scratch1)
    @. fs.scratch1 = weight * fs.F_τu
    scatter_polys_to_nodes!(CRExt.se_field_to_vec(fx.F_turb_ρτxz), eg, fs.scratch1)
    @. fs.scratch1 = weight * fs.F_τv
    scatter_polys_to_nodes!(CRExt.se_field_to_vec(fx.F_turb_ρτyz), eg, fs.scratch1)

    Utilities.apply_dss!(cov, remapping.flux_dss_buffer)
    FT = CC.Spaces.undertype(axes(cov))
    cutoff = FT(cov_cutoff)
    for field in values(fx)
        Utilities.apply_dss!(field, remapping.flux_dss_buffer)
        @. field = ifelse(cov > cutoff, field / max(cov, cutoff), FT(0))
    end

    @. remapping.temp_uv_vec = CC.Geometry.UVVector(fx.F_turb_ρτxz, fx.F_turb_ρτyz)
    cartesian_to_contravariant!(fx.F_turb_ρτxz, fx.F_turb_ρτyz, remapping.temp_uv_vec)
    return nothing
end
