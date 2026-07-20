#=
# Per-polygon turbulent fluxes on the exchange grid

Turbulent fluxes for the ocean (and sea ice) are computed with one
`SurfaceFluxes.surface_fluxes` evaluation per exchange-grid polygon, from
atmospheric state gathered off the SE nodes and surface state read directly
from the owning FV cell. The per-polygon fluxes are then aggregated
conservatively to both sides:

- to the FV grid with per-polygon sea-ice-concentration weighting (the
  coastline- and ice-edge-resolving improvement over remapping
  boundary-space fluxes);
- to the SE boundary space via the L2 scatter followed by `weighted_dss!`.
  The scatter result is coverage-weighted, so it is divided by the DSS'd
  nodal wet coverage (`node_cov_dss`), which — unlike normalizing before
  DSS — does not average real coastal fluxes with no-data zeros from dry
  neighbor elements. Momentum is handled in the UV (east/north) basis
  through gather, kernel, scatter and DSS (the UV basis is global, so its
  components are DSS-safe scalars), and converted to the local CT1/CT2
  components expected in `csf.F_turb_ρτxz/yz` only at the very end.

All per-step operations run on the compute device; the only per-step host
work is kernel launches and small array-view wrappers.
=#

import StaticArrays
import NVTX

# Relative nodal-coverage cutoff below which the SE-side flux is set to zero
# (nodes essentially not covered by wet ocean; their area fraction vanishes
# there too, so they never contribute to the coupler sums).
const EXCHANGE_COV_CUTOFF = 1e-3

"""
    ExchangeFluxState{FT, VF}

Device-resident per-polygon scratch for the exchange-grid flux computation of
one surface model (one instance for the ocean, one for sea ice).

Holds the gathered atmospheric state (`T_atmos`, `q_tot`, `q_liq`, `q_ice`,
`ρ_atmos`, `u_atmos`, `v_atmos` in the UV basis, `height_int`, `height_sfc`),
the gathered surface state (`T_sfc` [K], `sic`), the flux outputs (`F_sh`,
`F_lh`, `F_moisture`, `F_τu`, `F_τv` in the UV basis), running time
accumulators for the slow-surface path (`acc_*`, with `n_acc` counting
contributions), and two generic per-polygon scratch vectors.
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
    DEBUG_POLY_FLUX_NANS

Opt-in diagnostic toggle for [`check_poly_flux_nans`](@ref). Off by default
(the check costs one host `Bool` read per flux evaluation). Enable from the
REPL with

    Base.get_extension(ClimaCoupler, :ClimaCouplerCMIPExt).DEBUG_POLY_FLUX_NANS[] = true
"""
const DEBUG_POLY_FLUX_NANS = Ref(false)

_poly_flux_outputs(fs::ExchangeFluxState) =
    (; fs.F_sh, fs.F_lh, fs.F_moisture, fs.F_τu, fs.F_τv)
_poly_flux_outputs(is::IceExchangeState) =
    merge(_poly_flux_outputs(is.fluxes), (; is.T_sfc_new))

_poly_flux_inputs(fs::ExchangeFluxState) = (;
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
)
_poly_flux_inputs(is::IceExchangeState) =
    merge(_poly_flux_inputs(is.fluxes), (; is.R, is.T_i, is.SW_d, is.LW_d))

"""
    check_poly_flux_nans(state, eg::ExchangeGrid, label)

Purely diagnostic NaN guard for the per-polygon flux evaluation (no-op unless
[`DEBUG_POLY_FLUX_NANS`](@ref) is enabled; never modifies any values). Scans
the flux outputs of `state` (an [`ExchangeFluxState`](@ref) or
[`IceExchangeState`](@ref)); if any is NaN, logs — for up to 10 offending
polygons — the polygon index, its owning FV cell and SE element, all flux
outputs, and all gathered inputs of the SurfaceFluxes evaluation, then
`error`s so the NaN is never scattered to the models. The detailed report
copies the per-polygon vectors to the host, which is fine on this abort path.
"""
function check_poly_flux_nans(state, eg::ExchangeGrid, label)
    DEBUG_POLY_FLUX_NANS[] || return nothing
    outputs = _poly_flux_outputs(state)
    any(v -> any(isnan, v), values(outputs)) || return nothing

    outs = map(Array, outputs)
    inputs = map(Array, _poly_flux_inputs(state))
    oc_of_poly = Array(eg.oc_of_poly)
    elem_of_poly = Array(eg.elem_of_poly)
    bad = findall(k -> any(o -> isnan(o[k]), values(outs)), 1:eg.n_poly)
    for k in first(bad, 10)
        out_str = join(("$nm = $(o[k])" for (nm, o) in pairs(outs)), ", ")
        in_str = join(("$nm = $(v[k])" for (nm, v) in pairs(inputs)), ", ")
        @error "NaN $label fluxes at polygon $k (FV cell $(oc_of_poly[k]), SE elem $(elem_of_poly[k])): $out_str; inputs: $in_str"
    end
    return error(
        "NaN in $label exchange-grid polygon fluxes at $(length(bad)) polygons" *
        (length(bad) > 10 ? " (first 10 reported above)" : ""),
    )
end

"""
    momentum_basis_fields(boundary_space)

Precompute the per-node 2×2 map between the local CT1/CT2 unit basis and the
global UV (east/north) basis: `a = UV(ĈT1)`, `b = UV(ĈT2)` and its
determinant. Enables allocation-free per-step basis conversions
([`ct12_to_uv!`](@ref), [`uv_to_ct12!`](@ref)); numerically identical to
`contravariant_to_cartesian!` and its inverse.
"""
function momentum_basis_fields(boundary_space)
    FT = CC.Spaces.undertype(boundary_space)
    u_ct = CC.Fields.ones(boundary_space)
    v_ct = CC.Fields.zeros(boundary_space)
    uv = CC.Fields.Field(CC.Geometry.UVVector{FT}, boundary_space)
    contravariant_to_cartesian!(uv, u_ct, v_ct)
    au = copy(uv.components.data.:1)
    av = copy(uv.components.data.:2)
    u_ct .= FT(0)
    v_ct .= FT(1)
    contravariant_to_cartesian!(uv, u_ct, v_ct)
    bu = copy(uv.components.data.:1)
    bv = copy(uv.components.data.:2)
    det = @. au * bv - av * bu
    return (; au, av, bu, bv, det)
end

"""
    ct12_to_uv!(uv_field, u_ct, v_ct, basis)

Allocation-free equivalent of `contravariant_to_cartesian!` using the
precomputed [`momentum_basis_fields`](@ref).
"""
function ct12_to_uv!(uv_field, u_ct, v_ct, basis)
    @. uv_field = CC.Geometry.UVVector(
        u_ct * basis.au + v_ct * basis.bu,
        u_ct * basis.av + v_ct * basis.bv,
    )
    return nothing
end

"""
    uv_to_ct12!(u_ct, v_ct, uv_field, basis)

Allocation-free equivalent of [`cartesian_to_contravariant!`](@ref) using the
precomputed [`momentum_basis_fields`](@ref).
"""
function uv_to_ct12!(u_ct, v_ct, uv_field, basis)
    u_uv = uv_field.components.data.:1
    v_uv = uv_field.components.data.:2
    @. u_ct = (u_uv * basis.bv - v_uv * basis.bu) / basis.det
    @. v_ct = (basis.au * v_uv - basis.av * u_uv) / basis.det
    return nothing
end

"""
    gather_atmos_state_to_polys!(fs::ExchangeFluxState, eg::ExchangeGrid, csf,
                                 temp_uv_vec, momentum_basis)

Gather the atmospheric near-surface state from the coupler fields (SE nodal)
onto the exchange-grid polygons. The nodal velocity components `csf.u_int` /
`csf.v_int` are given in the local CT1/CT2 unit basis; they are converted to
the global UV (east/north) basis (via `temp_uv_vec` scratch and the
precomputed `momentum_basis`) before gathering, so that per-polygon wind
speeds and stresses are basis-consistent across the nodes a polygon touches.
"""
NVTX.@annotate function gather_atmos_state_to_polys!(
    fs::ExchangeFluxState,
    eg::ExchangeGrid,
    csf,
    temp_uv_vec,
    momentum_basis,
)
    CRExt = get_ConservativeRegriddingCCExt()
    ct12_to_uv!(temp_uv_vec, csf.u_int, csf.v_int, momentum_basis)
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

# Inverse of `contravariant_to_cartesian!`: local CT1/CT2 unit-basis
# components of a UVVector, by exactly inverting the 2×2 forward map
# `uv = τ1 * UV(ĈT1) + τ2 * UV(ĈT2)`.
@inline function _uv_to_ct12_components(uv, local_geometry)
    a = CC.Geometry.UVVector(
        CT12(CT1(unit_basis_vector_data(CT1, local_geometry)), local_geometry),
        local_geometry,
    )
    b = CC.Geometry.UVVector(
        CT12(CT2(unit_basis_vector_data(CT2, local_geometry)), local_geometry),
        local_geometry,
    )
    det = a.u * b.v - a.v * b.u
    τ1 = (uv.u * b.v - uv.v * b.u) / det
    τ2 = (a.u * uv.v - a.v * uv.u) / det
    return StaticArrays.SVector(τ1, τ2)
end

@inline _uv_to_ct1(uv, local_geometry) = _uv_to_ct12_components(uv, local_geometry)[1]
@inline _uv_to_ct2(uv, local_geometry) = _uv_to_ct12_components(uv, local_geometry)[2]

"""
    cartesian_to_contravariant!(ρτxz, ρτyz, uv_field)

Convert a `UVVector` field into the local CT1/CT2 unit-basis components
expected by the coupler flux fields (`csf.F_turb_ρτxz/yz`). Exact inverse of
`contravariant_to_cartesian!`.
"""
function cartesian_to_contravariant!(ρτxz, ρτyz, uv_field)
    local_geometry = CC.Fields.local_geometry_field(ρτxz)
    @. ρτxz = _uv_to_ct1(uv_field, local_geometry)
    @. ρτyz = _uv_to_ct2(uv_field, local_geometry)
    return nothing
end

"""
    scatter_poly_fluxes_to_boundary!(remapping, eg::ExchangeGrid,
                                     fs::ExchangeFluxState, weight)

Aggregate per-polygon fluxes onto the SE boundary space as a `weight`-weighted
average, filling the `remapping.flux_scratch` fields (`F_turb_ρτxz`,
`F_turb_ρτyz`, `F_sh`, `F_lh`, `F_turb_moisture`) ready for
`FluxCalculator.update_flux_fields!`.

`weight` is a per-polygon vector selecting the sub-surface the fluxes apply
to — `1 - sic` for the open ocean, `sic` for sea ice — so the nodal result is
a per-unit-*weighted*-area flux, consistent with the corresponding area
fraction the coupler multiplies it by. Steps:

1. raw L2 scatter of `weight * F` for each output (momentum in UV) and of
   `weight` itself (the nodal weighted coverage);
2. `weighted_dss!` on each scalar (UV components are DSS-safe; dividing only
   after DSS avoids averaging real coastal fluxes with no-data zeros);
3. division by the DSS'd weighted coverage (zero below
   `EXCHANGE_COV_CUTOFF`);
4. UV → CT1/CT2 conversion of the momentum pair.

Uses `fs.scratch1` internally; `weight` must not alias it.
"""
NVTX.@annotate function scatter_poly_fluxes_to_boundary!(
    remapping,
    eg::ExchangeGrid,
    fs::ExchangeFluxState,
    weight,
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
    cutoff = FT(EXCHANGE_COV_CUTOFF)
    for field in values(fx)
        Utilities.apply_dss!(field, remapping.flux_dss_buffer)
        @. field = ifelse(cov > cutoff, field / max(cov, cutoff), FT(0))
    end

    @. remapping.temp_uv_vec = CC.Geometry.UVVector(fx.F_turb_ρτxz, fx.F_turb_ρτyz)
    uv_to_ct12!(
        fx.F_turb_ρτxz,
        fx.F_turb_ρτyz,
        remapping.temp_uv_vec,
        remapping.momentum_basis,
    )
    return nothing
end
