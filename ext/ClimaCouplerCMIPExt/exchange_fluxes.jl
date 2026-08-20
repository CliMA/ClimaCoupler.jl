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
import ClimaComms

"""
    ExchangeFluxState{FT, VF}

Per-polygon scratch for ocean or ice exchange-grid fluxes: atmos + surface
state, fluxes `F_*`, accumulators `acc_*` / `n_acc`, and two scratch vectors.
Paired ocean/ice states share the nine atmos vectors
([`paired_exchange_flux_states`](@ref)).
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

# Leading atmos fields through `height_sfc` (shared by paired ocean/ice states).
@inline _atmos_fields(fs::ExchangeFluxState) = (
    fs.T_atmos,
    fs.q_tot,
    fs.q_liq,
    fs.q_ice,
    fs.ρ_atmos,
    fs.u_atmos,
    fs.v_atmos,
    fs.height_int,
    fs.height_sfc,
)

function _cpu_exchange_flux_state(FT, n_poly::Int; atmos = nothing)
    n_atmos = 9
    n_rest = fieldcount(ExchangeFluxState) - 1 - n_atmos
    atmos_vecs = isnothing(atmos) ? ntuple(_ -> zeros(FT, n_poly), n_atmos) : atmos
    @assert length(atmos_vecs) == n_atmos
    return ExchangeFluxState{FT, Vector{FT}}(
        atmos_vecs...,
        ntuple(_ -> zeros(FT, n_poly), n_rest)...,
        Ref(0),
    )
end

ExchangeFluxState{FT}(arch, n_poly::Int) where {FT} =
    on_device(arch, _cpu_exchange_flux_state(FT, n_poly))

"""
    IceExchangeState{FT, VF}

[`ExchangeFluxState`](@ref) plus ice skin-balance inputs (`R`, `T_i`, `SW_d`,
`LW_d`) and diagnosed `T_sfc_new`.
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
    n_ice = fieldcount(IceExchangeState) - 1
    state = IceExchangeState(
        _cpu_exchange_flux_state(FT, n_poly),
        ntuple(_ -> zeros(FT, n_poly), n_ice)...,
    )
    return on_device(arch, state)
end

"""
    paired_exchange_flux_states(FT, arch, n_poly)

Ocean + ice flux state sharing atmos polygon vectors. Ice gathers first in
`turbulent_fluxes!`; ocean reuses via `ensure_atmos_gathered!`.
"""
function paired_exchange_flux_states(::Type{FT}, arch, n_poly::Int) where {FT}
    ocean = ExchangeFluxState{FT}(arch, n_poly)
    AT = typeof(ocean.T_atmos)
    z() = similar(ocean.T_atmos)
    n_atmos = 9
    n_rest = fieldcount(ExchangeFluxState) - 1 - n_atmos
    ice_fluxes = ExchangeFluxState{FT, AT}(
        _atmos_fields(ocean)...,
        ntuple(_ -> z(), n_rest)...,
        Ref(0),
    )
    ice = IceExchangeState{FT, AT}(
        ice_fluxes,
        ntuple(_ -> z(), fieldcount(IceExchangeState) - 1)...,
    )
    return ocean, ice
end

# Gather atmos onto shared polygon vectors once per coupling step.
NVTX.@annotate function ensure_atmos_gathered!(remapping, csf)
    remapping.atmos_gathered[] && return nothing
    gather_atmos_state_to_polys!(
        remapping.ocean_flux_state,
        remapping.exchange_grid,
        csf,
        remapping.temp_uv_vec,
        remapping.uv_basis,
    )
    remapping.atmos_gathered[] = true
    return nothing
end

"""
    average_and_reset_exchange_accumulators!(fs) -> Bool

Average `acc_*` into `F_*`, then zero accumulators. Returns `false` if `n_acc == 0`.
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

# Gather atmos SE nodes → polygons. Convert winds CT1/CT2 → UV first so
# stresses stay basis-consistent across nodes that touch a polygon.
NVTX.@annotate function gather_atmos_state_to_polys!(
    fs::ExchangeFluxState,
    eg::ExchangeGrid,
    csf,
    temp_uv_vec,
    uv_basis,
)
    contravariant_to_cartesian!(temp_uv_vec, csf.u_int, csf.v_int, uv_basis)
    u_uv = temp_uv_vec.components.data.:1
    v_uv = temp_uv_vec.components.data.:2
    gather_nodes_to_polys!(
        (
            fs.u_atmos,
            fs.v_atmos,
            fs.T_atmos,
            fs.q_tot,
            fs.q_liq,
            fs.q_ice,
            fs.ρ_atmos,
            fs.height_int,
            fs.height_sfc,
        ),
        eg,
        (
            se_nodal_vec(u_uv),
            se_nodal_vec(v_uv),
            se_nodal_vec(csf.T_atmos),
            se_nodal_vec(csf.q_tot_atmos),
            se_nodal_vec(csf.q_liq_atmos),
            se_nodal_vec(csf.q_ice_atmos),
            se_nodal_vec(csf.ρ_atmos),
            se_nodal_vec(csf.height_int),
            se_nodal_vec(csf.height_sfc),
        ),
    )
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
    backend = KernelAbstractions.get_backend(fs.F_sh)
    launch_kernel!(
        _ocean_polygon_fluxes_kernel!,
        backend,
        length(fs.F_sh),
        _kernel_state(fs),
        surface_fluxes_params,
        thermo_params,
        config,
    )
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
    solver_opts,
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
    update_T_sfc_cb = UpdateTSfc(R, T_i, σ, ϵ, SW_d, LW_d, α_albedo, T_melt)
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
        update_T_sfc_cb;
        solver_opts,
    )
end

@kernel function _ice_polygon_fluxes_kernel!(
    vecs,
    surface_fluxes_params,
    thermo_params,
    config,
    solver_opts,
    σ,
    ϵ,
    α_albedo,
    T_melt,
)
    k = @index(Global)
    @inbounds out = _polygon_ice_fluxes(
        surface_fluxes_params,
        thermo_params,
        config,
        solver_opts,
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
    _store_ocean_polygon_fluxes!(vecs, k, out)
    @inbounds vecs.T_sfc_new[k] = out.T_sfc_new
end

"""
    compute_ice_polygon_fluxes!(is, …)

Per-polygon SurfaceFluxes + skin temperature (`sic > 0`); store fluxes / `T_sfc_new`
and accumulate.
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
    FT = eltype(is.fluxes.F_sh)
    solver_opts = ice_surface_flux_solver_opts(FT)
    backend = KernelAbstractions.get_backend(is.fluxes.F_sh)
    launch_kernel!(
        _ice_polygon_fluxes_kernel!,
        backend,
        length(is.fluxes.F_sh),
        _kernel_state(is),
        surface_fluxes_params,
        thermo_params,
        config,
        solver_opts,
        σ,
        ϵ,
        α_albedo,
        T_melt,
    )
    return nothing
end

"""
    cartesian_to_contravariant!(ρτxz, ρτyz, uv_field, uv_basis)
    cartesian_to_contravariant!(ρτxz, ρτyz, uv_field)

Convert a `UVVector` field into the local CT1/CT2 unit-basis components
expected by the coupler flux fields (`csf.F_turb_ρτxz/yz`). Exact inverse of
`contravariant_to_cartesian!`. The 4-argument form applies precomputed
[`uv_basis_coefficients`](@ref); the 3-argument form builds them on the fly.
"""
cartesian_to_contravariant!(ρτxz, ρτyz, uv_field) =
    cartesian_to_contravariant!(ρτxz, ρτyz, uv_field, uv_basis_coefficients(axes(ρτxz)))

function cartesian_to_contravariant!(ρτxz, ρτyz, uv_field, uv_basis)
    a = parent(ρτxz)
    b = parent(ρτyz)
    u = parent(uv_field.components.data.:1)
    v = parent(uv_field.components.data.:2)
    p11 = parent(uv_basis.u_to_ct1)
    p12 = parent(uv_basis.v_to_ct1)
    p21 = parent(uv_basis.u_to_ct2)
    p22 = parent(uv_basis.v_to_ct2)
    @. a = p11 * u + p12 * v
    @. b = p21 * u + p22 * v
    return nothing
end

"""
    scalar_weighted_dss!(field::CC.Fields.Field, buffer)

Weighted DSS for a scalar field, bit-identical to `Spaces.weighted_dss!` but
without loading `LocalGeometry` (which `dss_transform` does even for scalars,
boxing a 22-float struct per perimeter node at `-O0`). Multiplies by
`dss_weights` then reuses `buffer` for the unweighted topology DSS.
"""
function scalar_weighted_dss!(field::CC.Fields.Field, buffer)
    isnothing(buffer) && return nothing
    data = CC.Fields.field_values(field)
    space = axes(field)
    topo = CC.Spaces.topology(space)
    device = ClimaComms.device(topo)
    perimeter = CC.Topologies.Perimeter2D(size(data, 2))
    dp = parent(data)
    wp = parent(CC.Spaces.dss_weights(space))
    @. dp = dp * wp
    CC.Topologies.dss_load_perimeter_data!(device, buffer, data, perimeter)
    CC.Topologies.dss_local_ghost!(device, buffer.perimeter_data, perimeter, topo)
    CC.Topologies.fill_send_buffer!(device, buffer)
    ClimaComms.start(buffer.graph_context)
    CC.Topologies.dss_local!(device, buffer.perimeter_data, perimeter, topo)
    ClimaComms.finish(buffer.graph_context)
    CC.Topologies.load_from_recv_buffer!(device, buffer)
    CC.Topologies.dss_ghost!(device, buffer.perimeter_data, perimeter, topo)
    CC.Topologies.dss_unload_perimeter_data!(device, data, buffer, perimeter)
    return nothing
end

function _dss_and_normalize!(field, cov, buffer, cutoff, z)
    scalar_weighted_dss!(field, buffer)
    @. field = ifelse(cov > cutoff, field / max(cov, cutoff), z)
    return nothing
end

"""
    scatter_poly_fluxes_to_boundary!(remapping, eg, fs, weight; cov_cutoff=1e-3)

Weight-average polygon fluxes onto SE boundary scratch (`weight` = `1-sic` or
`sic`): fused L2 scatter of `weight` and `weight*F`, DSS, normalize by coverage,
then UV → CT1/CT2 for momentum. Low-coverage nodes (`< cov_cutoff`) get zero.
"""
NVTX.@annotate function scatter_poly_fluxes_to_boundary!(
    remapping,
    eg::ExchangeGrid,
    fs::ExchangeFluxState,
    weight;
    cov_cutoff = 1e-3,
)
    fx = remapping.flux_scratch
    cov = remapping.weight_cov_scratch
    scatter_polys_to_nodes_weighted!(
        (
            se_nodal_vec(fx.F_sh),
            se_nodal_vec(fx.F_lh),
            se_nodal_vec(fx.F_turb_moisture),
            se_nodal_vec(fx.F_turb_ρτxz),
            se_nodal_vec(fx.F_turb_ρτyz),
        ),
        se_nodal_vec(cov),
        eg,
        (fs.F_sh, fs.F_lh, fs.F_moisture, fs.F_τu, fs.F_τv),
        weight,
    )

    FT = CC.Spaces.undertype(axes(cov))
    cutoff = FT(cov_cutoff)
    z = zero(FT)
    buf = remapping.flux_dss_buffer
    scalar_weighted_dss!(cov, buf)
    for f in (fx.F_sh, fx.F_lh, fx.F_turb_moisture, fx.F_turb_ρτxz, fx.F_turb_ρτyz)
        _dss_and_normalize!(f, cov, buf, cutoff, z)
    end

    parent(remapping.temp_uv_vec.components.data.:1) .= parent(fx.F_turb_ρτxz)
    parent(remapping.temp_uv_vec.components.data.:2) .= parent(fx.F_turb_ρτyz)
    cartesian_to_contravariant!(
        fx.F_turb_ρτxz,
        fx.F_turb_ρτyz,
        remapping.temp_uv_vec,
        remapping.uv_basis,
    )
    return nothing
end
