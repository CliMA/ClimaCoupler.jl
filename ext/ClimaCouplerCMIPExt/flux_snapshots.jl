#=
# Flux snapshots

Host-side copies of every exchange-grid flux on each grid it passes through
(intersection polygons, SE boundary nodes, FV cells), for offline plotting with
`Plotting.plot_flux_snapshot`. Nothing here runs unless snapshots are requested
(`flux_snapshot_interval`) or a NaN appears (`flux_snapshot_on_nan`).
=#

"""
    exchange_grid_geometry(boundary_space, grid_oc, eg::ExchangeGrid)

Plain-array geometry for plotting exchange-grid data, as `(; poly, nodes, cells)`:

- `poly`: vertex rings in CSR form (`vert_ptr`, `vert_lon`, `vert_lat`, degrees,
  rings not closed), plus `area` [m²] and the owning `elem` and `cell`, in the
  same order as every per-polygon array.
- `nodes`: `lon`, `lat`, quadrature weight `Jw` [m²] and `elem` for each flat SE
  node, plus the node → polygon incidence of the scatter (`poly_ptr`, `poly`).
- `cells`: FV cell corners `corner_lon`/`corner_lat` (4 × `n_oc`, in index order
  `(i,j)`, `(i+1,j)`, `(i+1,j+1)`, `(i,j+1)`), centers `lon`/`lat`, `wet_area`
  [m²], and `Nx`, `Ny`. Cell `c = (j - 1) Nx + i`.

Recomputes the intersection polygons (seconds on a production grid), so call it once.
"""
function exchange_grid_geometry(boundary_space, grid_oc, eg::ExchangeGrid)
    CRExt = get_ConservativeRegriddingCCExt()
    GO, GI = CRExt.GO, CRExt.GI
    boundary_space_cpu = CC.Adapt.adapt(Array, boundary_space)
    eg_cpu = on_device(OC.CPU(), eg)

    # Polygons, in the kept numbering `build_exchange_grid` uses.
    (; polys, keep) = exchange_polygons(boundary_space_cpu, grid_oc)
    kept = polys[keep]
    @assert length(kept) == eg_cpu.n_poly
    to_lonlat = GO.GeographicFromUnitSphere()
    vert_ptr = Int32[1]
    vert_lon = Float64[]
    vert_lat = Float64[]
    for polygon in kept
        ring = GI.getexterior(polygon)
        for v in 1:(GI.npoint(ring) - 1) # rings repeat their first vertex
            lon, lat = to_lonlat(GI.getpoint(ring, v))
            push!(vert_lon, lon)
            push!(vert_lat, lat)
        end
        push!(vert_ptr, length(vert_lon) + 1)
    end
    poly = (;
        vert_ptr,
        vert_lon,
        vert_lat,
        area = Array(eg_cpu.area),
        elem = Array(eg_cpu.elem_of_poly),
        cell = Array(eg_cpu.oc_of_poly),
    )

    # SE nodes, flat in `se_nodal_vec` order.
    coords = CC.Fields.coordinate_field(boundary_space_cpu)
    Nq = CC.Quadratures.degrees_of_freedom(CC.Spaces.quadrature_style(boundary_space_cpu))
    lat = Float64.(CRExt.flat_nodal_data(CC.Fields.field_values(coords.lat)))
    lon = Float64.(CRExt.flat_nodal_data(CC.Fields.field_values(coords.long)))
    nodes = (;
        lon,
        lat,
        Jw = CRExt.se_node_weights(boundary_space_cpu),
        elem = Int32[div(n - 1, Nq^2) + 1 for n in eachindex(lat)],
        poly_ptr = Array(eg_cpu.snode_ptr),
        poly = Array(eg_cpu.spoly),
    )

    # FV cells.
    ug = OC.on_architecture(OC.CPU(), underlying_grid(grid_oc))
    Nx, Ny = size(ug, 1), size(ug, 2)
    corner_lon = Matrix{Float64}(undef, 4, Nx * Ny)
    corner_lat = similar(corner_lon)
    center_lon = Vector{Float64}(undef, Nx * Ny)
    center_lat = similar(center_lon)
    for j in 1:Ny, i in 1:Nx
        c = (j - 1) * Nx + i
        for (m, (di, dj)) in enumerate(((0, 0), (1, 0), (1, 1), (0, 1)))
            corner_lon[m, c] = ug.λᶠᶠᵃ[i + di, j + dj]
            corner_lat[m, c] = ug.φᶠᶠᵃ[i + di, j + dj]
        end
        center_lon[c] = ug.λᶜᶜᵃ[i, j]
        center_lat[c] = ug.φᶜᶜᵃ[i, j]
    end
    cells = (;
        corner_lon,
        corner_lat,
        lon = center_lon,
        lat = center_lat,
        wet_area = Array(eg_cpu.oc_wet_area),
        Nx,
        Ny,
    )
    return (; poly, nodes, cells)
end

# The exchange-grid fluxes a snapshot records, named as in `ExchangeFluxState`.
const SNAPSHOT_TURBULENT_FLUXES = (:F_sh, :F_lh, :F_moisture, :F_τu, :F_τv)
const SNAPSHOT_DOWNWARD_FLUXES = (:SW_d, :LW_d, :P_liq, :P_snow)

_to_host(x) = Array(x)
_nodes_to_host(field::CC.Fields.Field) = Array(se_nodal_vec(field))
_nodes_to_host(component) = Array(vec(parent(component))) # a `UVVector` component
_names_to_host(source, names) =
    NamedTuple{names}(map(name -> _to_host(getproperty(source, name)), names))

"""
    SimOutput.capture_flux_snapshot(sim::OceananigansSimulation, csf)

Snapshot of the atmosphere–ocean fluxes on the exchange grid, the SE boundary
nodes and the ocean cells; see the flux snapshots spec for the layout. Overwrites
scratch only (`flux_scratch`, `weight_cov_scratch`, `temp_uv_vec`,
`ocean_flux_state.scratch1/2`), all of which the next coupling step refills
before use. Returns `nothing` without the exchange grid.
"""
function SimOutput.capture_flux_snapshot(sim::OceananigansSimulation, csf)
    remapping = sim.remapping
    remapping.use_exchange_grid || return nothing
    eg = remapping.exchange_grid
    fs = remapping.ocean_flux_state
    model = sim.ocean.model
    grid = model.grid

    # Exchange grid. `fs` holds this step's values: downward fluxes from
    # `update_sim!`, turbulent ones from `compute_surface_fluxes!`.
    @. fs.scratch2 = 1 - fs.sic # open-water weight of each polygon
    poly = merge(
        _names_to_host(fs, (SNAPSHOT_TURBULENT_FLUXES..., SNAPSHOT_DOWNWARD_FLUXES...)),
        (;
            T_sfc = _to_host(fs.T_sfc),
            sic = _to_host(fs.sic),
            weight = _to_host(fs.scratch2),
        ),
    )

    # SE nodes, this surface only. Ocean and ice share `flux_scratch`, so redo the
    # scatter instead of trusting whichever surface wrote it last. Momentum comes
    # from `temp_uv_vec`, where the scatter leaves it in the east/north basis.
    scatter_poly_fluxes_to_boundary!(remapping, eg, fs, fs.scratch2)
    fx = remapping.flux_scratch
    uv = remapping.temp_uv_vec.components.data
    nodes = (;
        F_sh = _nodes_to_host(fx.F_sh),
        F_lh = _nodes_to_host(fx.F_lh),
        F_moisture = _nodes_to_host(fx.F_turb_moisture),
        F_τu = _nodes_to_host(uv.:1),
        F_τv = _nodes_to_host(uv.:2),
        cov = _nodes_to_host(remapping.weight_cov_scratch),
    )

    # SE nodes, what the atmosphere receives: every surface, area-weighted.
    contravariant_to_cartesian!(
        remapping.temp_uv_vec,
        csf.F_turb_ρτxz,
        csf.F_turb_ρτyz,
        remapping.uv_basis,
    )
    atmos = (;
        F_sh = _nodes_to_host(csf.F_sh),
        F_lh = _nodes_to_host(csf.F_lh),
        F_moisture = _nodes_to_host(csf.F_turb_moisture),
        F_τu = _nodes_to_host(uv.:1),
        F_τv = _nodes_to_host(uv.:2),
        SW_d = _nodes_to_host(csf.SW_d),
        LW_d = _nodes_to_host(csf.LW_d),
        P_liq = _nodes_to_host(csf.P_liq),
        P_snow = _nodes_to_host(csf.P_snow),
        ocean_fraction = _nodes_to_host(csf.ocean_area_fraction),
        ice_fraction = _nodes_to_host(csf.ice_area_fraction),
        land_fraction = _nodes_to_host(csf.land_area_fraction),
    )

    # Ocean cells: each flux weighted by open water and scattered like the BCs,
    # but in the flux's own units (before density/heat-capacity scaling).
    cell_buffer = similar(eg.oc_wet_area)
    cells = NamedTuple{SNAPSHOT_TURBULENT_FLUXES}(
        map(SNAPSHOT_TURBULENT_FLUXES) do name
            flux = getproperty(fs, name)
            @. fs.scratch1 = fs.scratch2 * flux
            scatter_polys_to_cells!(cell_buffer, eg, fs.scratch1)
            mirror_fold_partners!(cell_buffer, grid)
            _to_host(cell_buffer)
        end,
    )

    # What the ocean model receives this step, in its BC units (kinematic), and
    # the surface state it was computed from. `u` and `v` are staggered, so trim
    # every field to the tracer block.
    Nx, Ny, Nz = size(grid)
    surface(field, k) = _to_host(vec(OC.interior(field, 1:Nx, 1:Ny, k)))
    totals = (;
        T = surface(surface_flux(model.tracers.T), 1),
        S = surface(surface_flux(model.tracers.S), 1),
        u = surface(surface_flux(model.velocities.u), 1),
        v = surface(surface_flux(model.velocities.v), 1),
    )
    context =
        (; SST = surface(model.tracers.T, Nz), sic = surface(sim.ice_concentration, 1))

    return (; atmos, ocean = (; poly, nodes, cells, totals, context))
end

function SimOutput.flux_snapshot_geometry(sim::OceananigansSimulation, csf)
    sim.remapping.use_exchange_grid || return nothing
    return exchange_grid_geometry(
        axes(csf.F_sh),
        sim.ocean.model.grid,
        sim.remapping.exchange_grid,
    )
end
