import NCDatasets
import ClimaUtilities.ClimaArtifacts: @clima_artifact

using Oceananigans.BoundaryConditions:
    fill_halo_regions!, FPivotZipperBoundaryCondition, NoFluxBoundaryCondition, FieldBoundaryConditions
using Oceananigans.Grids: RightFaceFolded, generate_coordinate
using Oceananigans.OrthogonalSphericalShellGrids: Tripolar

const orca_topology = (OC.Periodic, RightFaceFolded, OC.Bounded)

const orca_metrics = (
    (:λᶜᶜᵃ, "lambda_cca", OC.Center, OC.Center), (:λᶠᶜᵃ, "lambda_fca", OC.Face, OC.Center),
    (:λᶜᶠᵃ, "lambda_cfa", OC.Center, OC.Face),   (:λᶠᶠᵃ, "lambda_ffa", OC.Face, OC.Face),
    (:φᶜᶜᵃ, "phi_cca", OC.Center, OC.Center),    (:φᶠᶜᵃ, "phi_fca", OC.Face, OC.Center),
    (:φᶜᶠᵃ, "phi_cfa", OC.Center, OC.Face),      (:φᶠᶠᵃ, "phi_ffa", OC.Face, OC.Face),
    (:Δxᶜᶜᵃ, "dx_cca", OC.Center, OC.Center),    (:Δxᶠᶜᵃ, "dx_fca", OC.Face, OC.Center),
    (:Δxᶜᶠᵃ, "dx_cfa", OC.Center, OC.Face),      (:Δxᶠᶠᵃ, "dx_ffa", OC.Face, OC.Face),
    (:Δyᶜᶜᵃ, "dy_cca", OC.Center, OC.Center),    (:Δyᶠᶜᵃ, "dy_fca", OC.Face, OC.Center),
    (:Δyᶜᶠᵃ, "dy_cfa", OC.Center, OC.Face),      (:Δyᶠᶠᵃ, "dy_ffa", OC.Face, OC.Face),
    (:Azᶜᶜᵃ, "area_cca", OC.Center, OC.Center),  (:Azᶠᶜᵃ, "area_fca", OC.Face, OC.Center),
    (:Azᶜᶠᵃ, "area_cfa", OC.Center, OC.Face),    (:Azᶠᶠᵃ, "area_ffa", OC.Face, OC.Face),
)

"""
    orca_one_grid_path()

Path to the eORCA1 mesh file inside the `orca_one_grid` artifact.
"""
orca_one_grid_path() = joinpath(@clima_artifact("orca_one_grid"), "orca_one_grid.nc")

"""
    read_orca_mesh(path = orca_one_grid_path())

Read the eORCA1 metric arrays, bottom height and conformal-mapping parameters from `path`.
"""
function read_orca_mesh(path = orca_one_grid_path())
    return NCDatasets.NCDataset(path) do ds
        metrics = NamedTuple(
            symbol => Array(ds[variable_name][:, :]) for (symbol, variable_name, _, _) in orca_metrics
        )
        (;
            metrics,
            bottom_height = Array(ds["bottom_height"][:, :]),
            Nx = Int(ds.attrib["Nx"]),
            Ny = Int(ds.attrib["Ny"]),
            radius = ds.attrib["radius"],
            north_poles_latitude = ds.attrib["north_poles_latitude"],
            first_pole_longitude = ds.attrib["first_pole_longitude"],
            southernmost_latitude = ds.attrib["southernmost_latitude"],
        )
    end
end

"""
    halo_filled_metric(data, helper_grid, boundary_conditions, LX, LY)

`data` placed in a field of the given stagger location and extended into the halos.

`fill_halo_regions!` resolves the periodic east/west fill and the north fold from the helper grid's
topology, so the halo values match what the mesh generator would have produced.
"""
function halo_filled_metric(data, helper_grid, boundary_conditions, LX, LY)
    field = OC.Field{LX, LY, OC.Center}(helper_grid; boundary_conditions)
    Nx, Ny = size(data)
    field.data[1:Nx, 1:Ny, 1] .= data
    fill_halo_regions!(field)
    return deepcopy(dropdims(field.data, dims = 3))
end

"""
    ORCAOneGrid(arch = OC.CPU(), FT = Float64; Nz, z, halo, active_cells_map)

Build the eORCA1 `ImmersedBoundaryGrid` from the `orca_one_grid` artifact.

`Nz`, `z` and `halo` are free: only the horizontal mesh comes from the artifact, and the bottom
height it stores is unsnapped, so `GridFittedBottom` discretizes it against the vertical grid
requested here.
"""
function ORCAOneGrid(
    arch = OC.CPU(),
    FT::DataType = Float64;
    Nz = 60,
    z = (-6000, 0),
    halo = (5, 5, 4),
    active_cells_map = true,
    mesh_path = orca_one_grid_path(),
)
    mesh = read_orca_mesh(mesh_path)
    Nx, Ny = mesh.Nx, mesh.Ny
    Hx, Hy, Hz = halo

    Lz, z_coordinate = generate_coordinate(FT, orca_topology, (Nx, Ny, Nz), halo, z, :z, 3, OC.CPU())

    helper_grid = OC.RectilinearGrid(;
        size = (Nx, Ny),
        halo = (Hx, Hy),
        x = (0, 1),
        y = (0, 1),
        topology = (OC.Periodic, RightFaceFolded, OC.Flat),
    )

    boundary_conditions = FieldBoundaryConditions(
        north = FPivotZipperBoundaryCondition(),
        south = NoFluxBoundaryCondition(),
        west = OC.PeriodicBoundaryCondition(),
        east = OC.PeriodicBoundaryCondition(),
        top = nothing,
        bottom = nothing,
    )

    to_architecture(data) = OC.Architectures.on_architecture(arch, map(FT, data))

    filled = NamedTuple(
        symbol => to_architecture(
            halo_filled_metric(mesh.metrics[symbol], helper_grid, boundary_conditions, LX, LY),
        ) for (symbol, _, LX, LY) in orca_metrics
    )

    underlying_grid = OC.OrthogonalSphericalShellGrid{orca_topology...}(
        arch,
        Nx,
        Ny,
        Nz,
        Hx,
        Hy,
        Hz,
        convert(FT, Lz),
        filled.λᶜᶜᵃ,
        filled.λᶠᶜᵃ,
        filled.λᶜᶠᵃ,
        filled.λᶠᶠᵃ,
        filled.φᶜᶜᵃ,
        filled.φᶠᶜᵃ,
        filled.φᶜᶠᵃ,
        filled.φᶠᶠᵃ,
        OC.Architectures.on_architecture(arch, z_coordinate),
        filled.Δxᶜᶜᵃ,
        filled.Δxᶠᶜᵃ,
        filled.Δxᶜᶠᵃ,
        filled.Δxᶠᶠᵃ,
        filled.Δyᶜᶜᵃ,
        filled.Δyᶠᶜᵃ,
        filled.Δyᶜᶠᵃ,
        filled.Δyᶠᶠᵃ,
        filled.Azᶜᶜᵃ,
        filled.Azᶠᶜᵃ,
        filled.Azᶜᶠᵃ,
        filled.Azᶠᶠᵃ,
        convert(FT, mesh.radius),
        Tripolar(mesh.north_poles_latitude, mesh.first_pole_longitude, mesh.southernmost_latitude),
    )

    bottom_field = OC.Field{OC.Center, OC.Center, Nothing}(underlying_grid)
    OC.set!(bottom_field, OC.Architectures.on_architecture(arch, mesh.bottom_height))

    return OC.ImmersedBoundaryGrid(
        underlying_grid,
        OC.GridFittedBottom(bottom_field);
        active_cells_map,
    )
end
