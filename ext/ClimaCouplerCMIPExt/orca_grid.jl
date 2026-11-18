import NCDatasets
import ClimaUtilities.ClimaArtifacts: @clima_artifact

using Oceananigans.BoundaryConditions:
    fill_halo_regions!,
    FPivotZipperBoundaryCondition,
    TPivotZipperBoundaryCondition,
    NoFluxBoundaryCondition,
    FieldBoundaryConditions,
    FPivot,
    TPivot
using Oceananigans.Grids:
    CenterFolded, FaceFolded, RightCenterFolded, RightFaceFolded, generate_coordinate
using Oceananigans.OrthogonalSphericalShellGrids: Tripolar, distribute_tripolar_grid

const orca_metrics = (
    (:λᶜᶜᵃ, "lambda_cca", OC.Center, OC.Center),
    (:λᶠᶜᵃ, "lambda_fca", OC.Face, OC.Center),
    (:λᶜᶠᵃ, "lambda_cfa", OC.Center, OC.Face),
    (:λᶠᶠᵃ, "lambda_ffa", OC.Face, OC.Face),
    (:φᶜᶜᵃ, "phi_cca", OC.Center, OC.Center),
    (:φᶠᶜᵃ, "phi_fca", OC.Face, OC.Center),
    (:φᶜᶠᵃ, "phi_cfa", OC.Center, OC.Face),
    (:φᶠᶠᵃ, "phi_ffa", OC.Face, OC.Face),
    (:Δxᶜᶜᵃ, "dx_cca", OC.Center, OC.Center),
    (:Δxᶠᶜᵃ, "dx_fca", OC.Face, OC.Center),
    (:Δxᶜᶠᵃ, "dx_cfa", OC.Center, OC.Face),
    (:Δxᶠᶠᵃ, "dx_ffa", OC.Face, OC.Face),
    (:Δyᶜᶜᵃ, "dy_cca", OC.Center, OC.Center),
    (:Δyᶠᶜᵃ, "dy_fca", OC.Face, OC.Center),
    (:Δyᶜᶠᵃ, "dy_cfa", OC.Center, OC.Face),
    (:Δyᶠᶠᵃ, "dy_ffa", OC.Face, OC.Face),
    (:Azᶜᶜᵃ, "area_cca", OC.Center, OC.Center),
    (:Azᶠᶜᵃ, "area_fca", OC.Face, OC.Center),
    (:Azᶜᶠᵃ, "area_cfa", OC.Center, OC.Face),
    (:Azᶠᶠᵃ, "area_ffa", OC.Face, OC.Face),
)

"""
    orca_grid_path(artifact = "orca_one_grid", context = ClimaComms.SingletonCommsContext())

Path to the eORCA mesh file inside `artifact`: `orca_one_grid` (eORCA1), `orca_half_grid` (eORCA05)
or `orca_quarter_grid` (eORCA025). The latter two are lazy: the root of `context` downloads them
while the other processes wait, so every process of `context` must call this.
"""
orca_grid_path(artifact = "orca_one_grid", context = ClimaComms.SingletonCommsContext()) =
    joinpath(@clima_artifact(artifact, context), artifact * ".nc")

"""
    orca_fold(north_fold_pivot)

The north-fold topology, `Tripolar` fold type and zipper boundary condition for an eORCA mesh whose
fold pivots about `north_fold_pivot`: `"FPivot"` (eORCA1, eORCA05) or `"TPivot"` (eORCA025).
"""
function orca_fold(north_fold_pivot)
    if north_fold_pivot == "FPivot"
        return (; topology = RightFaceFolded, folded = FaceFolded, pivot = FPivot,
                  zipper = FPivotZipperBoundaryCondition())
    elseif north_fold_pivot == "TPivot"
        return (; topology = RightCenterFolded, folded = CenterFolded, pivot = TPivot,
                  zipper = TPivotZipperBoundaryCondition())
    else
        error("Unsupported eORCA north-fold pivot $north_fold_pivot")
    end
end

"""
    read_orca_mesh(path = orca_grid_path())

Read the eORCA metric arrays, bottom height, conformal-mapping parameters and north-fold pivot from
`path`. The published `orca_one_grid` predates the `north_fold_pivot` attribute, so its absence means
the eORCA1 F pivot.
"""
function read_orca_mesh(path = orca_grid_path())
    return NCDatasets.NCDataset(path) do ds
        metrics = NamedTuple(
            symbol => Array(ds[variable_name][:, :]) for
            (symbol, variable_name, _, _) in orca_metrics
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
            north_fold_pivot = get(ds.attrib, "north_fold_pivot", "FPivot"),
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
    ORCAGrid(arch = OC.CPU(), FT = Float64; Nz, z, halo, active_cells_map, mesh_path)

Build an eORCA `ImmersedBoundaryGrid` from the mesh at `mesh_path` (default: the eORCA1
`orca_one_grid` artifact; see [`orca_grid_path`](@ref)).

`Nz`, `z` and `halo` are free: only the horizontal mesh comes from the artifact, and the bottom
height it stores is unsnapped, so `GridFittedBottom` discretizes it against the vertical grid
requested here.

On a `Distributed` `arch` the global mesh is built on the CPU and split with Oceananigans'
`distribute_tripolar_grid`, as `TripolarGrid` does; uneven splits give the remainder to the last rank.
"""
function ORCAGrid(
    arch = OC.CPU(),
    FT::DataType = Float64;
    Nz = 60,
    z = (-6000, 0),
    halo = (5, 5, 5),
    active_cells_map = true,
    mesh_path = orca_grid_path(),
)
    mesh = read_orca_mesh(mesh_path)
    Nx, Ny = mesh.Nx, mesh.Ny
    Hx, Hy, Hz = halo
    fold = orca_fold(mesh.north_fold_pivot)
    orca_topology = (OC.Periodic, fold.topology, OC.Bounded)

    Lz, z_coordinate =
        generate_coordinate(FT, orca_topology, (Nx, Ny, Nz), halo, z, :z, 3, OC.CPU())

    helper_grid = OC.RectilinearGrid(;
        size = (Nx, Ny),
        halo = (Hx, Hy),
        x = (0, 1),
        y = (0, 1),
        topology = (OC.Periodic, fold.topology, OC.Flat),
    )

    boundary_conditions = FieldBoundaryConditions(
        north = fold.zipper,
        south = NoFluxBoundaryCondition(),
        west = OC.PeriodicBoundaryCondition(),
        east = OC.PeriodicBoundaryCondition(),
        top = nothing,
        bottom = nothing,
    )

    global_arch = arch isa OC.Distributed ? OC.CPU() : arch
    to_architecture(data) = OC.Architectures.on_architecture(global_arch, map(FT, data))

    filled = NamedTuple(
        symbol => to_architecture(
            halo_filled_metric(
                mesh.metrics[symbol],
                helper_grid,
                boundary_conditions,
                LX,
                LY,
            ),
        ) for (symbol, _, LX, LY) in orca_metrics
    )

    global_grid = OC.OrthogonalSphericalShellGrid{orca_topology...}(
        global_arch,
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
        OC.Architectures.on_architecture(global_arch, z_coordinate),
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
        Tripolar(
            mesh.north_poles_latitude,
            mesh.first_pole_longitude,
            mesh.southernmost_latitude,
            fold.folded,
            fold.pivot,
        ),
    )

    underlying_grid =
        arch isa OC.Distributed ? distribute_tripolar_grid(arch, global_grid) : global_grid

    # `set!` takes the global array and keeps each rank's part on a distributed grid
    bottom_field = OC.Field{OC.Center, OC.Center, Nothing}(underlying_grid)
    OC.set!(bottom_field, mesh.bottom_height)

    return OC.ImmersedBoundaryGrid(
        underlying_grid,
        OC.GridFittedBottom(bottom_field);
        active_cells_map,
    )
end
