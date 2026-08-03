import NCDatasets
import Dates
import ClimaUtilities.ClimaArtifacts: @clima_artifact

# Ocean and sea-ice input data is read from pre-generated artifacts. The initial conditions are
# stored on their source datasets' own longitude-latitude(-depth) grids, together with the grid
# interfaces needed to rebuild those grids, so a single artifact serves every ocean grid. The
# generation scripts live in `ClimaArtifacts`.

const supported_initial_condition_dates = (Dates.Date(2010, 1, 1),)

"""
    check_supported_date(start_date)

Verify that ocean and sea-ice initial conditions exist for `start_date`.
"""
function check_supported_date(start_date)
    date = start_date isa Dates.DateTime ? Dates.Date(start_date) : start_date
    date in supported_initial_condition_dates && return nothing

    error("""
          Ocean and sea-ice initial conditions are available only for \
          $(join(supported_initial_condition_dates, ", ")); got $(date).
          """)
end

"""
    GriddedDataset

A field on a regular longitude-latitude(-depth) grid, read from an artifact, together with the grid
interfaces that reconstruct it. `z_interfaces` is empty for two-dimensional fields.
"""
struct GriddedDataset{D, L, P, Z}
    data::D
    longitude_interfaces::L
    latitude_interfaces::P
    z_interfaces::Z
end

"""
    read_gridded_dataset(path, variable_name; with_depth = true)

Read `variable_name` and its grid interfaces from the NetCDF file at `path`.
"""
function read_gridded_dataset(path, variable_name; with_depth = true)
    return NCDatasets.NCDataset(path) do ds
        data =
            with_depth ? Array(ds[variable_name][:, :, :]) : Array(ds[variable_name][:, :])
        GriddedDataset(
            data,
            Array(ds["longitude_interfaces"][:]),
            Array(ds["latitude_interfaces"][:]),
            with_depth ? Array(ds["z_interfaces"][:]) : Float64[],
        )
    end
end

"""
    source_field(dataset, arch)

Rebuild the dataset's native `LatitudeLongitudeGrid` and populate a field on it.

Halos are filled so that `interpolate!` can reach across the periodic seam in longitude.
"""
function source_field(dataset, arch)
    Nx, Ny = size(dataset.data, 1), size(dataset.data, 2)
    longitude = dataset.longitude_interfaces
    latitude = dataset.latitude_interfaces

    if isempty(dataset.z_interfaces)
        grid = OC.LatitudeLongitudeGrid(
            arch,
            Float64;
            size = (Nx, Ny),
            halo = (3, 3),
            longitude,
            latitude,
            topology = (OC.Periodic, OC.Bounded, OC.Flat),
        )
        field = OC.Field{OC.Center, OC.Center, Nothing}(grid)
    else
        Nz = size(dataset.data, 3)
        grid = OC.LatitudeLongitudeGrid(
            arch,
            Float64;
            size = (Nx, Ny, Nz),
            halo = (3, 3, 3),
            longitude,
            latitude,
            z = dataset.z_interfaces,
        )
        field = OC.CenterField(grid)
    end

    OC.set!(field, dataset.data)
    OC.BoundaryConditions.fill_halo_regions!(field)

    return field
end

"""
    interpolate_onto!(target_field, dataset)

Interpolate `dataset` from its native grid onto `target_field`.

Since the source data has no gaps, this is plain interpolation with no masking.
"""
function interpolate_onto!(target_field, dataset)
    arch = OC.Architectures.architecture(target_field.grid)
    OC.Fields.interpolate!(target_field, source_field(dataset, arch))
    return target_field
end

"""
    set_ocean_initial_conditions!(ocean_model, start_date)

Initialize ocean temperature and salinity from the EN4 state estimate for `start_date`.
"""
function set_ocean_initial_conditions!(ocean_model, start_date)
    check_supported_date(start_date)

    year = Dates.year(start_date)
    month = lpad(Dates.month(start_date), 2, '0') # e.g., 1 → "01"
    directory = @clima_artifact("en4_ocean_initial_conditions_$(year)_$month")
    path = joinpath(directory, "en4_ocean_initial_conditions_$(year)_$month.nc")
    @info "Initializing ocean temperature and salinity from $path"

    interpolate_onto!(ocean_model.tracers.T, read_gridded_dataset(path, "temperature"))
    interpolate_onto!(ocean_model.tracers.S, read_gridded_dataset(path, "salinity"))

    return nothing
end

"""
    set_sea_ice_initial_conditions!(sea_ice_model, start_date)

Initialize sea-ice concentration and thickness from the ECCO4 state estimate for `start_date`.
"""
function set_sea_ice_initial_conditions!(sea_ice_model, start_date)
    check_supported_date(start_date)

    year = Dates.year(start_date)
    month = lpad(Dates.month(start_date), 2, '0') # e.g., 1 → "01"
    directory = @clima_artifact("ecco4_seaice_initial_conditions_$(year)_$month")
    path = joinpath(directory, "ecco4_seaice_initial_conditions_$(year)_$month.nc")
    @info "Initializing sea-ice concentration and thickness from $path"

    concentration = read_gridded_dataset(path, "sea_ice_concentration"; with_depth = false)
    thickness = read_gridded_dataset(path, "sea_ice_thickness"; with_depth = false)

    interpolate_onto!(sea_ice_model.ice_concentration, concentration)
    interpolate_onto!(sea_ice_model.ice_thickness, thickness)

    return nothing
end

"""
    tripolar_one_degree_bottom_height(grid; minimum_depth, major_basins, interpolation_passes)

The ETOPO-derived bottom height for the one-degree tripolar grid.

The regridding parameters are baked into the artifact, so they are validated rather than applied: a
mismatch means the caller is asking for bathymetry this artifact does not describe.
"""
function tripolar_one_degree_bottom_height(
    grid;
    minimum_depth,
    major_basins,
    interpolation_passes,
)
    directory = @clima_artifact("tripolar_one_degree_bathymetry")
    path = joinpath(directory, "tripolar_one_degree_bathymetry.nc")

    bottom_height, attributes = NCDatasets.NCDataset(path) do ds
        Array(ds["bottom_height"][:, :]), Dict(ds.attrib)
    end

    for (name, requested) in (
        ("minimum_depth", minimum_depth),
        ("major_basins", major_basins),
        ("interpolation_passes", interpolation_passes),
    )
        baked = attributes[name]
        baked == requested || error("""
              The tripolar bathymetry artifact was generated with $name = $baked, but $requested was \
              requested. Regenerate it with staging/climaartifacts/tripolar_one_degree_bathymetry/ \
              or use the baked value.
              """)
    end

    size(bottom_height) == (grid.Nx, grid.Ny) || error(
        "The tripolar bathymetry artifact is $(size(bottom_height)) but the grid is ($(grid.Nx), $(grid.Ny))",
    )

    bottom_field = OC.Field{OC.Center, OC.Center, Nothing}(grid)
    OC.set!(bottom_field, bottom_height)

    return bottom_field
end
