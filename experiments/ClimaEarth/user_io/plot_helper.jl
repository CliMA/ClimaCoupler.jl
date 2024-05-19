import Adapt
import CUDA
import Plots
import ClimaCoupler: PostProcessor, TestHelper

"""
    Plots.plot(post_data::DataPackage; zmd_params = (;), hsd_params = (;))

Coordinates plotting based on parsed data types.
"""
function Plots.plot(post_data::PostProcessor.DataPackage; zmd_params = (;), hsd_params = (;))
    if post_data.tag isa PostProcessor.ZLatData
        plot_params = zmd_params
    elseif post_data.tag isa PostProcessor.LatLonData
        plot_params = hsd_params
    else
        plot_params = (;)
    end
    Plots.contourf(post_data.tag, post_data; plot_params...)
end

"""
    Plots.contourf(
        ::PostProcessor.ZLatData,
        p::PostProcessor.DataPackage;
        xlabel = "lat (deg N)",
        ylabel = "z (km)",
        yaxis = (:log,),
        yflip = false,
        clims = nothing,
        units = " ",
    )

Plots a filled contour plot on the latitude-level plane.
"""
function Plots.contourf(
    ::PostProcessor.ZLatData,
    p::PostProcessor.DataPackage;
    xlabel = "lat (deg N)",
    ylabel = "z (km)",
    yaxis = (:log,),
    yflip = false,
    clims = nothing,
    units = " ",
)
    clims = isnothing(clims) ? extrema(p.data) : clims
    Plots.contourf(
        p.coords.lat,
        p.coords.lev,
        p.data',
        color = :rainbow,
        title = string(p.name) * " [" * units * "]",
        xlabel = xlabel,
        ylabel = ylabel,
        yaxis = yaxis,
        yflip = yflip,
        clims = clims,
    )
end

"""
    Plots.contourf(
        ::PostProcessor.LatLonData,
        p::PostProcessor.DataPackage;
        xlabel = "lat (deg N)",
        ylabel = "z (km)",
        yaxis = (:log,),
        yflip = false,
        clims = nothing,
        units = " ",
    )

Plots a filled contour plot on the longitude-latitude plane.
"""
function Plots.contourf(
    ::PostProcessor.LatLonData,
    p::PostProcessor.DataPackage;
    xlabel = "lon (deg E)",
    ylabel = "lat (deg N)",
    clims = nothing,
    units = " ",
)
    clims = isnothing(clims) ? extrema(p.data) : clims
    Plots.contourf(
        p.coords.lon,
        p.coords.lat,
        p.data',
        color = :rainbow,
        title = string(p.name) * " [" * units * "]",
        xlabel = xlabel,
        ylabel = ylabel,
        clims = clims,
    )
end

"""
    get_ne(grid)

Return the number of horizontal elements in a grid.
"""
function get_ne(grid::CC.Grids.SpectralElementGrid2D)
    return grid.topology.mesh.ne
end
function get_ne(grid::CC.Grids.LevelGrid)
    return get_ne(grid.full_grid.horizontal_grid)
end
function get_ne(grid::CC.Grids.ExtrudedFiniteDifferenceGrid)
    return get_ne(grid.horizontal_grid)
end

"""
    get_R(grid)

Return the radius of a grid.
"""
function get_R(grid)
    return CC.Grids.global_geometry(grid).radius
end

"""
    get_height(grid)

Return the height of a grid if it is 3D, or nothing otherwise.
"""
function get_height(grid::CC.Grids.ExtrudedFiniteDifferenceGrid)
    return grid.vertical_grid.topology.mesh.domain.coord_max.z
end
function get_height(grid)
    return nothing # 2d case
end

"""
    to_cpu(field::CC.Fields.Field)

For a CPU field, return the field unchanged.
For a GPU field, copy the field and its underlying space onto the CPU.

Note that this function allocates a new space and field,
and should only be used for debugging and testing.
"""
function to_cpu(field::CC.Fields.Field)
    if parent(field) isa Array
        return field
    else
        # Copy field onto cpu space
        space = axes(field)
        FT = CC.Spaces.undertype(space)
        R = get_R(space.grid)
        ne = get_ne(space.grid)
        polynomial_degree = CC.Quadratures.polynomial_degree(CC.Spaces.quadrature_style(space.grid))
        nz = CC.Spaces.nlevels(space)
        height = get_height(space.grid)

        cpu_comms_ctx = ClimaComms.SingletonCommsContext(ClimaComms.CPUSingleThreaded())
        cpu_space = TestHelper.create_space(
            FT,
            comms_ctx = cpu_comms_ctx,
            R = R,
            ne = ne,
            polynomial_degree = polynomial_degree,
            nz = nz,
            height = height,
        )
        cpu_field = CC.Fields.ones(cpu_space)

        parent(cpu_field) .= Array(parent(field))
        return cpu_field
    end
end

to_cpu(arr::Array) = arr
to_cpu(arr::CUDA.CuArray) = Adapt.adapt(Array, arr)
