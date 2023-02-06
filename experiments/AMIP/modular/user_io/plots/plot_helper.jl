using Glob
using Printf
using ClimaCoupler.PostProcessor: PostProcessedData, ZLatData, LatLonData, DataPackage, ZLatLonData
using Plots

"""
    plot(post_data::DataPackage; zmd_params = (;), hsd_params = (;))

Coordinates plotting based on parsed data types.
"""
function plot(post_data::DataPackage; zmd_params = (;), hsd_params = (;))

    if post_data isa ZLatData
        plot_params = zmd_params
    elseif post_data isa LatLonData
        plot_params = hsd_params
    else
        plot_params = (;)
    end
    contourf(post_data.tag, post_data; plot_params...)
end

"""
    function contourf(
        ::ZLatData,
        p::DataPackage;
        xlabel = "lat (deg N)",
        ylabel = "z (km)",
        yaxis = (:log,),
        yflip = false,
        clims = nothing,
        units = " ",
    )

Plots a filled contour plot on the latitude-level plane.
"""
function contourf(
    ::ZLatData,
    p::DataPackage;
    xlabel = "lat (deg N)",
    ylabel = "z (km)",
    yaxis = (:log,),
    yflip = false,
    clims = nothing,
    units = " ",
)
    clims = clims == nothing ? extrema(p.data) : clims
    plot_p = Plots.contourf(
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
    function contourf(
        ::LatLonData,
        p::DataPackage;
        xlabel = "lat (deg N)",
        ylabel = "z (km)",
        yaxis = (:log,),
        yflip = false,
        clims = nothing,
        units = " ",
    )

Plots a filled contour plot on the longitude-latitude plane.
"""
function contourf(
    ::LatLonData,
    p::DataPackage;
    xlabel = "lon (deg E)",
    ylabel = "lat (deg N)",
    clims = nothing,
    units = " ",
)
    clims = clims == nothing ? extrema(p.data) : clims
    plot_p = Plots.contourf(
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
