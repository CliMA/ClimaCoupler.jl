import Plots
import ClimaCoupler: PostProcessor

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
