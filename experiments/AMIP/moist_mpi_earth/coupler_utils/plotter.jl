using Glob
using Printf

function plot(post_data::PostProcessedData; zmd_params = (;), hsd_params = (;))

    if post_data isa ZonalMeanData
        plot_params = zmd_params
    elseif post_data isa HorizontalSliceData
        plot_params = hsd_params
    else
        plot_params = (;)
    end
    contourf(post_data; plot_params...)
end

function contourf(
    p::ZonalMeanData;
    xlabel = "lat (deg N)",
    ylabel = "z (km)",
    yaxis = (:log,),
    yflip = false,
    clims = nothing,
    units = " ",
)
    clims = clims == nothing ? extrema(p.data_2d) : clims
    plot_p = Plots.contourf(
        p.lat,
        p.lev,
        p.data_2d',
        color = :rainbow,
        title = string(p.name) * " [" * units * "]",
        xlabel = xlabel,
        ylabel = ylabel,
        yaxis = yaxis,
        yflip = yflip,
        clims = clims,
    )
end

function contourf(p::HorizontalSliceData; xlabel = "lon (deg E)", ylabel = "lat (deg N)", clims = nothing, units = " ")
    clims = clims == nothing ? extrema(p.data_2d) : clims
    plot_p = Plots.contourf(
        p.lon,
        p.lat,
        p.data_2d',
        color = :rainbow,
        title = string(p.name) * " [" * units * "]",
        xlabel = xlabel,
        ylabel = ylabel,
        clims = clims,
    )
end
