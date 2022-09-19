using Glob
using Printf

function plot(post_data::PostProcessedData, name)
    info = post_data.info
    if info.slice_type == :zonal_mean
        p = contourf(ZonalMeanData(post_data, string(name)))
    elseif info.slice_type == :horizontal_2d
        p = contourf(HorizonalData(post_data, string(name)))
    end
    return p
end

function contourf(p::ZonalMeanData)
    plot_p = Plots.contourf(
        p.lat,
        p.lev,
        p.data_2d',
        color = :rainbow,
        title = string(p.title),
        xlabel = "lat (deg N)",
        ylabel = "z (km)",
        yaxis = (:log,),
        yflip = true,
    )
end

function contourf(p::HorizonalData)
    plot_p = Plots.contourf(
        p.lon,
        p.lat,
        p.data_2d',
        color = :rainbow,
        title = string(p.title),
        xlabel = "lon (deg N)",
        ylabel = "lat (deg N)",
    )
end

function amip_paperplots(plot_spec, files_dir; output_dir = ".", files_root = ".hdf5", fig_name = "amip_paperplots")

    diags_names = propertynames(plot_spec)

    all_plots = []
    for name in diags_names

        # extract data
        diag_data = read_latest_model_data(name, files_dir, files_root)

        # postprocess
        post_data = postprocess(name, diag_data, getproperty(plot_spec, name))

        # create individual plots
        p = plot(post_data, name)

        push!(all_plots, p)
    end

    # combine plots and save figure
    save_fig = Plots.plot(
        all_plots...,
        size = (1500, 800),
        right_margin = 12Plots.mm,
        left_margin = 12Plots.mm,
        bottom_margin = 12Plots.mm,
        top_margin = 12Plots.mm,
    )

    Plots.png(save_fig, joinpath(output_dir, fig_name * ".png"))

    return all_plots
end

function read_latest_model_data(name, filedir, root)
    @show name
    varfile_root = @sprintf "%s%s" string(name) root
    filename = glob("*" * varfile_root * "*", filedir)[end]

    hdfreader = InputOutput.HDF5Reader(filename)
    hdf5_data = InputOutput.read_field(hdfreader, string(name))
end

print_formatted(fmt, args...) = @eval @printf($fmt, $(args...))
