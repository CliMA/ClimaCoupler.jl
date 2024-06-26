import Glob
import Plots
import Printf
import ClimaComms
import ClimaCore as CC
import ClimaCoupler: PostProcessor

include("plot_helper.jl")

"""
    function amip_paperplots(
        post_spec::NamedTuple,
        plot_spec::NamedTuple,
        files_dir::String;
        output_dir = ".",
        files_root = ".hdf5",
        fig_name = "amip_paperplots",
    )

Coordinates the postprocessing and plotting of sample fields (specified in `post_spec`)
of the last monthly mean file. Any specific plot customization should be done here.
"""
function amip_paperplots(
    post_spec::NamedTuple,
    plot_spec::NamedTuple,
    files_dir::String;
    output_dir = ".",
    files_root = ".hdf5",
    fig_name = "amip_paperplots",
    nlat = 180,
    nlon = 360,
)
    diags_names = propertynames(post_spec)

    all_plots = []
    all_data = (;)
    for name in diags_names
        @info name

        # extract data
        diag_data = read_latest_model_data(name, files_dir, files_root)

        # postprocess
        post_data = PostProcessor.postprocess(
            name,
            diag_data,
            getproperty(post_spec, name),
            REGRID_DIR = files_dir,
            nlat = nlat,
            nlon = nlon,
        )
        post_data.data[1] = sum(post_data.data) == 0 ? post_data.data[1] + eps() : post_data.data[1] # avoids InexactError

        # create individual plots
        p = Plots.plot(
            post_data,
            zmd_params = (; getproperty(plot_spec, name)...),
            hsd_params = (; getproperty(plot_spec, name)...),
        )

        push!(all_plots, p)

        # create a named tuple with data
        data = post_data.data
        all_data = merge(all_data, [name => data])
    end

    # combine plots and save figure
    layout = Plots.@layout([A{0.05h}; [B C D; E F G; H I J]])
    title =
        Plots.plot(plot_title = "AMIP Monthly Mean Fields", grid = false, showaxis = false, bottom_margin = -50Plots.px)
    save_fig = Plots.plot(
        title,
        all_plots...,
        size = (1500, 1200),
        layout = layout,
        titlefont = Plots.font(12),
        margin = 2Plots.mm,
    )

    Plots.png(save_fig, joinpath(output_dir, fig_name * ".png"))

    return all_data, save_fig
end

"""
    read_latest_model_data(name::Symbol, filedir::String, root::String)

Reads in a variable from a HDF5 file, as outputted by `ClimaCoupler.Dignostics`.
"""
function read_latest_model_data(name::Symbol, filedir::String, root::String)
    varfile_root = Printf.@sprintf "%s%s" string(name) root
    filename = Glob.glob("*" * varfile_root * "*", filedir)[end]

    # Ensure file gets read onto CPU for postprocessing
    cpu_singleton_context = ClimaComms.SingletonCommsContext(ClimaComms.CPUSingleThreaded())
    hdfreader = CC.InputOutput.HDF5Reader(filename, cpu_singleton_context)
    var = CC.InputOutput.read_field(hdfreader, string(name))
    close(hdfreader)
    return var
end
