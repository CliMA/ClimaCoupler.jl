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
)

    diags_names = propertynames(post_spec)

    all_plots = []
    for name in diags_names
        @info name

        # extract data
        diag_data = read_latest_model_data(name, files_dir, files_root)

        # postprocess
        post_data = postprocess(name, diag_data, getproperty(post_spec, name))

        # create individual plots
        p = plot(
            post_data,
            zmd_params = (; getproperty(plot_spec, name)...),
            hsd_params = (; getproperty(plot_spec, name)...),
        )

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

"""
    read_latest_model_data(name::Symbol, filedir::String, root::String)

Reads in a variable from a HDF5 file, as outputted by `ClimaCoupler.Dignostics`. 
"""
function read_latest_model_data(name::Symbol, filedir::String, root::String)

    varfile_root = @sprintf "%s%s" string(name) root
    filename = glob("*" * varfile_root * "*", filedir)[end]

    hdfreader = InputOutput.HDF5Reader(filename)
    var = InputOutput.read_field(hdfreader, string(name))
    close(hdfreader)
    return var
end
