# using Plots
# import ClimaCorePlots
# import ClimaCore
# using ClimaCorePlots
# using ClimaCore.InputOutput # https://clima.github.io/ClimaCore.jl/dev/api/#InputOutput
# datafile = ClimaCore.InputOutput.HDF5Reader("experiments/AMIP/modular/output/amip/coarse_single_modular_cheng_ft64/u.monthly_mean_3d_.1979-02-01T00:00:00.hdf5")
# variable = ClimaCore.InputOutput.read_field(datafile, "u")
# u |> propertynames # Explore variables available in HDF5 Dataset
# Plots.plot(variable)  # Plots on an unwrapped cubed-sphere 

# ## TODO: reconfigure plots from amip pipeline run


## sample animations
if !is_distributed && parsed_args["anim"]
    @info "Animations"
    include("user_io/viz_explorer.jl")
    plot_anim(cs, COUPLER_ARTIFACTS_DIR)
end

## plotting AMIP results
@info "AMIP plots"


include("experiments/AMIP/modular/user_io/plot_helper.jl")

"""
    function diff_paperplots(
        post_spec::NamedTuple,
        plot_spec::NamedTuple,
        files_dir::String;
        output_dir = ".",
        files_root = ".hdf5",
        fig_name = "diff_paperplots",
    )
Coordinates the postprocessing and plotting of sample fields (specified in `post_spec`)
of the last monthly mean file. Any specific plot customization should be done here.
"""
function diff_paperplots(
    post_spec1::NamedTuple,
    post_spec2::NamedTuple
    plot_spec::NamedTuple,
    files_dir1::String
    files_dir2::String;
    output_dir = ".",
    files_root = ".hdf5",
    fig_name = "diff_paperplots",
)

    diags_names = propertynames(post_spec)

    all_plots = []
    for name in diags_names
        @info name

        # extract data
        diag_data1 = read_latest_model_data(name, files_dir1, files_root)
        diag_data2 = read_latest_model_data(name, files_dir2, files_root)

        # postprocess
        post_data1 = postprocess(name, diag_data1, getproperty(post_spec, name))
        post_data1.data[1] = sum(post_data1.data) == 0 ? post_data1.data[1] + eps() : post_data1.data[1] # avoids InexactError

        post_data2 = postprocess(name, diag_data2, getproperty(post_spec, name))
        post_data2.data[1] = sum(post_data2.data) == 0 ? post_data2.data[1] + eps() : post_data2.data[1] # avoids InexactError

        diff_data = @. post_data2 - post_data1

        # create individual plots
        p = plot(
            diff_data,
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


post_spec = (;
    T = (:regrid, :zonal_mean),
    u = (:regrid, :zonal_mean),
    q_tot = (:regrid, :zonal_mean),
    toa = (:regrid, :horizontal_slice),
    precipitation = (:regrid, :horizontal_slice),
    T_sfc = (:regrid, :horizontal_slice),
)

plot_spec = (;
    T = (; clims = (190, 320), units = "K"),
    u = (; clims = (-50, 50), units = "m/s"),
    q_tot = (; clims = (0, 50), units = "g/kg"),
    toa = (; clims = (-250, 210), units = "W/m^2"),
    precipitation = (clims = (0, 1e-6), units = "kg/m^2/s"),
    T_sfc = (clims = (225, 310), units = "K"),
)

amip_paperplots(
    post_spec,
    plot_spec,
    COUPLER_OUTPUT_DIR,
    files_root = ".monthly",
    output_dir = COUPLER_ARTIFACTS_DIR,
)

## NCEP reanalysis
ncep_post_spec = (;
    T = (:zonal_mean,),
    u = (:zonal_mean,),
    q_tot = (:zonal_mean,),
    toa = (:horizontal_slice,),
    precipitation = (:horizontal_slice,),
    T_sfc = (:horizontal_slice,),
)
ncep_plot_spec = plot_spec
ncep_paperplots(
    ncep_post_spec,
    ncep_plot_spec,
    COUPLER_OUTPUT_DIR,
    output_dir = COUPLER_ARTIFACTS_DIR,
    month_date = cs.dates.date[1],
) ## plot data that correspond to the model's last save_hdf5 call (i.e., last month)



