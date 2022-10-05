"""
    ncep_paperplots(
        post_spec,
        plot_spec,
        files_dir;
        month_date = Dates.DateTime(1979, 01, 01),
        output_dir = ".",
        fig_name = "ncep_paperplots",
    )

Coordinates the postprocessing and plotting of sample fields (specified in `post_spec`) 
of a particular monthly mean dataset (specified by `month_date`). Any plot NCEP- specific
customization should be done here. 
"""
function ncep_paperplots(
    post_spec,
    plot_spec,
    files_dir;
    month_date = Dates.DateTime(1979, 01, 01),
    output_dir = ".",
    fig_name = "ncep_paperplots",
)

    @info month_date

    tmp_dir = joinpath(files_dir, "ncep_tmp")
    isdir(tmp_dir) ? nothing : mkpath(tmp_dir)

    ncep_src = NCEPMonthlyDataSource(tmp_dir, [month_date])
    diags_vnames = propertynames(post_spec)

    all_plots = []
    for vname in diags_vnames
        @info vname

        # download and read data of this month
        diag_data = get_var(ncep_src, Val(vname))

        # postprocess
        post_data = postprocess(vname, diag_data, getproperty(post_spec, vname))

        # create individual plots
        zonal_mean_params = (; ylabel = "p (hPa)", yaxis = (;), yflip = true, getproperty(plot_spec, vname)...) # overwrite defaults
        p = plot(post_data, zmd_params = zonal_mean_params, hsd_params = (; getproperty(plot_spec, vname)...))

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

using Downloads
using NCDatasets
abstract type DataSource end
struct NCEPMonthlyDataSource <: DataSource
    tmp_dir::String
    month_date::Array
end

"""
    download_read_nc(https, tmp_dir, ncep_vname )

Downloads and reads nc datafile of a specified NCEP variable
"""
function download_read_nc(data_source::NCEPMonthlyDataSource, https, ncep_vname)
    local_file = joinpath(data_source.tmp_dir, ncep_vname * ".nc")
    Downloads.download(https, local_file)
    NCDataset(local_file) do ds
        t_i = findall(x -> Dates.yearmonth(x) == Dates.yearmonth(data_source.month_date[1]), ds["time"][:]) # time index of month in file
        d_i = length(size(ds[ncep_vname][:])) # index of time in the dimension list
        lev = "level" in keys(ds) ? ds["level"][:] : [Float64(-999)]
        data = dropdims(selectdim(ds[ncep_vname][:], d_i, t_i), dims = d_i)
        if length(size(data)) == 3
            data .= data[:, end:-1:1, end:-1:1]
        else
            data .= data[:, end:-1:1]
        end
        (; lon = ds["lon"][:], lat = ds["lat"][:][end:-1:1], lev = lev[end:-1:1], data = data)
    end
end

# specification of variables
function get_var(data_source::NCEPMonthlyDataSource, ::Val{:T})
    https = "https://downloads.psl.noaa.gov/Datasets/ncep.reanalysis2/Monthlies/pressure/air.mon.mean.nc"
    ncep_vname = "air"
    download_read_nc(data_source, https, ncep_vname)

end

function get_var(data_source::NCEPMonthlyDataSource, ::Val{:u})
    https = "https://downloads.psl.noaa.gov/Datasets/ncep.reanalysis2/Monthlies/pressure/uwnd.mon.mean.nc"
    ncep_vname = "uwnd"
    download_read_nc(data_source, https, ncep_vname)
end

function get_var(data_source::NCEPMonthlyDataSource, ::Val{:q_tot})
    https = "https://downloads.psl.noaa.gov/Datasets/ncep.reanalysis/Monthlies/pressure/shum.mon.mean.nc" # Note that not ncep.reanalysis2
    ncep_vname = "shum"
    download_read_nc(data_source, https, ncep_vname)
end

function get_var(data_source::NCEPMonthlyDataSource, ::Val{:toa})
    https_root = "https://downloads.psl.noaa.gov/Datasets/ncep.reanalysis2/Monthlies/gaussian_grid/"
    https_suffix = ".ntat.mon.mean.nc"

    # sum, with upward positive. NB: dlwrf not available
    toa_tuple = download_read_nc(data_source, https_root * "uswrf" * https_suffix, "uswrf")
    toa_tuple.data .-= download_read_nc(data_source, https_root * "dswrf" * https_suffix, "dswrf").data
    toa_tuple.data .+= download_read_nc(data_source, https_root * "ulwrf" * https_suffix, "ulwrf").data

    return toa_tuple
end

function get_var(data_source::NCEPMonthlyDataSource, ::Val{:precipitation})
    https = "https://downloads.psl.noaa.gov/Datasets/ncep.reanalysis2/Monthlies/gaussian_grid/prate.sfc.mon.mean.nc"
    ncep_vname = "prate"
    download_read_nc(data_source, https, ncep_vname)
end

function get_var(data_source::NCEPMonthlyDataSource, ::Val{:T_sfc})
    https = "https://downloads.psl.noaa.gov/Datasets/ncep.reanalysis/Monthlies/surface/air.sig995.mon.mean.nc" # NB: strictly speaking this is air temperature not T_sfc
    ncep_vname = "air"
    t_tuple = download_read_nc(data_source, https, ncep_vname)
    t_tuple.data .= t_tuple.data .+ 273.15
    return t_tuple
end
