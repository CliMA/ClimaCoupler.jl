using Downloads
using NCDatasets
import ClimaCoupler.Diagnostics: get_var
import ClimaCoupler.PostProcessor: postprocess

include("plot_helper.jl")

"""
    ncep_paperplots(
        post_spec::NamedTuple,
        plot_spec::NamedTuple,
        files_dir::String;
        month_date = Dates.DateTime(1979, 01, 01),
        output_dir = ".",
        fig_name = "ncep_paperplots",
    )

Coordinates the postprocessing and plotting of sample fields (specified in `post_spec`)
of a particular monthly mean dataset (specified by `month_date`). Any plot NCEP- specific
customization should be done here.
"""
function ncep_paperplots(
    post_spec::NamedTuple,
    plot_spec::NamedTuple,
    files_dir::String;
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
    all_data = (;)
    for vname in diags_vnames
        @info vname

        # download and read data of this month
        data, coords = get_var(ncep_src, Val(vname))
        raw_tag = length(size(data)) == 3 ? ZLatLonData() : LatLonData()

        # postprocess
        post_data = postprocess(vname, data, getproperty(post_spec, vname), coords = coords, raw_tag = raw_tag)

        # create individual plots
        zonal_mean_params = (; ylabel = "p (hPa)", yaxis = (;), yflip = true, getproperty(plot_spec, vname)...) # overwrite defaults
        p = plot(post_data, zmd_params = zonal_mean_params, hsd_params = (; getproperty(plot_spec, vname)...))

        push!(all_plots, p)

        # create a named tuple with data
        data = post_data.data
        all_data = merge(all_data, [vname => data])
    end

    # combine plots and save figure
    save_fig = Plots.plot(
        all_plots...,
        size = (1500, 1200),
        right_margin = 3Plots.mm,
        left_margin = 3Plots.mm,
        bottom_margin = 3Plots.mm,
        top_margin = 3Plots.mm,
    )

    Plots.png(save_fig, joinpath(output_dir, fig_name * ".png"))

    return all_data
end

abstract type DataSource end
struct NCEPMonthlyDataSource <: DataSource
    tmp_dir::String
    month_date::Array
end

"""
    download_read_nc(data_source::NCEPMonthlyDataSource, https::String, ncep_vname::String)

Downloads and reads nc datafile of a specified NCEP variable
"""
function download_read_nc(data_source::NCEPMonthlyDataSource, https::String, ncep_vname::String)
    local_file = joinpath(data_source.tmp_dir, ncep_vname * ".nc")
    Downloads.download(https, local_file)
    NCDataset(local_file) do ds
        t_i = findall(x -> Dates.yearmonth(x) == Dates.yearmonth(data_source.month_date[1]), Array(ds["time"])) # time index of month in file
        d_i = length(size(Array(ds[ncep_vname]))) # index of time in the dimension list
        lev = "level" in keys(ds) ? Array(ds["level"]) : [Float64(-999)]
        data = dropdims(selectdim(Array(ds[ncep_vname]), d_i, t_i), dims = d_i)
        if length(size(data)) == 3
            data .= data[:, end:-1:1, end:-1:1]
        else
            data .= data[:, end:-1:1]
        end
        coords = (; lon = Array(ds["lon"]), lat = Array(ds["lat"])[end:-1:1], lev = lev[end:-1:1])
        (Array(data), coords)
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

function get_var(data_source::NCEPMonthlyDataSource, ::Val{:toa_fluxes})
    https_root = "https://downloads.psl.noaa.gov/Datasets/ncep.reanalysis2/Monthlies/gaussian_grid/"
    https_suffix = ".ntat.mon.mean.nc"

    # sum, with upward positive. NB: dlwrf not available
    # toa_tuple = download_read_nc(data_source, https_root * "uswrf" * https_suffix, "uswrf")
    # toa_tuple.data .-= download_read_nc(data_source, https_root * "dswrf" * https_suffix, "dswrf").data
    # toa_tuple.data .+= download_read_nc(data_source, https_root * "ulwrf" * https_suffix, "ulwrf").data

    data_uswrf, coords = download_read_nc(data_source, https_root * "uswrf" * https_suffix, "uswrf")
    data_dswrf, _ = download_read_nc(data_source, https_root * "dswrf" * https_suffix, "dswrf")
    data_ulwrf, _ = download_read_nc(data_source, https_root * "ulwrf" * https_suffix, "ulwrf")

    return (data_uswrf .- data_dswrf .+ data_ulwrf, coords)
end

function get_var(data_source::NCEPMonthlyDataSource, ::Val{:precipitation_rate})
    https = "https://downloads.psl.noaa.gov/Datasets/ncep.reanalysis2/Monthlies/gaussian_grid/prate.sfc.mon.mean.nc"
    ncep_vname = "prate"
    download_read_nc(data_source, https, ncep_vname)
end

function get_var(data_source::NCEPMonthlyDataSource, ::Val{:T_sfc})
    https = "https://downloads.psl.noaa.gov/Datasets/ncep.reanalysis/Monthlies/surface/air.sig995.mon.mean.nc" # NB: strictly speaking this is air temperature not T_sfc
    ncep_vname = "air"
    t_celsius, coords = download_read_nc(data_source, https, ncep_vname)
    return (t_celsius .+ 273.15, coords)
end

function get_var(data_source::NCEPMonthlyDataSource, ::Val{:tubulent_energy_fluxes})
    https = "https://downloads.psl.noaa.gov/Datasets/ncep.reanalysis2/Monthlies/gaussian_grid/lhtfl.sfc.mon.mean.nc"
    ncep_vname = "lhtfl"
    lhtfl, coords = download_read_nc(data_source, https, ncep_vname)
    https = "https://downloads.psl.noaa.gov/Datasets/ncep.reanalysis2/Monthlies/gaussian_grid/shtfl.sfc.mon.mean.nc"
    ncep_vname = "shtfl"
    shtfl, coords = download_read_nc(data_source, https, ncep_vname)
    return (shtfl .+ lhtfl, coords)
end
