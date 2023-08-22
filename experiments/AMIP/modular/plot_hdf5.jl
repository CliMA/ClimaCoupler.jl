# This is a custom script


import SciMLBase: step!, reinit!
using OrdinaryDiffEq
using OrdinaryDiffEq: ODEProblem, solve, SSPRK33, savevalues!, Euler
using LinearAlgebra
import Test: @test
using Dates
using UnPack
using Plots
using Statistics: mean

using ClimaCore.Utilities: half, PlusHalf
using ClimaCore: InputOutput, Fields
import ClimaCore.Spaces as Spaces

import ClimaCoupler
import ClimaCoupler.Regridder
import ClimaCoupler.Regridder:
    update_surface_fractions!, combine_surfaces!, combine_surfaces_from_sol!, dummmy_remap!, binary_mask
import ClimaCoupler.ConservationChecker:
    EnergyConservationCheck, WaterConservationCheck, check_conservation!, plot_global_conservation
import ClimaCoupler.Utilities: CoupledSimulation, float_type, swap_space!
import ClimaCoupler.BCReader:
    bcfile_info_init, float_type_bcf, update_midmonth_data!, next_date_in_file, interpolate_midmonth_to_daily
import ClimaCoupler.TimeManager: current_date, datetime_to_strdate, trigger_callback, Monthly, EveryTimestep
import ClimaCoupler.Diagnostics: get_var, init_diagnostics, accumulate_diagnostics!, save_diagnostics, TimeMean
import ClimaCoupler.PostProcessor: postprocess

import ClimaCoupler.Interfacer:
    AtmosModelSimulation,
    SurfaceModelSimulation,
    SurfaceStub,
    SeaIceModelSimulation,
    LandModelSimulation,
    OceanModelSimulation,
    get_field,
    update_field!,
    update_sim!
import ClimaCoupler.FluxCalculator:
    PartitionedStateFluxes,
    CombinedStateFluxes,
    combined_turbulent_fluxes!,
    MoninObukhovScheme,
    partitioned_turbulent_fluxes!
import ClimaCoupler.FieldExchanger:
    import_atmos_fields!,
    import_combined_surface_fields!,
    update_sim!,
    update_model_sims!,
    reinit_model_sims!,
    step_model_sims!
import ClimaCoupler.Checkpointer: checkpoint_model_state, get_model_state_vector, restart_model_state!

include("user_io/user_diagnostics.jl")

HOME_DIR = "experiments/AMIP/modular/"
COUPLER_OUTPUT_DIR = HOME_DIR*"output/amip/coarse_single_modular_businger_ft64"
COUPLER_ARTIFACTS_DIR = COUPLER_OUTPUT_DIR*"_artifacts"

@info COUPLER_OUTPUT_DIR
@info COUPLER_ARTIFACTS_DIR

include("user_io/plot_helper.jl")
## ClimaESM
@info "AMIP plots"
include("user_io/amip_visualizer.jl")
amip_post_spec = (;
    T = (:regrid, :zonal_mean),
    u = (:regrid, :zonal_mean),
    q_tot = (:regrid, :zonal_mean),
    toa  = (:regrid, :horizontal_slice),
    precipitation  = (:regrid, :horizontal_slice),
    T_sfc = (:regrid, :horizontal_slice),
    )

amip_plot_spec = (;
    T = (; clims = (190, 320), units = "K"),
    u = (; clims = (-50, 50), units = "m/s"),
    q_tot = (; clims = (0, 30), units = "g/kg"),
    toa = (; clims = (-250, 250), units = "W/m^2"),
    precipitation = (clims = (0, 1e-4), units = "kg/m^2/s"),
    T_sfc = (clims = (225, 310), units = "K"),
    )

# amip_data = amip_paperplots(
#     amip_post_spec,
#     amip_plot_spec,
#     COUPLER_OUTPUT_DIR,
#     files_root = ".monthly",
#     output_dir = COUPLER_ARTIFACTS_DIR,
#     )


## NCEP reanalysis
@info "NCEP plots"
include("user_io/ncep_visualizer.jl")
ncep_post_spec = (;
    T = (:zonal_mean,),
    u = (:zonal_mean,),
    q_tot = (:zonal_mean,),
    toa = (:horizontal_slice,),
    precipitation = (:horizontal_slice,),
    T_sfc = (:horizontal_slice,),
)

ncep_plot_spec = plot_spec
# ncep_data = ncep_paperplots(
#     ncep_post_spec,
#     ncep_plot_spec,
#     COUPLER_OUTPUT_DIR,
#     output_dir = COUPLER_ARTIFACTS_DIR,
#     month_date = Dates.DateTime(1979, 01, 01),
# ) ## plot data that correspond to the model's last save_hdf5 call (i.e., last month)



function diff_paperplots(;
    files_dir = COUPLER_OUTPUT_DIR,
    ncep_post_spec = ncep_post_spec,
    ncep_plot_spec = ncep_plot_spec,
    ncep_month_date = Dates.DateTime(1979, 02, 01),
    amip_post_spec = amip_post_spec,
    amip_plot_spec = amip_plot_spec,
    output_dir = COUPLER_ARTIFACTS_DIR,
    files_root = ".monthly",
    nlat = 180,
    nlon = 360,
    fig_name = "diff_paperplots",
    )

    # NCEP Data extraction
    tmp_dir = joinpath(files_dir, "ncep_tmp")
    isdir(tmp_dir) ? nothing : mkpath(tmp_dir)

    ncep_src = NCEPMonthlyDataSource(tmp_dir, [ncep_month_date])
    diags_vnames = propertynames(ncep_post_spec)

    # AMIP Data extraction
    diags_names = propertynames(amip_post_spec)

    all_plots = []
    all_data = (;)
    for name in diags_names, vname in diags_vnames

        @info name
        @info vname

        # AMIP: extract data
        diag_data = read_latest_model_data(name, files_dir, files_root)

        # NCEP: download and read data of this month
        data, coords = get_var(ncep_src, Val(vname))
        raw_tag = length(size(data)) == 3 ? ZLatLonData() : LatLonData()

        # AMIP: post processes
        amip_post_data = postprocess(name, diag_data, getproperty(amip_post_spec, name), nlat = nlat, nlon = nlon)
        amip_post_data.data[1] = sum(amip_post_data.data) == 0 ? amip_post_data.data[1] + eps() : amip_post_data.data[1] # avoids InexactError

        # NCEP: post processes
        ncep_post_data = postprocess(vname, data, getproperty(ncep_post_spec, vname), coords = coords, raw_tag = raw_tag)
        # @show ncep_post_data
        # @show amip_post_data

        # NCEP - AMIP Difference
        post_data = @. ncep_post_data.data - amip_post_data.data

        # AMIP/NCEP create individual plots
        p = plot(
            post_data,
            zmd_params = (; getproperty(amip_plot_spec, name)...),
            hsd_params = (; getproperty(amip_plot_spec, name)...),
        )

        push!(all_plots, p)

        # create a named tuple with data
        data = post_data.data
        all_data = merge(all_data, [name => data])
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

diff_paperplots()



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
        t_i = findall(x -> Dates.yearmonth(x) == Dates.yearmonth(data_source.month_date[1]), ds["time"][:]) # time index of month in file
        d_i = length(size(ds[ncep_vname][:])) # index of time in the dimension list
        lev = "level" in keys(ds) ? ds["level"][:] : [Float64(-999)]
        data = dropdims(selectdim(ds[ncep_vname][:], d_i, t_i), dims = d_i)
        if length(size(data)) == 3
            data .= data[:, end:-1:1, end:-1:1]
        else
            data .= data[:, end:-1:1]
        end
        coords = (; lon = ds["lon"][:], lat = ds["lat"][:][end:-1:1], lev = lev[end:-1:1])
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

function get_var(data_source::NCEPMonthlyDataSource, ::Val{:toa})
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

function get_var(data_source::NCEPMonthlyDataSource, ::Val{:precipitation})
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
