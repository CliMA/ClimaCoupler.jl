import ClimaAnalysis
import ClimaAnalysis: OutputVar
import ClimaCalibrate
import Dates
import EnsembleKalmanProcesses as EKP
import JLD2
import ClimaDiagnostics
import ClimaCore
# Access CalibrateConfig
include(joinpath(@__DIR__, "run_calibration.jl"))

"""
    load_vars(obsdir, short_names)

Load NetCDF files belonging to `short_names` in `obsdir` as `OutputVar`s.
"""
function load_vars()
    # TODO: Use era5_monthly_averages_surface_single_level_1979_2024 artifact

    ta_900hpa = OutputVar("/glade/u/home/nefrathe/clima/ClimaCoupler.jl/era5_oct_2017_2024_ta_900hpa.nc")
    ta_900hpa = reverse_dim(ta_900hpa, latitude_name(ta_900hpa))

    tas = OutputVar("/glade/u/home/nefrathe/clima/ClimaCoupler.jl/era5_oct_2017_2024_tas.nc")
    tas = reverse_dim(tas, latitude_name(tas))
    tas.attributes["short_name"] = "tas"

    ta_900hpa = slice(ta_900hpa; pressure = 900)

    tas_minus_ta_900hpa = tas - ta_900hpa
    tas_minus_ta_900hpa.attributes["short_name"] = "tas - ta"

    lhf = OutputVar("era5_oct_2017_2024_surface_fluxes.nc", "avg_slhtf")
    lhf = reverse_dim(lhf, latitude_name(lhf))
    lhf = -lhf # Flip to positive upward to match simulation
    lhf.attributes["short_name"] = "hfls"

    lhf = ClimaAnalysis.Var._shift_by(lhf, date -> date - Dates.Hour(6))

    shf = OutputVar("era5_oct_2017_2024_surface_fluxes.nc", "avg_ishf")
    shf = reverse_dim(shf, latitude_name(shf))
    shf = -shf # Flip to positive upward to match simulation
    shf.attributes["short_name"] = "hfss"
    shf = ClimaAnalysis.Var._shift_by(shf, date -> date - Dates.Hour(6))

    rsns = OutputVar("era5_oct_2017_2024_surface_fluxes.nc", "avg_snswrf")
    rsns = reverse_dim(rsns, latitude_name(rsns))
    rsns = -rsns # Flip to positive upward to match simulation
    rsns.attributes["short_name"] = "rsns"
    rsns = ClimaAnalysis.Var._shift_by(rsns, date -> date - Dates.Hour(6))

    rlns = OutputVar("era5_oct_2017_2024_surface_fluxes.nc", "avg_snlwrf")
    rlns = reverse_dim(rlns, latitude_name(ta_900hpa))
    rlns = -rlns # Flip to positive upward to match simulation
    rlns.attributes["short_name"] = "rlns"
    rlns = ClimaAnalysis.Var._shift_by(rlns, date -> date - Dates.Hour(6))

    return [tas_minus_ta_900hpa, tas, lhf, shf, rsns, rlns]
end

"""
    preprocess_vars(vars)

Preprocess each OutputVar in `vars` by keeping the relevant dates in
`sample_date_ranges`.
"""
function preprocess_vars(vars)
    diagnostic_var = "/glade/derecho/scratch/nefrathe/tmp/output_quick_1/iteration_000/member_001/wxquest_diagedmf/output_active/clima_atmos/mslp_1week_average.nc"
    if isfile(diagnostic_var)
        out_var = OutputVar(diagnostic_var)
        resample_var(x) = ClimaAnalysis.resampled_as(x, out_var ; dim_names = ["lon", "lat"])
    else
        resample_var(x) = resampled_lonlat(CALIBRATE_CONFIG.config_file)
    end

    vars = map(vars) do var
        var = resample_var(var)
        var = set_units(var, var_units[short_name(var)])
        var = apply_landmask(var)
        remove_global_mean(var)
    end

    return vars
end

var_units = Dict(
    "pr" => "kg m^-2 s^-1",
    "mslp" => "Pa",
    "tas" => "K",
    "tas - ta" => "K",
    "hfls" => "W m^-2",
    "hfss" => "W m^-2",
    "rsns" => "W m^-2",
    "rlns" => "W m^-2",
    )


function make_observation_vector(vars, sample_date_ranges)
    covar_estimator =  ClimaCalibrate.ObservationRecipe.ScalarCovariance(; scalar = 5.0, use_latitude_weights = true)
    obs_vec = map(sample_date_ranges) do sample_date_range
        ClimaCalibrate.ObservationRecipe.observation(
            covar_estimator,
            vars,
            first(sample_date_range),
            last(sample_date_range),
        )
    end
    return obs_vec
end

"""
    resampled_lonlat(config_file)

Return a function to resample longitude and latitudes according to the model
grid specified by `config_file`.
"""
function resampled_lonlat(config_file)
    config_dict = get_coupler_config_dict(CALIBRATE_CONFIG.config_file)
    if !isnothing(config_dict["netcdf_interpolation_num_points"])
        (nlon, nlat, nlev) = tuple(config_dict["netcdf_interpolation_num_points"]...)
    else
        cs = CoupledSimulation(config_file)
        center_space = cs.model_sims.atmos_sim.domain.center_space
        (nlon, nlat, nlev) =
            ClimaDiagnostics.Writers.default_num_points(center_space)
        stretch = center_space.grid.vertical_grid.topology.mesh.stretch
        dz_bottom = center_space.grid.vertical_grid.topology.mesh.faces[2].z
        z = range(dz_bottom, ClimaCore.Spaces.z_max(center_space), nlev)
    end
    lon = range(-180, 180, nlon)
    lat = range(-90, 90, nlat)
    # TODO: Account for stretch for 3D variables
    return var -> resampled_as(var; lon, lat)
end

if abspath(PROGRAM_FILE) == @__FILE__
    sample_date_ranges = CALIBRATE_CONFIG.sample_date_ranges
    short_names = CALIBRATE_CONFIG.short_names
    config_file = CALIBRATE_CONFIG.config_file
    @info "Generating observations for $short_names"
    @info "The number of samples is $(length(sample_date_ranges)) over $sample_date_ranges"

    unprocessed_vars = load_vars()
    preprocessed_vars = preprocess_vars(unprocessed_vars)

    JLD2.save_object(
        joinpath(pkgdir(ClimaCoupler),"experiments/calibration/era5_preprocessed_vars.jld2"),
        preprocessed_vars,
    )
    observation_vector = make_observation_vector(preprocessed_vars, sample_date_ranges)
    JLD2.save_object(
        joinpath(pkgdir(ClimaCoupler),"experiments/calibration/weatherquest_obs_vec.jld2"),
        observation_vector,
    )
end
