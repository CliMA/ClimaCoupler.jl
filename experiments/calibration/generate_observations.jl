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
    ta_900hpa = OutputVar("/glade/u/home/nefrathe/clima/ClimaCoupler.jl/ta_900hpa.nc")
    ta_900hpa = reverse_dim(ta_900hpa, latitude_name(ta_900hpa))

    tas = OutputVar("/glade/u/home/nefrathe/clima/ClimaCoupler.jl/tas.nc")
    tas = reverse_dim(tas, latitude_name(tas))
    tas.attributes["short_name"] = "tas"

    ta_900hpa = slice(ta_900hpa; pressure = 900)

    tas_minus_ta_900hpa = tas - ta_900hpa
    tas_minus_ta_900hpa.attributes["short_name"] = "tas - ta"
    return [tas_minus_ta_900hpa, tas]
end

"""
    preprocess_vars(vars)

Preprocess each OutputVar in `vars` by keeping the relevant dates in
`sample_date_ranges`.
"""
function preprocess_vars(vars)
    out_var = OutputVar("/glade/derecho/scratch/nefrathe/tmp/output_quick_1/iteration_000/member_001/wxquest_diagedmf/output_active/clima_atmos/mslp_1week_average.nc")
    resample_var(x) = ClimaAnalysis.resampled_as(x, out_var ; dim_names = ["lon", "lat"])
    vars = resample_var.(vars)

    vars = map(vars) do var
        set_units(var, var_units[short_name(var)])
    end
    return vars
end

var_units = Dict(
    "pr" => "kg m^-2 s^-1",
    "mslp" => "Pa",
    "tas" => "K",
    "tas - ta" => "K",
)

resample_var() = resampled_lonlat()

function make_observation_vector(vars, sample_date_ranges)
    covar_estimator = ClimaCalibrate.ObservationRecipe.SVDplusDCovariance(
        sample_date_ranges;
        model_error_scale = 0.05,
        regularization = 1e-3,
        use_latitude_weights = true,
    )
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
    find_weekly_dates(sample_date_range)

Find weekly dates for the purpose of preprocessing the `OutputVar`s.
"""
function find_weekly_dates(vars, sample_date_ranges)
    # Find the earliest and latest years that can be used
    earliest_year = min((first(ClimaAnalysis.dates(var)) for var in vars)...)
    latest_year = max((last(ClimaAnalysis.dates(var)) for var in vars)...)

    # Check conversion is possible
    dates = Dates.DateTime[]
    for range in sample_date_ranges
        if ((last(range) - first(range)) % Dates.Week(1)).value != 0
            error(
                "The dates in $sample_date_ranges should differ by weeks",
            )
        end
        append!(dates, collect(first(range):Dates.Week(1):last(range)))
    end

    # Find all weekly dates
    min_year, max_year = extrema(year.(dates))

    alldates = Dates.DateTime[]
    for curr_date in dates
        for curr_year = year(earliest_year):year(latest_year)
            diff_year = year(curr_date) - curr_year
            push!(alldates, curr_date - Year(diff_year))
        end
    end
    return sort!(alldates)
end

"""
    find_weekly_ranges(sample_date_range)

Find the weekly ranges for the purpose of constructing the TSVD covariance
matrix.
"""
function find_weekly_ranges(vars, sample_date_range)
    var_dates = union([ClimaAnalysis.dates(var) for var in vars]...)
    min_year, max_year = year.(extrema(var_dates))
    curr_year = minimum(year.(extrema(var_dates)))

    date_ranges = typeof(sample_date_range)[]
    for curr_year = min_year:max_year
        delta_year = year(first(sample_date_range)) - curr_year
        date_range = sample_date_range .- Year(delta_year)
        if first(date_range) in var_dates && last(date_range) in var_dates
            push!(date_ranges, date_range)
        end
    end
    return sort!(date_ranges)
end


"""
    resampled_lonlat(config_file)

Return a function to resample longitude and latitudes according to the model
grid specified by `config_file`.
"""
function resampled_lonlat(config_file)
    cs = CoupledSimulation(config_file)
    center_space = cs.model_sims.atmos_sim.domain.center_space
    (lon_nlevels, lat_nlevels, z_nlevels) =
        ClimaDiagnostics.Writers.default_num_points(center_space)
    longitudes = range(-180, 180, lon_nlevels)
    latitudes = range(-90, 90, lat_nlevels)
    stretch = center_space.grid.vertical_grid.topology.mesh.stretch
    # TODO: Account for stretch for 3D variables and interpolate to pressure?
    dz_bottom = center_space.grid.vertical_grid.topology.mesh.faces[2].z
    z_levels = range(dz_bottom, ClimaCore.Spaces.z_max(center_space), z_nlevels)
    return var -> resampled_to(var; lon = longitudes, lat = latitudes)
end

if abspath(PROGRAM_FILE) == @__FILE__
    obsdir = "weekly"
    sample_date_ranges = CALIBRATE_CONFIG.sample_date_ranges
    short_names = CALIBRATE_CONFIG.short_names
    config_file = CALIBRATE_CONFIG.config_file
    @info "Generating observations for $short_names"
    @info "The number of samples is $(length(sample_date_ranges)) over $sample_date_ranges"

    # unprocessed_vars = load_vars(obsdir, short_names)
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
