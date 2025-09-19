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
include(joinpath(@__DIR__, "observation_utils.jl"))

function make_observation_vector(
    short_names,
    sample_date_ranges,
    obsdir,
    config_file,
)
    vars = load_vars(obsdir, short_names)
    vars = preprocess_vars(vars, sample_date_ranges, config_file)
    return make_observation_vector(vars, sample_date_ranges)
end

"""
    load_vars(obsdir, short_names)

Load NetCDF files belonging to `short_names` in `obsdir` as `OutputVar`s.
"""
function load_vars(obsdir, short_names)
    # Missing years at the end
    available_short_names = ("pr", "mslp", "tas") # pr and tas over land
    all(short_name -> short_name in available_short_names, short_names) ||
        error("Variable names must be in $(keys(varname_to_filepath))")
    allfiles = readdir(obsdir, join = true)
    vars = map(short_names) do short_name
        files = sort(filter(file -> occursin(short_name, file), allfiles))
        var = ClimaAnalysis.OutputVar(files, short_name)
        var.attributes["short_name"] = short_name
        var
    end
    return vars
end

"""
    preprocess_vars(vars, sample_date_ranges, config_file)

Preprocess each OutputVar in `vars` by keeping the relevant dates in
`sample_date_ranges`.


"""
function preprocess_vars(vars, sample_date_ranges, config_file)
    weekly_dates = find_weekly_dates(vars, sample_date_ranges)
    weekly_dates = unique(weekly_dates)
    vars = map(vars) do var
        var = ClimaAnalysis.select(var; by = ClimaAnalysis.MatchValue(), time = weekly_dates)
        # TODO: Do any additional preprocessing here for units...
        if !issorted(var.dims[ClimaAnalysis.latitude_name(var)])
            ClimaAnalysis.reverse_dim!(var, ClimaAnalysis.latitude_name(var))
        end
        @assert issorted(var.dims[ClimaAnalysis.latitude_name(var)])
        var = resampled_as(shift_longitude(var, -180.0, 180.0), diagnostic_var2d, dim_names = ["longitude", "latitude"])
        # apply ocean mask
        if ClimaAnalysis.short_name(var) in ("pr", "tas")
            var = ClimaAnalysis.apply_oceanmask(var)
        end

        if ClimaAnalysis.short_name(var) == "pr"
            # Change sign and convert
            # 1 mm / week= (1 kg m-2) / 604_800 s
            var = ClimaAnalysis.convert_units(var, "mm", conversion_function = x -> -x/604_800)
        end

        var.attributes["units"] = var_units[ClimaAnalysis.short_name(var)]
        var
    end
    return vars
end

var_units = Dict(
    "pr" => "kg m^-2 s^-1",
    "mslp" => "Pa",
    "tas" => "K",
)

resample_var() = resampled_lonlat()

function make_observation_vector(vars, sample_date_ranges)
    obs_vec = map(sample_date_ranges) do sample_date_range
        weekly_range = find_weekly_ranges(vars, sample_date_range)
        covar_estimator = ClimaCalibrate.ObservationRecipe.SVDplusDCovariance(
            weekly_range;
            model_error_scale = 0.05,
            regularization = 1e-5,
        )
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

    @info "The number of samples is $(length(sample_date_ranges))"

    diagnostic_var2d = OutputVar("/glade/derecho/scratch/nefrathe/tmp/output_quick_1/iteration_000/member_001/wxquest_diagedmf/output_active/clima_atmos/mslp_1week_average.nc")

    unprocessed_vars = load_vars(obsdir, short_names)
    preprocessed_vars = preprocess_vars(unprocessed_vars, sample_date_ranges, config_file)
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

#= TODO: 
- add ocean mask for tas, pr: ClimaAnalysis.apply_oceanmask(...)
- check sign and units of era5 precip data
- use 2+ weeks of data
- verify land radiation start date
- plot spread of loss
- exclude precip
- compute cov weighted bias per variable
- change plotted loss
- use is_complete function for gensemblebuilder, overwrite failed columns with nans
- sanity check: simplify observations to use spatial avg, e.g. mean surface temp in northern hemisphere
=#

# Different Septembers over the years for a single iteration
# Geopotential height or specific humidity or clouds
# Use average but for pr, convert average to sum
# precipation time scale
