import ClimaAnalysis
import ClimaAnalysis: OutputVar
import ClimaCalibrate
import Dates
import EnsembleKalmanProcesses as EKP
import JLD2

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
    vars = map(vars) do var
        # resample_var(select_dates(var, weekly_dates))
        var = select_dates(var, weekly_dates)
        # TODO: Do any additional preprocessing here for units...
    end
    return vars
end

function make_observation_vector(vars, sample_date_ranges)
    obs_vec = map(sample_date_ranges) do sample_date_range
        weekly_range = find_weekly_ranges(vars, sample_date_range)
        covar_estimator = ClimaCalibrate.ObservationRecipe.SVDplusDCovariance(
            weekly_range;
            model_error_scale = 0.05,
            regularization = 1e-2,
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
    select_dates(var, dates)

Select `dates` in `var`.
"""
function select_dates(var::OutputVar, dates)
    # TODO: This function can be simplfied with select
    # Metadata comes from observational data and var is simulation data
    # when calibrating
    sim_dates = ClimaAnalysis.dates(var)
    common_date_indices = indexin(dates, sim_dates)

    any(isnothing, common_date_indices) && error(
        "There are dates ($(setdiff(dates, sim_dates))) that are not present in the dates in var",
    )

    var_time_name = ClimaAnalysis.time_name(var)
    dims = deepcopy(var.dims)
    dims[var_time_name] = dims[var_time_name][common_date_indices]

    time_idx = var.dim2index[var_time_name]
    time_indices =
        ntuple(x -> ifelse(x == time_idx, common_date_indices, Colon()), ndims(var.data))
    data = view(var.data, time_indices...)
    return ClimaAnalysis.remake(var, dims = dims, data = data)
end

if abspath(PROGRAM_FILE) == @__FILE__
    obsdir = joinpath(@__DIR__, "../../weekly")
    sample_date_ranges = CALIBRATE_CONFIG.sample_date_ranges
    short_names = CALIBRATE_CONFIG.short_names
    config_file = CALIBRATE_CONFIG.config_file

    @info "The number of samples is $(length(sample_date_ranges))"

    vars = load_vars(obsdir, short_names)
    vars = preprocess_vars(vars, sample_date_ranges, config_file)
    observation_vector = make_observation_vector(vars, sample_date_ranges)
    JLD2.save_object(
        joinpath(@__DIR__, "weatherquest_obs_vec.jld2"),
        observation_vector,
    )
end


# Different Septembers over the years for a single iteration
# Geopotential height or specific humidity or clouds
# Use average but for pr, convert average to sum
# precipation time scale
