import Dates

"""
    struct CalibrateConfig{SPINUP <: Dates.Period, EXTEND <: Dates.Period}
        short_names::Vector{String}
        minibatch_size::Int64
        n_iterations::Int64
        sample_date_ranges::Vector{NTuple{2, DATE}}
        extend::EXTEND
        spinup::SPINUP
        nelements::Tuple{Int64, Int64}
        output_dir::String
        rng_seed::Int64
    end

A configuration struct for keeping track of multiple fields that are of interest
to a user running calibration, or that are needed in multiple places (e.g., for
ensemble members and generating observations).
"""
struct CalibrateConfig{SPINUP <: Dates.Period, EXTEND <: Dates.Period}
    "Configuration file to use for ClimaCoupler simulation"
    config_file::String

    "The short names of the observations used for calibration. The short names
    should match the same names used for the diagnostics."
    short_names::Vector{String}

    "The size of the minibatch for each iteration"
    minibatch_size::Int64

    "The number of iterations to run the calibration for"
    n_iterations::Int64

    "The date ranges of the samples for calibration and used to determine the
    start and end dates of a simulation for each iteration of calibration"
    sample_date_ranges::Vector{NTuple{2, Dates.DateTime}}

    "The amount of time to run a simulation after the last date of the
    minibatch"
    extend::EXTEND

    "The amount of time to run a simulation before the first date of the
    minibatch"
    spinup::SPINUP

    "The directory to store the iterations and members of the calibration."
    output_dir::String

    "An integer value for ensuring calibrations are the same between multiple
    calibrations with the same settings"
    rng_seed::Int64
end

"""
    CalibrateConfig(;
        config_file,
        short_names,
        sample_date_ranges,
        extend,
        spinup = Dates.Month(3),
        minibatch_size,
        n_iterations,
        output_dir = "calibration/weatherquest",,
        rng_seed = 42,
    )

Initializes a CalibrateConfig, which is of interest to a user running
calibration or contains values needed in multiple places during calibration.

Keyword arguments
=====================

- `config_file`: Configuration file to use for ClimaCoupler simulation.

- `short_names`: Short names of the observations. The currently supported short
  names are `pr`, `tas`, and `mslp`.

- `minibatch_size`: The size of the minibatch for each iteration.

- `n_iterations`: The number of iterations to run the calibration for.

- `sample_date_ranges`: The date ranges for each sample. The dates should be the
  same as found in the time series data of the observations.

- `extend`: The amount of time to run the simulation after the end date
  determined by `sample_date_ranges`. For seasonal averages, `extend` should be
  `Dates.Month(3)` and for monthly averages, `extend` should be
  `Dates.Month(1)`.

- `spinup`: The amount of time to run the simulation before the start date
  determined by `sample_date_ranges`.

- `nelements`: The resolution of the model. This is also used to determine the
  mask of the observations.

- `output_dir`: The location to save the calibration at.

- `rng_seed`: An integer to ensure that calibration runs with the same settings
  are the same.
"""
function CalibrateConfig(;
    config_file,
    short_names,
    minibatch_size,
    n_iterations,
    sample_date_ranges,
    extend,
    spinup = Dates.Month(3),
    output_dir = "calibration/weatherquest",
    rng_seed = 42,
)
    isempty(short_names) && error("Cannot run calibration with no short names")
    isempty(sample_date_ranges) &&
        error("Cannot run calibration with no date ranges for the samples")

    sample_date_ranges = [
        (Dates.DateTime(date_pair[1]), Dates.DateTime(date_pair[2])) for
        date_pair in sample_date_ranges
    ]

    for (start_date, stop_date) in sample_date_ranges
        start_date <= stop_date || error(
            "The start date ($start_date) should be before the stop date ($stop_date)",
        )
    end
    issorted(sample_date_ranges) ||
        error("The samples in $sample_date_ranges should be sorted")

    minibatch_size > 0 || error("The minibatch size ($minibatch_size) should be positive")
    n_iterations > 0 || error("The number of iterations ($n_iterations) should be positive")

    num_samples = length(sample_date_ranges)
    minibatch_size > num_samples && error(
        "The minibatch size is $minibatch_size, but the number of samples is $num_samples",
    )

    remaining = num_samples % minibatch_size
    remaining == 0 || @warn(
        "Number of samples is not divisible by the minibatch size; the last $remaining samples may be missing when running the calibration"
    )

    return CalibrateConfig(
        config_file,
        short_names,
        minibatch_size,
        n_iterations,
        sample_date_ranges,
        extend,
        spinup,
        output_dir,
        rng_seed,
    )

end
