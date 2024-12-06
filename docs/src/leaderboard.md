# Leaderboard

## AMIP Driver

### Add a new variable to compare against observations
Computing errors against observations are all contained in the `leaderboard` folder. The
files in the leaderboard folder are `data_sources.jl` and `leaderboard.jl`. Loading and
preprocessing variables of interest are done in `data_sources.jl` and computing the errors
and plotting are done in `leaderboard.jl`. To add a new variable, you ideally only need to
modify `data_sources.jl`.

#### Add a new 3D variable to the bias plots
If you want to add a new 3D variable defined over latitude, longitude, and time to the bias
plots, you add the variable to `sim_var_dict`, `obs_var_dict`, `compare_vars_biases_groups`,
and optionally `compare_vars_biases_plot_extrema`. The variables `sim_var_dict`,
`obs_var_dict`, `compare_vars_biases_groups`, `compare_vars_biases_plot_extrema` are in the
function `get_sim_var_dict`, `get_obs_var_dict`, `get_compare_vars_biases_groups`, and
`get_compare_vars_biases_plot_extrema` respectively.

The dictionaries `sim_var_dict` and `obs_var_dict` map short names to an anonymous function
that returns a [`OutputVar`](https://clima.github.io/ClimaAnalysis.jl/dev/var/). Both
dictionaries must use the same short names as the keys so that the right simulation and
observational data are compared.

!!! tip "Preprocessing"
    Observational and simulational data should be preprocessed for dates and units. For
    simulation data, monthly averages correspond to the first day following the month.
    For instance, the monthly average corresponding to January 2010 is on the date
    2/1/2010. Preprocessing is done to shift this date to 1/1/2010. When preprocessing
    data, we follow the convention that the first day corresponds to the monthly average
    for that month. For observational data, you should check the convention being followed
    and preprocess the dates if necessary.

    For `obs_var_dict`, the anonymous function must take in a start date. The start date is
    used in `leaderboard.jl` to adjust the seconds in the `OutputVar` to match between start
    date in the simulation data.

    Units should be the same between the simulation and observational data.

The variable `compare_vars_biases_groups` is an array of arrays of short names that control
which variables are plotted together. You can add the variable to an existing array or make
a new array. The dictionary `compare_vars_biases_plot_extrema` maps short names to tuples.
The dictionary sets the lower and upper bounds of the bias plots.

#### Add a new variable to the leaderboard
If you want to add a new variable to the leaderboard, you add the variable to
`rmse_var_names` and `rmse_var_dict`. The array `rmse_var_names` is a list of short names.
The dictionary `rmse_var_dict` maps short name to
[`RMSEVariable`](https://clima.github.io/ClimaAnalysis.jl/dev/rmse_var/). A `RMSEVariable`
must be initialized for each variable of interest. The CliMA model is added with units to
the `RMSEVariable`. It is assumed that the `RMSEVariable` contains only the columns "DJF",
"MAM", "JJA", "SON", and "ANN" in that order. The file `leaderboard.jl` will load the
appropriate data into the `RMSEVariable`.

### Add a new variable to compare against observations in pressure coordinates
To add a new variable, you only need to modify the variable `sim_var_pfull_dict` in the
function `get_sim_var_in_pfull_dict`, the variable `obs_var_dict` in the function
`get_obs_var_in_pfull_dict`, and the variable `compare_vars_biases_plot_extrema` in the
function `get_compare_vars_biases_plot_extrema_pfull`. The variables and functions are
defined exactly the same as their analogous versions in the section above.

It is expected that the dimensions of the variables are time, latitude, longitude, and
pressure in no particular order and the units for the pressure dimension is expected to be
`hPa`.
