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
include("user_io/user_diagnostics.jl")


COUPLER_OUTPUT_DIR = "experiments/AMIP/modular/output/amip/coarse_single_modular_businger_ft64"
COUPLER_ARTIFACTS_DIR = COUPLER_OUTPUT_DIR*"_artifacts"

@info COUPLER_OUTPUT_DIR
@info COUPLER_ARTIFACTS_DIR

include("user_io/plot_helper.jl")


# ## ClimaESM
# @info "AMIP plots"
# include("user_io/amip_visualizer.jl")
# post_spec = (;
# T = (:regrid, :zonal_mean),
# u = (:regrid, :zonal_mean),
# q_tot = (:regrid, :zonal_mean),
# toa  = (:regrid, :horizontal_slice),
# precipitation  = (:regrid, :horizontal_slice),
# T_sfc = (:regrid, :horizontal_slice),
# )

# plot_spec = (;
# T = (; clims = (190, 320), units = "K"),
# u = (; clims = (-50, 50), units = "m/s"),
# q_tot = (; clims = (0, 30), units = "g/kg"),
# toa = (; clims = (-250, 250), units = "W/m^2"),
# precipitation = (clims = (0, 1e-4), units = "kg/m^2/s"),
# T_sfc = (clims = (225, 310), units = "K"),
# )
# amip_data = amip_paperplots(
# post_spec,
# plot_spec,
# COUPLER_OUTPUT_DIR,
# files_root = ".monthly",
# output_dir = COUPLER_ARTIFACTS_DIR,
# )

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
ncep_data = ncep_paperplots(
    ncep_post_spec,
    ncep_plot_spec,
    COUPLER_OUTPUT_DIR,
    output_dir = COUPLER_ARTIFACTS_DIR,
    month_date = Dates.DateTime(1979, 01, 01),
) ## plot data that correspond to the model's last save_hdf5 call (i.e., last month)

