redirect_stderr(IOContext(stderr, :stacktrace_types_limited => Ref(false)))

using ClimaComms
comms_ctx = ClimaComms.context()
const pid, nprocs = ClimaComms.init(comms_ctx)


import SciMLBase: ODEProblem, solve, step!, init, reinit!
using LinearAlgebra
import Test: @test
using Dates
using Plots
using Statistics: mean
import ClimaAtmos as CA
import YAML

using ClimaCore.Utilities: half, PlusHalf
using ClimaCore: InputOutput, Fields
import ClimaCore.Spaces as Spaces

## coupler specific imports
import ClimaCoupler
import ClimaCoupler.Regridder
import ClimaCoupler.Regridder:
    update_surface_fractions!, combine_surfaces!, combine_surfaces_from_sol!, dummmy_remap!, binary_mask
import ClimaCoupler.ConservationChecker:
    EnergyConservationCheck, WaterConservationCheck, check_conservation!, plot_global_conservation
import ClimaCoupler.Utilities: swap_space!
import ClimaCoupler.BCReader:
    bcfile_info_init, float_type_bcf, update_midmonth_data!, next_date_in_file, interpolate_midmonth_to_daily
import ClimaCoupler.TimeManager:
    current_date,
    datetime_to_strdate,
    trigger_callback,
    Monthly,
    EveryTimestep,
    HourlyCallback,
    MonthlyCallback,
    update_firstdayofmonth!,
    trigger_callback!
import ClimaCoupler.Diagnostics: get_var, init_diagnostics, accumulate_diagnostics!, save_diagnostics, TimeMean
import ClimaCoupler.PostProcessor: postprocess

import ClimaCoupler.Interfacer:
    CoupledSimulation,
    float_type,
    AtmosModelSimulation,
    SurfaceModelSimulation,
    SurfaceStub,
    SeaIceModelSimulation,
    LandModelSimulation,
    OceanModelSimulation,
    get_field,
    update_field!
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

## helpers for component models
include("../experiments/AMIP/components/atmosphere/climaatmos_init.jl")
include("../experiments/AMIP/components/land/bucket_init.jl")
include("../experiments/AMIP/components/land/bucket_utils.jl")
include("../experiments/AMIP/components/ocean/slab_ocean_init.jl")
include("../experiments/AMIP/components/ocean/prescr_seaice_init.jl")
include("../experiments/AMIP/user_io/user_diagnostics.jl")
include("../experiments/AMIP/user_io/user_logging.jl")
