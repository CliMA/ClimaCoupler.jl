"""
    ClimaCoupler

Module for atmos-ocean-land coupled simulations.
"""
module ClimaCoupler

include("Interfacer.jl")
include("Utilities.jl")
include("TimeManager.jl")
include("ConservationChecker.jl")
include("FluxCalculator.jl")
include("FieldExchanger.jl")
include("Checkpointer.jl")
include("Input.jl")
include("SimOutput/SimOutput.jl")
include("Plotting.jl")
include("SimCoordinator.jl")
include("Models.jl")
include("CalibrationTools.jl")

# Import key functions from submodules to re-export at top level
import ..Interfacer: CoupledSimulation
import ..SimCoordinator: run!, step!, setup_and_run
import ..Plotting: postprocess

# Export all modules and key functions
export CalibrationTools,
    ConservationChecker,
    Checkpointer,
    FieldExchanger,
    FluxCalculator,
    Input,
    Interfacer,
    Models,
    Plotting,
    SimCoordinator,
    SimOutput,
    TimeManager,
    Utilities,
    CoupledSimulation,
    run!,
    step!,
    setup_and_run,
    postprocess

# ---------------------------------------------------------------------------
# Precompile workload
#
# Cache host-side type inference and native code for coupler-owned construction
# paths, so a cold process skips re-inferring them on first use. `ClimaCoupler`
# does not depend on `ClimaAtmos` (it is a weak dependency behind an extension),
# and `PrecompileTools` can only cache specializations owned by this package or
# its dependencies. This workload therefore exercises only ClimaAtmos-free paths:
# it builds a column boundary space (pure ClimaCore / ClimaComms) and allocates
# the default coupler fields on it. The heavy `get_simulation` / `step!`
# compilation lives in ClimaAtmos and must be precompiled there.
#
# Keep this block fast and free of artifacts, network, and GPU: it runs on every
# `Pkg.precompile()`.
# ---------------------------------------------------------------------------
import PrecompileTools
import ClimaComms

PrecompileTools.@setup_workload begin
    comms_ctx = ClimaComms.SingletonCommsContext(ClimaComms.CPUSingleThreaded())
    PrecompileTools.@compile_workload begin
        for FT in (Float32, Float64)
            boundary_space = Utilities.create_boundary_space(
                FT,
                "column",
                nothing,
                comms_ctx;
                column_latlon = (FT(0), FT(0)),
            )
            field_names = Interfacer.default_coupler_fields()
            Interfacer.init_coupler_fields(FT, field_names, boundary_space)
        end
    end
end

end
