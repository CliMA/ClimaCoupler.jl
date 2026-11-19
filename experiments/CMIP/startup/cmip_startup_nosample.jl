# Phase-resolved breakdown of CMIP launch -> first coupled step.
#
#   CLIMACOMMS_DEVICE=CUDA julia -O0 --project=experiments/CMIP \
#       <this file> --config_file config/ci_configs/cmip_oceananigans_climaseaice.yml

import Profile
Profile.init(n = 4 * 10^7, delay = 0.01)
Base.cumulative_compile_timing(true)

include(joinpath(@__DIR__, "phase_probe.jl"))
using .PhaseProbe: @phase, @phase_nosample, report
PhaseProbe.DUMPDIR[] = get(ENV, "COUPLER_PERF_OUT", "")

# ----------------------------------------------------------------- load phases
# `using` is only legal at top level, so each import goes through `@eval Main`.
# These are not sampled: the profiler interacts badly with package loading.
# The warm half of an A/B. `CMIPWarmup` sits at the top of the layered stack in
# `experiments/warmup`, so loading it pulls in every warmup image; run with
# `--project=experiments/warmup`. Loading them is itself a cost worth attributing, which is
# why it is its own phase rather than folded into the imports below.
if get(ENV, "WARMUP_LOAD", "") == "1"
    @phase_nosample "load: warmup images" @eval Main using CMIPWarmup
end

@phase_nosample "load: ClimaCoupler" @eval Main using ClimaCoupler
@phase_nosample "load: Makie stack" @eval Main using CairoMakie,
    GeoMakie, Makie, Poppler_jll, Printf
@phase_nosample "load: Oceananigans+SeaIce" @eval Main import Oceananigans,
    ClimaSeaIce, KernelAbstractions, ConservativeRegridding, Adapt
@phase_nosample "load: ClimaLand" @eval Main import ClimaLand
@phase_nosample "load: ClimaAtmos" @eval Main import ClimaAtmos
@phase_nosample "load: ClimaComms backends" begin
    @eval Main import ClimaComms
    @eval Main ClimaComms.@import_required_backends
end

# ----------------------------------------------------------------- setup phase
config_file = Main.Input.parse_commandline(Main.Input.argparse_settings())["config_file"]
println("config_file = ", config_file)

cs = @phase_nosample "setup: CoupledSimulation" Main.CoupledSimulation(config_file)

# -------------------------------------------------------------- stepping phases
@phase_nosample "step 1" Main.step!(cs)

report()

# How many distinct GPU kernels were compiled, and how many host methods?
try
    CUDA = getfield(Main, :CUDA)
    n = length(CUDA.compiler_cache(CUDA.context()))
    println("GPU kernels compiled: ", n)
catch err
    println("GPU kernel count unavailable: ", err)
end
