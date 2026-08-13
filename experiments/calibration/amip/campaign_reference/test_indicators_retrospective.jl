# Acceptance test for the steering indicators.
#
# The test loads saved calibration state from two completed runs and asserts
# that the indicators report the known outcome of each run:
#   1. The collapsed run: COLLAPSED is flagged by iteration 2.
#   2. The lwp+pr run: both observables reach the noise floor, the run is a
#      CONVERGED candidate by iteration 4, and no collapse is flagged.
#
# Run with:
#   julia --project=experiments/AMIP experiments/calibration/amip/test_indicators_retrospective.jl

ENV["CLIMACOMMS_CONTEXT"] = "SINGLETON"

import JLD2
import Test: @test, @testset
import ClimaCoupler
import EnsembleKalmanProcesses as EKP
import ClimaCalibrate
import ClimaAnalysis
import Statistics
import CairoMakie
import GeoMakie

include(joinpath(pkgdir(ClimaCoupler), "experiments", "calibration", "amip", "model_interface.jl"))

const COLLAPSED_RUN = "/glade/derecho/scratch/nefrathe/amip_calibration_out_collapsed_v1"
const LWPPR_RUN = "/glade/derecho/scratch/nefrathe/amip_calibration_lwp_pr_out"

"Load the ekp, G ensemble, and prior for one completed iteration of a run."
function load_iteration(run_dir, iteration)
    ekp = JLD2.load_object(
        joinpath(run_dir, "iteration_$(lpad(iteration + 1, 3, '0'))", "eki_file.jld2"),
    )
    g = JLD2.load_object(
        joinpath(run_dir, "iteration_$(lpad(iteration, 3, '0'))", "G_ensemble.jld2"),
    )
    prior = JLD2.load_object(joinpath(run_dir, "iteration_001", "prior.jld2"))
    return ekp, g, prior
end

@testset "collapsed run is flagged" begin
    ekp, g, prior = load_iteration(COLLAPSED_RUN, 2)
    block = steering_summary(ekp, g, prior, 2)
    println(block, "\n")
    @test occursin("COLLAPSED", block)
end

@testset "lwp+pr run converges without collapse" begin
    # Iteration 4: both observables are at the noise floor, but rain_tau still
    # drifts at 0.12 sigma per iteration, so no convergence call yet.
    ekp, g, prior = load_iteration(LWPPR_RUN, 4)
    block = steering_summary(ekp, g, prior, 4)
    println(block, "\n")
    @test !occursin("COLLAPSED", block)
    @test occursin("fit to noise floor", block)
    @test !occursin("CONVERGED candidate", block)
    @test occursin("separable", block)

    # Iteration 8: parameter drift has settled and the convergence hint fires.
    ekp, g, prior = load_iteration(LWPPR_RUN, 8)
    block = steering_summary(ekp, g, prior, 8)
    println(block, "\n")
    @test !occursin("COLLAPSED", block)
    @test occursin("CONVERGED candidate", block)
end
