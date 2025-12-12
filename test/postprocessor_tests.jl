#=
    Unit tests for ClimaCoupler Postprocessor module
=#
using Test
import ClimaComms
ClimaComms.@import_required_backends
import ClimaCoupler: Postprocessor

# Initialize MPI context
ClimaComms.init(ClimaComms.context())

@testset "Postprocessor tests" begin
    @testset "simulated_years_per_day" begin
        # Test case 1: 0.5 year simulated in 1 day of walltime
        cs1 = (; tspan = [0.0, 0.5 * 365.25 * 86400.0])  # 0.5 year in seconds
        walltime1 = 86400.0  # 1 day in seconds
        sypd1 = Postprocessor.simulated_years_per_day(cs1, walltime1)
        @test sypd1 ≈ 0.5

        # Test case 2: 1 year simulated in 2 days of walltime
        cs2 = (; tspan = [0.0, 365.25 * 86400.0])
        walltime2 = 2 * 86400.0  # 2 days in seconds
        sypd2 = Postprocessor.simulated_years_per_day(cs2, walltime2)
        @test sypd2 ≈ 0.5

        # Test case 3: 2 years simulated in 1 day of walltime
        cs3 = (; tspan = [0.0, 2 * 365.25 * 86400.0])
        walltime3 = 86400.0
        sypd3 = Postprocessor.simulated_years_per_day(cs3, walltime3)
        @test sypd3 ≈ 2.0

        # Test case 4: Non-zero start time
        cs4 = (; tspan = [1000.0, 1000.0 + 365.25 * 86400.0])
        walltime4 = 86400.0
        sypd4 = Postprocessor.simulated_years_per_day(cs4, walltime4)
        @test sypd4 ≈ 1.0
    end

    @testset "walltime_per_coupling_step" begin
        # Test case 1: 1 year with 1 day coupling steps, 1 day walltime
        cs1 = (; tspan = [0.0, 365.25 * 86400.0], Δt_cpl = 86400.0)  # 1 day steps
        walltime1 = 86400.0  # 1 day walltime
        n_steps = 365  # approximately 365 steps
        expected_walltime_per_step = walltime1 / n_steps
        actual_walltime_per_step = Postprocessor.walltime_per_coupling_step(cs1, walltime1)
        @test actual_walltime_per_step ≈ expected_walltime_per_step

        # Test case 2: 1 year with 1 hour coupling steps, 1 day walltime
        cs2 = (; tspan = [0.0, 365.25 * 86400.0], Δt_cpl = 3600.0)  # 1 hour steps
        walltime2 = 86400.0
        n_steps2 = 365 * 24  # approximately 365*24 steps
        expected_walltime_per_step2 = walltime2 / n_steps2
        actual_walltime_per_step2 = Postprocessor.walltime_per_coupling_step(cs2, walltime2)
        @test actual_walltime_per_step2 ≈ expected_walltime_per_step2

        # Test case 3: Exact calculation - 10 steps, 100 seconds walltime
        cs3 = (; tspan = [0.0, 100.0], Δt_cpl = 10.0)  # 10 steps
        walltime3 = 50.0  # 50 seconds walltime
        expected_walltime_per_step3 = 50.0 / 10.0  # 5 seconds per step
        actual_walltime_per_step3 = Postprocessor.walltime_per_coupling_step(cs3, walltime3)
        @test actual_walltime_per_step3 ≈ expected_walltime_per_step3

        # Test case 4: Non-zero start time
        cs4 = (; tspan = [1000.0, 1000.0 + 100.0], Δt_cpl = 10.0)
        walltime4 = 50.0
        expected_walltime_per_step4 = 50.0 / 10.0
        actual_walltime_per_step4 = Postprocessor.walltime_per_coupling_step(cs4, walltime4)
        @test actual_walltime_per_step4 ≈ expected_walltime_per_step4
    end

    @testset "save_sypd_walltime_to_disk" begin
        # Create a temporary directory for artifacts
        import Base.Filesystem: mktempdir, rm
        artifacts_dir = mktempdir()
        dir_paths = (; artifacts_dir = artifacts_dir)

        # Create a NamedTuple that mimics the minimal interface needed
        comms_ctx = ClimaComms.context()
        cs = (; tspan = [0.0, 365.25 * 86400.0], Δt_cpl = 86400.0, dir_paths = dir_paths)

        # Extend ClimaComms.context for NamedTuples in this test scope
        ClimaComms.context(cs::NamedTuple) = comms_ctx

        walltime = 86400.0

        # Only test if we're on root process
        if ClimaComms.iamroot(comms_ctx)
            Postprocessor.save_sypd_walltime_to_disk(cs, walltime)

            # Check that files were created
            sypd_file = joinpath(artifacts_dir, "sypd.txt")
            walltime_file = joinpath(artifacts_dir, "walltime_per_step.txt")

            @test isfile(sypd_file)
            @test isfile(walltime_file)

            # Read and verify contents
            sypd_value = parse(Float64, readchomp(sypd_file))
            expected_sypd = Postprocessor.simulated_years_per_day(cs, walltime)
            @test sypd_value ≈ expected_sypd

            walltime_value = parse(Float64, readchomp(walltime_file))
            expected_walltime = Postprocessor.walltime_per_coupling_step(cs, walltime)
            @test walltime_value ≈ expected_walltime
        end

        # Clean up: remove the temporary directory and its contents
        # This ensures cleanup happens regardless of whether we're on root process
        rm(artifacts_dir, recursive = true)
    end
end
