using Test
import ClimaCoupler: SimOutput
import ArgParse

@testset "get_benchmark_args" begin
    # Test with empty ARGS (should use defaults)
    # We need to temporarily set ARGS
    old_args = ARGS
    try
        empty!(ARGS)
        push!(ARGS, "--job_id_coupled", "test_coupled")
        args = SimOutput.get_benchmark_args()
        @test args["job_id_coupled"] == "test_coupled"
        @test isnothing(args["job_id_atmos"])
        @test args["coupler_output_dir"] == "experiments/ClimaEarth/output"
    finally
        empty!(ARGS)
        append!(ARGS, old_args)
    end
end

@testset "get_run_info" begin
    test_output_dir = mktempdir()
    parsed_args = Dict(
        "job_id_coupled" => "test_coupled",
        "job_id_atmos" => "test_atmos",
        "job_id_coupled_io" => "test_coupled_io",
        "job_id_atmos_diagedmf" => "test_atmos_diagedmf",
        "job_id_coupled_progedmf_coarse" => "test_coarse",
        "job_id_coupled_progedmf_fine" => "test_fine",
        "coupler_output_dir" => test_output_dir,
    )

    # Test each run type
    for (run_type, expected_job_id) in [
        ("coupled", "test_coupled"),
        ("atmos", "test_atmos"),
        ("coupled_io", "test_coupled_io"),
        ("atmos_diagedmf", "test_atmos_diagedmf"),
        ("coupled_progedmf_coarse", "test_coarse"),
        ("coupled_progedmf_fine", "test_fine"),
    ]
        job_id, artifacts_dir = SimOutput.get_run_info(parsed_args, run_type)
        @test job_id == expected_job_id
        @test artifacts_dir == joinpath(test_output_dir, expected_job_id, "artifacts")
    end

    # Test invalid run type
    @test_throws ErrorException SimOutput.get_run_info(parsed_args, "invalid_type")

    # Test missing job_id
    parsed_args_no_job =
        Dict("job_id_coupled" => nothing, "coupler_output_dir" => test_output_dir)
    @test_throws ErrorException SimOutput.get_run_info(parsed_args_no_job, "coupled")
end

@testset "append_table_data" begin
    # Create a temporary directory with sypd.txt
    test_dir = mktempdir()
    artifacts_dir = joinpath(test_dir, "artifacts")
    mkpath(artifacts_dir)

    # Write test SYPD value
    sypd_value = 123.456789
    open(joinpath(artifacts_dir, "sypd.txt"), "w") do f
        println(f, sypd_value)
    end

    # Test appending data
    initial_data = [
        ["Column1" "Column2" "Column3"]
        ["" "" ""]
    ]
    setup_id = "Test Setup"
    job_id = "test_job"

    result = SimOutput.append_table_data(initial_data, setup_id, job_id, artifacts_dir)

    @test size(result) == (4, 3) # 4 rows, 3 columns
    @test result[1, :] == ["Column1", "Column2", "Column3"]
    @test result[3, :] == [setup_id, "job ID:", job_id]
    @test result[:, 1] == ["Column1", "", "Test Setup", ""]
    @test result[4, 2] == "SYPD:"
    @test result[4, 3] == round(sypd_value, digits = 4)
end
