#=
    Unit tests for ClimaCoupler Input module
=#
using Test
import ArgParse
import Dates
import ClimaCoupler: Input, Utilities
import ClimaCoupler
import YAML

@testset "test argparse_settings and parse_commandline" begin
    settings = Input.argparse_settings()

    # Test with empty ARGS (should use defaults)
    empty!(ARGS)

    parsed = Input.parse_commandline(settings)
    @test parsed["config_file"] ==
          joinpath(pkgdir(ClimaCoupler), "config/ci_configs/amip_default.yml")
    @test parsed["FLOAT_TYPE"] == "Float64"
    @test parsed["mode_name"] == "amip"  # default value

    # Test with custom arguments
    # We use try/finally to ensure ARGS is restored even if the test fails,
    # since ARGS is a global variable that could affect other tests
    try
        push!(ARGS, "--mode_name", "slabplanet")
        push!(ARGS, "--t_end", "1000secs")
        parsed = Input.parse_commandline(settings)
        @test parsed["mode_name"] == "slabplanet"
        @test parsed["t_end"] == "1000secs"
    finally
        empty!(ARGS)
    end
end

@testset "get_coupler_config_dict" begin
    # Use a simple test config file
    config_file = joinpath(pkgdir(ClimaCoupler), "test", "config", "input_test_config.yml")

    # Check that CLI arguments are parsed
    try
        push!(ARGS, "--dt_cpl", "600secs")
        config_dict = Input.get_coupler_config_dict(config_file)
        @test config_dict["dt_cpl"] == "120secs" # value in the config file
        @test config_dict["mode_name"] == "amip"
        @test config_dict["job_id"] == "input_test_config" # default to file name

        # Test that atmos config file is overwritten by coupler config file
        @test config_dict["atmos_config_file"] == "test/config/input_test_atmos_config.yml"
        @test config_dict["h_elem"] == 6 # 6 in coupler config, 16 in atmos config
    finally
        empty!(ARGS)
    end
end

@testset "get_coupler_args" begin
    # Create a minimal config dict for testing
    config_dict = Dict(
        "job_id" => "test_job",
        "mode_name" => "amip",
        "use_itime" => true,
        "unique_seed" => false,
        "FLOAT_TYPE" => "Float64",
        "t_end" => "800secs",
        "t_start" => "0secs",
        "start_date" => "20000101",
        "dt_cpl" => "400secs",
        "dt" => "400secs",
        "share_surface_space" => true,
        "checkpoint_dt" => "90days",
        "detect_restart_files" => false,
        "restart_dir" => nothing,
        "restart_t" => nothing,
        "restart_cache" => true,
        "save_cache" => true,
        "use_coupler_diagnostics" => true,
        "use_land_diagnostics" => true,
        "evolving_ocean" => true,
        "energy_check" => false,
        "conservation_softfail" => false,
        "rmse_check" => false,
        "coupler_output_dir" => "test_output",
        "land_model" => "bucket",
        "land_temperature_anomaly" => "aquaplanet",
        "land_spun_up_ic" => true,
        "bucket_albedo_type" => "map_static",
        "bucket_initial_condition" => "",
        "coupler_toml" => [],
        "era5_initial_condition_dir" => nothing,
        "ocean_model" => "prescribed",
        "ice_model" => "prescribed",
        "land_fraction_source" => "etopo",
        "binary_area_fraction" => true,
        "component_dt_dict" => Dict(
            "dt_atmos" => 400.0,
            "dt_land" => 400.0,
            "dt_ocean" => 400.0,
            "dt_seaice" => 400.0,
        ),
        "print_config_dict" => true,
    )

    args = Input.get_coupler_args(config_dict)
    # Test that expected fields are present
    @test haskey(args, :job_id)
    @test haskey(args, :sim_mode)
    @test haskey(args, :FT)
    @test haskey(args, :t_end)
    @test haskey(args, :t_start)
    @test haskey(args, :Δt_cpl)
    @test haskey(args, :component_dt_dict)

    # Test some values
    @test args.job_id == "test_job"
    @test args.FT == Float64
    @test args.sim_mode == ClimaCoupler.Interfacer.AMIPMode
    @test args.land_model == Val(:bucket)
    @test args.ocean_model == Val(:prescribed)
    @test args.ice_model == Val(:prescribed)
    @test args.land_fraction_source == "etopo"

    # Test that component_dt_dict is preserved
    @test args.component_dt_dict isa Dict
    @test haskey(args.component_dt_dict, "dt_atmos")
    @test !haskey(args.component_dt_dict, "dt")
end

@testset "get_diag_period" begin
    secs_per_day = 86400

    # Test for simulation longer than 90 days (monthly means)
    t_start = 0.0
    t_end = 100.0 * secs_per_day  # 100 days
    period, diagnostics_dt = Input.get_diag_period(t_start, t_end)
    @test period == "1months"
    @test diagnostics_dt == Dates.Month(1)

    # Test for simulation between 30 and 90 days (10-day means)
    t_end = 50.0 * secs_per_day  # 50 days
    period, diagnostics_dt = Input.get_diag_period(t_start, t_end)
    @test period == "10days"
    @test diagnostics_dt == Dates.Day(10)

    # Test for simulation between 1 and 30 days (daily means)
    t_end = 10.0 * secs_per_day  # 10 days
    period, diagnostics_dt = Input.get_diag_period(t_start, t_end)
    @test period == "1days"
    @test diagnostics_dt == Dates.Day(1)

    # Test for simulation shorter than 1 day (hourly means)
    t_end = 0.5 * secs_per_day  # 12 hours
    period, diagnostics_dt = Input.get_diag_period(t_start, t_end)
    @test period == "1hours"
    @test diagnostics_dt == Dates.Hour(1)

    # Test boundary cases
    # Exactly 90 days should be monthly
    t_end = 90.0 * secs_per_day
    period, diagnostics_dt = Input.get_diag_period(t_start, t_end)
    @test period == "1months"

    # Exactly 30 days should be 10-day means
    t_end = 30.0 * secs_per_day
    period, diagnostics_dt = Input.get_diag_period(t_start, t_end)
    @test period == "10days"

    # Exactly 1 day should be daily
    t_end = 1.0 * secs_per_day
    period, diagnostics_dt = Input.get_diag_period(t_start, t_end)
    @test period == "1days"
end

@testset "parse_component_dts!" begin
    # Test case 1: All component dt's are specified
    config_dict = Dict{String, Any}(
        "dt_cpl" => "400secs",
        "dt_atmos" => "200secs",
        "dt_land" => "200secs",
        "dt_ocean" => "400secs",
        "dt_seaice" => "400secs",
        "dt" => "300secs",  # Should be removed
    )
    @test_logs (:warn, "Removing dt in favor of individual component dt's") Input.parse_component_dts!(
        config_dict,
    )

    @test haskey(config_dict, "component_dt_dict")
    @test !haskey(config_dict, "dt")  # Should be removed
    @test config_dict["component_dt_dict"]["dt_atmos"] == 200.0
    @test config_dict["component_dt_dict"]["dt_land"] == 200.0
    @test config_dict["component_dt_dict"]["dt_ocean"] == 400.0
    @test config_dict["component_dt_dict"]["dt_seaice"] == 400.0

    # Test case 2: Only generic dt is specified
    config_dict = Dict{String, Any}("dt_cpl" => "400secs", "dt" => "200secs")
    Input.parse_component_dts!(config_dict)

    @test haskey(config_dict, "component_dt_dict")
    @test haskey(config_dict, "dt")  # Should remain
    @test config_dict["component_dt_dict"]["dt_atmos"] == 200.0
    @test config_dict["component_dt_dict"]["dt_land"] == 200.0
    @test config_dict["component_dt_dict"]["dt_ocean"] == 200.0
    @test config_dict["component_dt_dict"]["dt_seaice"] == 200.0

    # Test case 3: Some (but not all) component dt's specified - should use generic dt
    config_dict = Dict{String, Any}(
        "dt_cpl" => "400secs",
        "dt" => "200secs",
        "dt_atmos" => "100secs",  # Should be removed with warning
        "dt_land" => nothing,
        "dt_ocean" => nothing,
        "dt_seaice" => nothing,
    )
    @test_logs (
        :warn,
        "Removing dt_atmos from config in favor of dt because not all component dt's are specified",
    ) Input.parse_component_dts!(config_dict)

    @test haskey(config_dict, "component_dt_dict")
    @test !haskey(config_dict, "dt_atmos")  # Should be removed
    @test config_dict["component_dt_dict"]["dt_atmos"] == 200.0
    @test config_dict["component_dt_dict"]["dt_land"] == 200.0
    @test config_dict["component_dt_dict"]["dt_ocean"] == 200.0
    @test config_dict["component_dt_dict"]["dt_seaice"] == 200.0

    # Test case 4: Error when dt is missing and not all component dt's are specified
    config_dict = Dict{String, Any}(
        "dt_cpl" => "400secs",
        "dt_atmos" => "200secs",
        "dt_land" => nothing,
        "dt_ocean" => nothing,
        "dt_seaice" => nothing,
    )
    @test_throws "dt or (dt_atmos, dt_land, dt_ocean, and dt_seaice) must be specified" Input.parse_component_dts!(
        config_dict,
    )

    # Test case 5: Error when component dt is not divisible by coupler dt
    config_dict = Dict{String, Any}(
        "dt_cpl" => "400secs",
        "dt_atmos" => "300secs",  # 400 is not divisible by 300
        "dt_land" => "200secs",
        "dt_ocean" => "200secs",
        "dt_seaice" => "200secs",
    )
    @test_throws "Coupler dt must be divisible by all component dt's" Input.parse_component_dts!(
        config_dict,
    )

    # Test case 6: Warning when dt is removed in favor of component dt's
    config_dict = Dict{String, Any}(
        "dt_cpl" => "400secs",
        "dt" => "200secs",
        "dt_atmos" => "200secs",
        "dt_land" => "200secs",
        "dt_ocean" => "200secs",
        "dt_seaice" => "200secs",
    )
    @test_logs (:warn, "Removing dt in favor of individual component dt's") Input.parse_component_dts!(
        config_dict,
    )
    @test !haskey(config_dict, "dt")
end
