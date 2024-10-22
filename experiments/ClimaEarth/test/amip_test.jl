## AMIP testing
# This script runs a coarse AMIP simulation for a short period of time.
# This is meant to be used by upstream packages to test if changes to those
# packages affect the interfaces used by ClimaCoupler.jl.
# Note that this test doesn't check stability or correctness of the model.

# Specify the config file and job id
push!(ARGS, "--config_file", joinpath(@__DIR__, "amip_test.yml"))
push!(ARGS, "--job_id", "amip_test")

# Run the AMIP test
include("../run_amip.jl")
