## AMIP testing
# This script runs a coarse AMIP simulation for a short period of time.
# This is meant to be used by upstream packages to test if changes to those
# packages affect the interfaces used by ClimaCoupler.jl.
# Note that this test doesn't check stability or correctness of the model.

driver_path = joinpath("..", "run_amip.jl")
config_file = "amip_test.yml"

# Run the AMIP test
run(`julia --project $driver_path --config_file=$config_file --job_id="amip_test"`)
