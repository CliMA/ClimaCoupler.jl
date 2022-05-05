# atmos_init: for ClimaAtmos pre-AMIP interface

# clone locally ClimaAtmos and checkout the checkpoint branch (this is necessary because the driver is currently not a module, and the checkpoint branch avoids the necessity for frequent updates)
CWD = pwd()
ATMOS_DIR = CWD * "/../../../../ClimaAtmos.jl_cpl/"

run(`rm -rf $ATMOS_DIR`)
ClimaComms.barrier(comms_ctx)


if !is_distributed || (is_distributed && ClimaComms.iamroot(comms_ctx))
    run(`git clone https://github.com/CliMA/ClimaAtmos.jl.git $ATMOS_DIR`)
    run(`chmod 775 checkout_ClimaAtmos.sh`)
    run(`."/"checkout_ClimaAtmos.sh $ATMOS_DIR`)
end

ClimaComms.barrier(comms_ctx)

Pkg.develop(path = ATMOS_DIR)
Pkg.activate(joinpath(ATMOS_DIR, "examples"))
Pkg.instantiate()
Pkg.add(PackageSpec(name = "ClimaCore", version = "0.10.3"))

driver_orig = ATMOS_DIR * "examples/hybrid/driver.jl"
driver_new = ATMOS_DIR * "examples/hybrid/driver_new.jl"

# get only the chunks of the driver.jl we need (i.e. before solve!)
if !is_distributed || (is_distributed && ClimaComms.iamroot(comms_ctx))

    run(
        pipeline(
            `awk '$1=="sol" {f=0;next} f{print;next} NR==1 {f=1}'`;
            stdin = driver_orig,
            stdout = driver_new,
            append = false,
        ),
    )

    # remove the hard coded tendency specifications (so that we can define a diffusion tendency that allows for the correct coupled boundary conditions)
    run(`sed -i.bak -e '94,123d' $driver_new`) # yep!  

    # remove hard coded MPI init (optional but neater)
    run(`sed -i.bak -e '98,119d' $driver_new`) # yep!
end

# init model using the modified driver
include(joinpath(ATMOS_DIR, "examples/hybrid/cli_options.jl"))

(s, parsed_args) = parse_commandline()

# add coupler-taylored CA functions
TEST_NAME = "coupled_atmos"

coupler_atmos_file = CWD * "/atmos/" * TEST_NAME * ".jl"

ClimaComms.barrier(comms_ctx)
include(coupler_atmos_file)

# specify sim parameters
const FT = parsed_args["FLOAT_TYPE"] == "Float64" ? Float64 : Float32 # parsed_args["FLOAT_TYPE"] = FT
parsed_args["dt"] = string(Δt_cpl) * "secs"
parsed_args["t_end"] = string(t_end) * "secs"
parsed_args["enable_threading"] = true
parsed_args["dt_save_to_sol"] = string(saveat) * "secs"

atoms_setup_dir = joinpath(ATMOS_DIR, "examples/hybrid/sphere/")

if !is_distributed || (is_distributed && ClimaComms.iamroot(comms_ctx))
    run(`cp $coupler_atmos_file $atoms_setup_dir`)
end

# init model using the modified driver
ClimaComms.barrier(comms_ctx)
include(driver_new) # this stops just before `solve!`

spaces = (; center_space = center_space, face_space = face_space)

struct AtmosSimulation{P, Y, D, I}
    params::P
    Y_init::Y
    domain::D
    integrator::I
end

function atmos_init(::Type{FT}, Y, spaces, integrator; params = nothing) where {FT}
    AtmosSimulation(params, Y, spaces, integrator)
end
