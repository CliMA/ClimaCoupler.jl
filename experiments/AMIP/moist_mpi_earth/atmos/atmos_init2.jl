
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
Pkg.add(PackageSpec(name = "ClimaCore", version = "0.10.3"))
Pkg.pin("ClimaCore")
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
