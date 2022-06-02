CWD = pwd()
ATMOS_DIR = CWD * "../ClimaAtmos.jl_cpl/"
driver_new = ATMOS_DIR * "examples/hybrid/driver_new.jl"
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
