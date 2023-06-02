# atmos_init: for ClimaAtmos pre-AMIP interface
import ClimaAtmos

driver_file = joinpath(pkgdir(ClimaAtmos), "examples", "hybrid", "driver.jl")
ENV["CI_PERF_SKIP_RUN"] = true
try
    include(driver_file)
catch err
    if err.error !== :exit_profile
        rethrow(err.error)
    end
end
# the clima atmos `integrator` is now defined
struct AtmosSimulation{P, Y, D, I} <: AtmosModelSimulation
    params::P
    Y_init::Y
    domain::D
    integrator::I
end

function atmos_init(::Type{FT}, Y, integrator; params = nothing) where {FT}
    center_space = axes(Y.c.ρe_tot)
    face_space = axes(Y.f.w)
    spaces = (; center_space = center_space, face_space = face_space)
    if :ρe_int in propertynames(Y.c)
        @warn("Running with ρe_int in coupled mode is not tested yet.")
    end

    AtmosSimulation(params, Y, spaces, integrator)
end
