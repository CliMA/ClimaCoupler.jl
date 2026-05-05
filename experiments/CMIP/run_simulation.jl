# # CMIP Driver

#=
## Overview

CMIP is a standard experimental protocol of the World Climate Research Programme (WCRP).
It is used as a model benchmark for coupled, prognostic atmospheric, land, ocean, and sea ice model components.

For more information, see the WCRP's specifications for [CMIP](https://wcrp-cmip.org).

## Running the CMIP configuration
To run a coupled simulation in the default CMIP configuration, run the
following command from the root directory of the repository:
```bash
julia --project=experiments/CMIP experiments/CMIP/run_simulation.jl --config_file config/ci_configs/cmip_oceananigans_climaseaice.yml
```

## Configuration
You can also specify a custom configuration file to run the coupled simulation
in a different setup. The configuration file should be a TOML file that overwrites
the input fields specified in the ClimaCoupler Input module.
A set of example configuration files can be found in the `config/ci_configs/` directory.

To run the coupled simulation interactively with a different configuration file,
set the `config_file` variable in this script to be the path to that file.

For more details about running a coupled simulation, including how to run in a
Slabplanet configuration, please see our documentation.
=#

# Load the necessary modules to run the coupled simulation
include("code_loading.jl")

## using GeometryOps
## using ConservativeRegridding
## 
## @eval GeometryOps.UnitSpherical.RobustCrossProduct begin
##     function robust_cross_product(a::AbstractVector, b::AbstractVector)
##         result, was_stable = stable_cross_product(a, b)
##         was_stable && return normalize(result)
##         a == b && return find_orthogonal(a)
##         if isDoubleFloatsAvailable()
##             result, was_stable = stable_cross_product(to_doublefloat.(a), to_doublefloat.(b))
##             was_stable && return normalize(Float64.(result))
##         end
##         result = exact_cross_product(a, b)
##         return normalize(result)
##     end
## end
## 
## @eval ConservativeRegridding begin
##    function regrid!(dst_field::DenseVector, regridder::Regridder, src_field::DenseVector)
##        T = eltype(regridder.intersections)
## 
##        # Ensure src eltype matches the matrix so CUSPARSE.mul! can dispatch
##        src = if eltype(src_field) == T
##            src_field
##        else
##            regridder.src_temp .= src_field   # in-place convert via the existing work buffer
##            regridder.src_temp
##        end
## 
##        if eltype(dst_field) == T
##            LinearAlgebra.mul!(dst_field, regridder.intersections, src)
##            dst_field ./= regridder.dst_areas
##        else
##            # Run mul! into the matching-eltype work buffer, then convert back
##            LinearAlgebra.mul!(regridder.dst_temp, regridder.intersections, src)
##            regridder.dst_temp ./= regridder.dst_areas
##            dst_field .= regridder.dst_temp
##        end
##        return dst_field
##    end
## end

# Get the configuration file from the command line (or manually set it here)
config_file = "../../config/longrun_configs/cmip_edonly_bucket.yml"

# Set up and run the coupled simulation
cs = CoupledSimulation(config_file)

import Oceananigans as OC
wall_time = Ref(time_ns())
function progress(sim) 
    ocean = sim.model

    (Tmax, Tmin) = extrema(ocean.tracers.T)
    Smax = maximum(ocean.tracers.S)
    Smin = minimum(ocean.tracers.S)
    umax = maximum(ocean.velocities.u)
    vmax = maximum(ocean.velocities.v)
    wmax = maximum(ocean.velocities.w)
    ηmax = maximum(ocean.free_surface.displacement)
    ηmin = minimum(ocean.free_surface.displacement)
    step_time = 1e-9 * (time_ns() - wall_time[])
    @info "time: $(OC.prettytime(sim)), iteration: $(OC.iteration(sim)), Δt: $(OC.prettytime(sim.Δt)), " *
          "extrema(η): ($(round(ηmin, sigdigits=2)), $(round(ηmax, sigdigits=2))) " *
          "extrema(T, S): ($(round(Tmin, digits=2)), $(round(Tmax, digits=2))) ᵒC, " *
          "($(round(Smin, digits=2)), $(round(Smax, digits=2))) psu " *
          "maximum(u): ($(round(umax, sigdigits=2)), $(round(vmax, sigdigits=2)), $(round(wmax, sigdigits=2))) m/s, " *
          "wall time: $(OC.prettytime(step_time))"
    
    wall_time[] = time_ns()

    return nothing
end

# Attaching a progress function to the ocean
OC.add_callback!(cs.model_sims.ocean_sim.ocean, progress, OC.IterationInterval(1))

run!(cs)

# Postprocessing
conservation_softfail = Input.get_coupler_config_dict(config_file)["conservation_softfail"]
rmse_check = Input.get_coupler_config_dict(config_file)["rmse_check"]
postprocess(cs; conservation_softfail, rmse_check)
