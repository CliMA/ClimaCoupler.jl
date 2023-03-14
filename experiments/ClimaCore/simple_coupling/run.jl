import SciMLBase: step!
using OrdinaryDiffEq
using OrdinaryDiffEq: ODEProblem, solve, SSPRK33, savevalues!, Euler
using LinearAlgebra
import Test: @test
using Dates
using UnPack
using Plots

using ClimaCore.Utilities: half, PlusHalf
using ClimaCore: InputOutput, Fields, Geometry, Spaces, Meshes, Domains, Operators

using ClimaComms

import ClimaCoupler
import ClimaCoupler.Regridder
import ClimaCoupler.ConservationChecker
import ClimaCoupler.BCReader
import ClimaCoupler.TimeManager
import ClimaCoupler.Diagnostics
import ClimaCoupler.PostProcessor

pkg_dir = pkgdir(ClimaCoupler)
COUPLER_OUTPUT_DIR = joinpath(pkg_dir, "experiments/ClimaCore/simple_coupling/output")
mkpath(COUPLER_OUTPUT_DIR)

const FT = Float64
comms_context = ClimaComms.SingletonCommsContext()

# init component models. These are adapter functions for
# - init outputting ModelSimulation
# - setting parameters
include("components/column/init.jl")
col_sim = init(DiffusiveColumn{FT}())
boundary_space = axes(col_sim.Y_init.integrated_flux)

include("components/surface/init.jl")
slab_sim = init(ThermalSlab{FT}(), space = boundary_space)

# read in exchange instructions
include("components/field_exchange.jl")

land_mask = 0.5

#=
### Online Diagnostics
User can write custom diagnostics in the `user_diagnostics.jl`.
=#
sample_diags = init_diagnostics(
    (:T_sfc),
    boundary_space;
    output_dir = COUPLER_OUTPUT_DIR,
    name_tag = "sample_",
)

## coupler simulation
cs = _CoupledSimulation(comms_context, boundary_space,
    model_sims = (;col = col_sim, slab = slab_sim), # model simulations
    coupler_field_names = (:T_S, :F_A), # coupler exchange fields (defined in field exchanger)
    clock = _Clock((1,10)..., 0, dt),
    # diagnostics = (;sample_diags, ),
    conservation_checks = nothing
    )

coupler_build_cache!(cs)

#=
## Initial States Exchange
=#
## share states between models
include("components/push_pull.jl")
atmos_pull!(cs)
parsed_args["ode_algo"] == "ARS343" ? step!(atmos_sim.integrator, Δt_cpl, true) : nothing
atmos_push!(cs)
land_pull!(cs)

## reinitialize (TODO: avoid with interfaces)
reinit!(atmos_sim.integrator)
reinit!(land_sim.integrator)
mode_name == "amip" ? (ice_pull!(cs), reinit!(ice_sim.integrator)) : nothing
mode_name == "slabplanet" ? (ocean_pull!(cs), reinit!(ocean_sim.integrator)) : nothing

#=
## Coupling Loop
=#
coupler_step()