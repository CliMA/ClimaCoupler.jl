# coupler_driver
# don't forget to run with threading: julia --project --threads 8 (MPI not that useful for debugging coarse runs)
using Pkg
import SciMLBase: step!
using OrdinaryDiffEq
using OrdinaryDiffEq: ODEProblem, solve, SSPRK33, savevalues!, Euler
using LinearAlgebra
import Test: @test
using ClimaCore.Utilities: half, PlusHalf
Pkg.add(PackageSpec(name = "ClimaCore", version = "0.10.3"))

# Get the paths to the necessary data files
include("artifacts.jl")

# Load up a bunch of utils
include("coupler_utils/flux_calculator.jl")
include("coupler_utils/conservation_checker.jl")
include("coupler_utils/regridder.jl")
include("coupler_utils/masker.jl")
include("coupler_utils/general_helper.jl")

# initiate spatial and temporal info
debug_mode = true
t_end = debug_mode ? 100e2 : 2592000 * 1
Δt_cpl = 2e2
saveat = debug_mode ? Δt_cpl * 1 : Δt_cpl * 100

tspan = (0, t_end)

# init MPI
include("mpi/mpi_init.jl")

# init atmos model component
include("atmos/atmos_init1.jl")
# To change FT - we think, go to driver_new.jl and change it by hand, then continue on!
# To turn radiation to gray - coupler_atmos and in driver_new
include("atmos/atmos_init2.jl")
atmos_sim = atmos_init(FT, Y, spaces, integrator, params = params);


@show debug_mode

# init a 2D bounary space at the surface, assuming the same instance (and MPI distribution if applicable) as the atmos domain above
boundary_space = ClimaCore.Fields.level(atmos_sim.domain.face_space, half) # global surface grid

# init land-sea mask
mask = LandSeaMask(FT, mask_data, "LSMASK", boundary_space) # TODO: split up the nc file to individual times for faster computation
# init land model components
Pkg.add( url = "https://github.com/CliMA/ClimaLSM.jl", rev = "albedo_models")
include("bucket/bucket_init.jl")
bucket_sim = bucket_init(FT, FT.(tspan); dt = FT(Δt_cpl), space = boundary_space, saveat = FT(saveat));

# Include ocean and ice models.
include("slab/slab_utils.jl")
include("slab_ocean/slab_init.jl")
prescribed_sst = false
if prescribed_sst == true
    SST = ncreader_rll_to_cgll_from_space(FT, sst_data, "SST", boundary_space)  # a sample SST field from https://gdex.ucar.edu/dataset/158_asphilli.html
    SST = swap_space!(SST, axes(mask)) .* (abs.(mask .- 1)) .+ FT(273.15) # TODO: avoids the "space not the same instance" error
    ocean_params = OceanSlabParameters(FT(20), FT(1500.0), FT(800.0), FT(280.0), FT(1e-3), FT(1e-5))
else
    slab_ocean_sim = slab_ocean_init(FT, tspan, dt = Δt_cpl, space = boundary_space, saveat = saveat)
    ocean_params = nothing
    SST = nothing
end

include("slab_ice/slab_init.jl")
prescribed_sic = true
SIC = ncreader_rll_to_cgll_from_space(FT, sic_data, "SEAICE", boundary_space)
SIC = swap_space!(SIC, axes(mask)) .* (abs.(mask .- 1))

slab_ice_sim = slab_ice_init(FT, tspan, dt = Δt_cpl, space = boundary_space, saveat = saveat, prescribed_sic = SIC)

#load push/pull fields and methods
include("./push_pull.jl")

# Get ready for the integration

# init conservation info collector
CS = OnlineConservationCheck([], [],[],[], [],[],[],[])
# init coupling
coupler_sim = CouplerSimulation(Δt_cpl, integrator.t, boundary_space, FT, mask)
atmos_pull!(atmos_sim, slab_ice_sim, bucket_sim, slab_ocean_sim, mask, boundary_space, prescribed_sst, z0m_S,  z0b_S, T_S, ocean_params, SST)
atmos_push!(atmos_sim, boundary_space, F_A, F_E, F_R, parsed_args)
bucket_pull!(bucket_sim, F_A, F_E, F_R, ρ_sfc)
reinit!(atmos_sim.integrator)
reinit!(bucket_sim.integrator)
if (prescribed_sst !== true) && (prescribed_sic == true)
    ocean_pull!(slab_ocean_sim, F_A, F_R)
    reinit!(slab_ocean_sim.integrator)
end
ice_pull!(slab_ice_sim, F_A, F_R)
reinit!(slab_ice_sim.integrator)

if !is_distributed && (@isdefined CS)
        check_conservation(CS, coupler_sim, atmos_sim, bucket_sim, slab_ocean_sim, slab_ice_sim, F_A .+ F_R)
end
# At this stage, the integrators all have dY(0) computed based p(0), Y(0), t(0)
@show "Starting coupling loop"
walltime = @elapsed for t in (tspan[1]+Δt_cpl:Δt_cpl:tspan[end])
    @show t
    ## Atmos
    atmos_pull!(atmos_sim, slab_ice_sim, bucket_sim, slab_ocean_sim, mask, boundary_space, prescribed_sst, z0m_S,  z0b_S, T_S, ocean_params, SST);
    step!(atmos_sim.integrator, t - atmos_sim.integrator.t, true); # NOTE: instead of Δt_cpl, to avoid accumulating roundoff error
 
    #clip TODO: this is bad!! > limiters
    parent(atmos_sim.integrator.u.c.ρq_tot) .= heaviside.(parent(atmos_sim.integrator.u.c.ρq_tot)); # negligible for total energy cons

    atmos_push!(atmos_sim, boundary_space, F_A, F_E, F_R, dF_A, parsed_args);

    ## Bucket Land
    bucket_pull!(bucket_sim, F_A, F_E, F_R, ρ_sfc);
    step!(bucket_sim.integrator, t - bucket_sim.integrator.t, true);

    ## Slab ocean
    if (prescribed_sst !== true) && (prescribed_sic == true)
        ocean_pull!(slab_ocean_sim, F_A, F_R)
        step!(slab_ocean_sim.integrator, t - slab_ocean_sim.integrator.t, true)
     end

    ## Slab ice
    ice_pull!(slab_ice_sim, F_A, F_R)
    step!(slab_ice_sim.integrator, t - slab_ice_sim.integrator.t, true)

    if !is_distributed && (@isdefined CS)
        check_conservation(CS, coupler_sim, atmos_sim, bucket_sim, slab_ocean_sim, slab_ice_sim, F_A .+ F_R)
    end

end

@show walltime
@show "Postprocessing"
diff_ρe_tot_atmos = CS.ρe_tot_atmos .- CS.ρe_tot_atmos[1]
diff_ρe_tot_slab = (CS.ρe_tot_land .- CS.ρe_tot_land[1])
diff_ρe_tot_slab_ocean = (CS.ρe_tot_ocean .- CS.ρe_tot_ocean[1])
diff_ρe_tot_slab_seaice = (CS.ρe_tot_seaice .- CS.ρe_tot_seaice[1])
times = tspan[1]:coupler_sim.Δt:tspan[end]
plot1 = Plots.plot(times,diff_ρe_tot_atmos, label = "atmos")
Plots.plot!(times,diff_ρe_tot_slab, label = "land")
Plots.plot!(times,diff_ρe_tot_slab_ocean, label = "ocean")
Plots.plot!(times,diff_ρe_tot_slab_seaice, label = "seaice")
tot = CS.ρe_tot_atmos .+ CS.ρe_tot_ocean .+ CS.ρe_tot_land .+ CS.ρe_tot_seaice

Plots.plot!(times, 
    tot .- tot[1],
    label = "tot",
    xlabel = "time [s]",
    ylabel = "energy(t) - energy(t=0) [J]",
        xticks = (collect(1:length(times))[1:50:end], times[1:50:end]),
            )
plot2 = Plots.plot(times, (tot .- tot[1]) ./ tot[1], ylabel = "|dE_earth|/E_earth", label = "",xlabel = "time [s]")

plot3 =Plots.plot(times[1:end-1], abs.((CS.ρe_tot_atmos[2:end] .- CS.ρe_tot_atmos[1:end-1]) ./ Δt_cpl .- (CS.F_energy_ocean[1:end-1] .+ CS.F_energy_land[1:end-1].+ CS.F_energy_ice[1:end-1])), label = "atmos")
Plots.plot!(times[1:end-1], abs.((CS.ρe_tot_land[2:end] .- CS.ρe_tot_land[1:end-1]) ./ Δt_cpl .+ CS.F_energy_land[1:end-1]), label = "land")
Plots.plot!(times[1:end-1], abs.((CS.ρe_tot_ocean[2:end] .- CS.ρe_tot_ocean[1:end-1]) ./ Δt_cpl .+ CS.F_energy_ocean[1:end-1]), label = "ocean")
#Plots.plot!(times[1:end-1], (CS.ρe_tot_seaice[2:end] .- CS.ρe_tot_seaice[1:end-1]) ./ Δt_cpl .+ CS.F_energy_ice[1:end-1], label = "ice")

plot!(yaxis = :log)
plot!(title = "abs{[E(t+dt) - E(t)]/dt + ∑FdA}")
plot(plot1,plot2, plot3)
savefig("conservation_land_ocean_atmos.png")

# collect solutions

# animations
if debug == false
    sol_atm = atmos_sim.integrator.sol
    sol_slab = bucket_sim.integrator.sol
    sol_slab_ice = slab_ice_sim.integrator.sol
    sol_slab_ocean = prescribed_sst !== true ? slab_ocean_sim.integrator.sol : nothing
    
    include("mpi/mpi_postprocess.jl")


    include("coupler_utils/viz_explorer.jl")
    plot_anim()
end


# TODO:
# - update MPI, conservation plots 
# - performance checks, incl threading (--threads=8)
# - add prescribable albedo to RRTMGP + ClimaAtmos
# - add in oupler specific abstractions
# - replace heavisides with smooth functions
# - test dynamical ea ice model 
# - test sea ice for conservation with slab ocean

# Next PR
# - keep adding coupler specific interface
# - fluxes: re-enable different ways to calculate / accumulate fluxes (at overy coupler timestep; at every atmos timestep via callback; via specification of an additional variable) 
# - formalize conservation/physical/performance tests: add error threshold and exception, interval, show option, and make a general interface for it
# - SurfaceFluxes: combine LHF and SHF into enthalpy flux formulation to avoid division by zero
# - Temporally varying SSTs/sea ice: https://pcmdi.llnl.gov/mips/amip/details/
