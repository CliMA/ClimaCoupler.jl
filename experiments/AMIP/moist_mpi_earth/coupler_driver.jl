# coupler_driver
# don't forget to run with threading: julia --project --threads 8 (MPI not that useful for debugging coarse runs)
using Pkg
#Pkg.rm("ClimaAtmos")
import SciMLBase: step!
using OrdinaryDiffEq
using OrdinaryDiffEq: ODEProblem, solve, SSPRK33, savevalues!, Euler
using LinearAlgebra
import Test: @test
using ClimaCore.Utilities: half, PlusHalf
Pkg.add(PackageSpec(name = "ClimaCore", version = "0.10.3"))

# Get the paths to the necessary data files - land sea mask, sst map, sea ice concentration
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
Δt_cpl = 2e2 # should divide dt_rad with no remainder?
saveat = debug_mode ? Δt_cpl * 1 : Δt_cpl * 100

tspan = (0, t_end)

# init MPI
include("mpi/mpi_init.jl")

# init atmos model component
include("atmos/atmos_init.jl") # FT defined in here
atmos_sim = atmos_init(FT, Y, spaces, integrator, params = params);

# init a 2D bounary space at the surface, assuming the same instance (and MPI distribution if applicable) as the atmos domain above
boundary_space = ClimaCore.Fields.level(atmos_sim.domain.face_space, half) # global surface grid

# read in the land sea mask
#mask = LandSeaMask(FT, mask_data, "LSMASK", boundary_space) # TODO: split up the nc file to individual times for faster computation
mask = zeros(boundary_space)
# Currently, land, sea, and ice all solve their equations everywhere on the globe, and we apply a mask
# to handle boundary conditions correctly, and to handle plotting surface quantities correctly and to
# track conservation quantities.

# This is fine as long as their computations are cheap compared to atmos and as long as there is no
# horizontal flow.

# init land model components
Pkg.add(url = "https://github.com/CliMA/ClimaLSM.jl")
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
    slab_ocean_sim = nothing
else
    slab_ocean_sim = slab_ocean_init(FT, tspan, dt = Δt_cpl, space = boundary_space, saveat = saveat)
    ocean_params = nothing
    SST = nothing
end

include("slab_ice/slab_init.jl")

# read in sea ice concentration and turn into map of where sea ice is present vs not, based on a
# threshold.
#SIC = ncreader_rll_to_cgll_from_space(FT, sic_data, "SEAICE", boundary_space)
#SIC = swap_space!(SIC, axes(mask)) .* (abs.(mask .- 1))
SIC = zeros(boundary_space)
ice_mask = get_ice_mask.(SIC .- FT(25))# here 25% and lower is considered ice free # TODO: generalize to a smaoot function of ice fraction
slab_ice_sim = slab_ice_init(FT, tspan, dt = Δt_cpl, space = boundary_space, saveat = saveat, ice_mask = ice_mask)

#load push/pull fields and methods
include("./push_pull.jl")

# Get ready for the integration

# init conservation info collector
CS = OnlineConservationCheck([], [],[],[], [],[],[],[])
# Reinit after computing the fluxes at t=0.
coupler_sim = CouplerSimulation(Δt_cpl, integrator.t, boundary_space, FT, mask)

atmos_pull!(coupler_sim, atmos_sim, slab_ice_sim, bucket_sim, slab_ocean_sim, boundary_space, prescribed_sst, z0m_S,  z0b_S, T_S, ocean_params, SST, univ_mask)
atmos_push!(atmos_sim, boundary_space, F_A, F_E, F_R, parsed_args)
bucket_pull!(bucket_sim, F_A, F_E, F_R, ρ_sfc)
reinit!(atmos_sim.integrator)
reinit!(bucket_sim.integrator)
if prescribed_sst !== true
    ocean_pull!(slab_ocean_sim, F_A, F_R)
    reinit!(slab_ocean_sim.integrator)
end
ice_pull!(slab_ice_sim, F_A, F_R)
reinit!(slab_ice_sim.integrator)
if !is_distributed && (@isdefined CS)
    check_conservation(CS, coupler_sim, atmos_sim, bucket_sim, slab_ocean_sim, slab_ice_sim, F_A .+ F_R, univ_mask)
end

# At this stage, the integrators all have dY(0) computed based cache(0), Y(0), t(0)
@show "Starting coupling loop"
walltime = @elapsed for t in (tspan[1]+Δt_cpl:Δt_cpl:tspan[end])
    @show t
    ## Atmos
    # sets p = p(0)
    atmos_pull!(coupler_sim, atmos_sim, slab_ice_sim, bucket_sim, slab_ocean_sim, boundary_space, prescribed_sst, z0m_S,  z0b_S, T_S, ocean_params, SST, univ_mask);

    #Y(0) -> Y(1):  Y(1) - Y(0) = dY(0) *dt, then computes dY(1) from Y(1).
    # Note, this would be using p(0) still
    step!(atmos_sim.integrator, t - atmos_sim.integrator.t, true); # NOTE: instead of Δt_cpl, to avoid accumulating roundoff error
 
    #clip TODO: this is bad!! > limiters
    # Maybe we can remove - atmos does not include
    parent(atmos_sim.integrator.u.c.ρq_tot) .= heaviside.(parent(atmos_sim.integrator.u.c.ρq_tot)); # negligible for total energy cons

    atmos_push!(atmos_sim, boundary_space, F_A, F_E, F_R, parsed_args);

    ## Bucket Land
    bucket_pull!(bucket_sim, F_A, F_E, F_R, ρ_sfc);
    step!(bucket_sim.integrator, t - bucket_sim.integrator.t, true);

    ## Slab ocean
    if prescribed_sst !== true
        ocean_pull!(slab_ocean_sim, F_A, F_R)
        step!(slab_ocean_sim.integrator, t - slab_ocean_sim.integrator.t, true)
     end

    ## Slab ice
    ice_pull!(slab_ice_sim, F_A, F_R)
    step!(slab_ice_sim.integrator, t - slab_ice_sim.integrator.t, true)

    if !is_distributed && (@isdefined CS)
        check_conservation(CS, coupler_sim, atmos_sim, bucket_sim, slab_ocean_sim, slab_ice_sim, F_A .+ F_R, univ_mask)
    end

end

@show walltime
@show "Postprocessing"
plot_global_energy(CS, coupler_sim)

diff_ρe_tot_atmos = CS.ρe_tot_atmos .- CS.ρe_tot_atmos[1]
diff_ρe_tot_slab = (CS.ρe_tot_land .- CS.ρe_tot_land[1])
diff_ρe_tot_slab_ocean = (CS.ρe_tot_ocean .- CS.ρe_tot_ocean[1])
diff_ρe_tot_slab_seaice = (CS.ρe_tot_seaice .- CS.ρe_tot_seaice[1])
diff_toa_net_source = (CS.toa_net_source .- CS.toa_net_source[1])
#times = tspan[1]+coupler_sim.Δt:coupler_sim.Δt:tspan[end]
#plot1 = Plots.plot(times,diff_ρe_tot_atmos, label = "atmos")
#Plots.plot!(times,diff_ρe_tot_slab, label = "land")
#Plots.plot!(times,diff_ρe_tot_slab_ocean, label = "ocean")
#Plots.plot!(times,diff_ρe_tot_slab_seaice, label = "seaice")
#Plots.plot!(times, diff_toa_net_source, label = "toa")
tot = CS.ρe_tot_atmos .+ CS.ρe_tot_ocean .+ CS.ρe_tot_land .+ CS.ρe_tot_seaice .+ CS.toa_net_source

#Plots.plot!(times, 
#    tot .- tot[1],
#    label = "tot",
#    xlabel = "time [s]",
#    ylabel = "energy(t) - energy(t=0) [J]")
#savefig("conservation_land_ocean_atmos_ice1_dt.png")
#plot2 = Plots.plot(times, (tot .- tot[1]) ./ tot[1], ylabel = "|dE_earth|/E_earth", label = "",xlabel = "time [s]")

#savefig("conservation_land_ocean_atmos_ice2_dt.png")
#plot3 =Plots.plot(times[1:end-1], abs.((CS.ρe_tot_atmos[2:end] .- CS.ρe_tot_atmos[1:end-1]) ./ Δt_cpl .- (CS.F_energy_ocean[1:end-1] .+ CS.F_energy_land[1:end-1].+ CS.F_energy_ice[1:end-1])), label = "atmos")
#plot3 = Plots.plot(times[1:end-1], abs.((CS.ρe_tot_land[2:end] .- CS.ρe_tot_land[1:end-1]) ./ Δt_cpl .+ CS.dE_energy_#land[1:end-1]), label = "land")
#Plots.plot!(times[1:end-1], abs.((CS.ρe_tot_ocean[2:end] .- CS.ρe_tot_ocean[1:end-1]) ./ Δt_cpl .+ CS.dE_energy_ocean#[1:end-1]), label = "ocean")
#Plots.plot!(times[1:end-1], abs.((CS.ρe_tot_seaice[2:end] .- CS.ρe_tot_seaice[1:end-1]) ./ Δt_cpl .+ CS.dE_energy_ice#[1:end-1]), label = "ice")

#plot!(yaxis = :log, yticks = [1e8,1e10,1e12,1e14, 1e16, 1e18,1e20])
#plot!(title = "abs{[E(t+dt) - E(t)]/dt + ∑FdA}")
#savefig("conservation_land_ocean_atmos_ice3_dt.png")



# collect solutions

# animations
if debug_mode == false
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
