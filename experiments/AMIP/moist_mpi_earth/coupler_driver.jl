# coupler_driver
# don't forget to run with threading: julia --project --threads 8 (MPI not that useful for debugging coarse runs)

# import packages
# Note we are setting SIC = 0 because we get energy errors when it is not.
using Pkg
import SciMLBase: step!
using OrdinaryDiffEq
using OrdinaryDiffEq: ODEProblem, solve, SSPRK33, savevalues!, Euler
using LinearAlgebra
import Test: @test
using ClimaCore.Utilities: half, PlusHalf
Pkg.add(PackageSpec(name = "ClimaCore", version = "0.10.3"))
## When is FT set? It is already being used (as a global) in flux_calculator
# import coupler utils
include("coupler_utils/flux_calculator.jl")
include("coupler_utils/conservation_checker.jl")
include("coupler_utils/regridder.jl")
include("coupler_utils/masker.jl")
include("coupler_utils/general_helper.jl")

# # initiate spatial and temporal info
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
infile = "data/seamask.nc"
mask = LandSeaMask(FT, infile, "LSMASK", boundary_space) # TODO: split up the nc file to individual times for faster computation
#mask = ones(boundary_space)
# init surface (slab) model components
# ClimaLSM unregistered:
Pkg.add( url = "https://github.com/CliMA/ClimaLSM.jl", rev = "albedo_models")
#Pkg.develop(path="../../../../ClimaLSM.jl")

include("bucket/bucket_init.jl")
bucket_sim = bucket_init(FT, FT.(tspan); dt = FT(Δt_cpl), space = boundary_space, saveat = FT(saveat));

include("slab/slab_utils.jl")

include("slab_ocean/slab_init.jl")
prescribed_sst = false
if prescribed_sst == true
    SST = ncreader_rll_to_cgll_from_space(FT, "data/sst.nc", "SST", boundary_space)  # a sample SST field from https://gdex.ucar.edu/dataset/158_asphilli.html
    SST = swap_space!(SST, axes(mask)) .* (abs.(mask .- 1)) .+ FT(273.15) # TODO: avoids the "space not the same instance" error
    ocean_params = OceanSlabParameters(FT(20), FT(1500.0), FT(800.0), FT(280.0), FT(1e-3), FT(1e-5))
else
    slab_ocean_sim = slab_ocean_init(FT, tspan, dt = Δt_cpl, space = boundary_space, saveat = saveat, mask = mask)
    SST = nothing
    ocean_params = nothing
end

include("slab_ice/slab_init.jl")
prescribed_sic = true
if prescribed_sic == true
    SIC = zeros(boundary_space)
    # SIC = ncreader_rll_to_cgll_from_space(FT, "data/sic.nc", "SEAICE", boundary_space)
#    SIC = swap_space!(SIC, axes(mask)) .* (abs.(mask .- 1)) # zero over land, -# over ocean
     slab_ice_sim = slab_ice_init(FT, tspan, dt = Δt_cpl, space = boundary_space, saveat = saveat, prescribed_sic = SIC)
else
    slab_ice_sim = slab_ice_init(
        FT,
        tspan,
        dt = Δt_cpl,
        space = boundary_space,
        saveat = saveat,
        ocean_params = slab_ocean_sim.integrator.p.params,
    )
end

# init coupler's boundary fields for regridding (TODO: technically this can be bypassed by directly rigridding on model grids)
T_S = ClimaCore.Fields.zeros(boundary_space) # temperature
z0m_S = ClimaCore.Fields.zeros(boundary_space)
z0b_S = ClimaCore.Fields.zeros(boundary_space)
ρ_sfc = ClimaCore.Fields.zeros(boundary_space)
q_sfc = ClimaCore.Fields.zeros(boundary_space)

F_A = ClimaCore.Fields.zeros(boundary_space) # aerodynamic turbulent fluxes
F_E = ClimaCore.Fields.zeros(boundary_space) # evaporation due to turbulent fluxes
F_R = ClimaCore.Fields.zeros(boundary_space) # radiative fluxes
dF_A = ClimaCore.Fields.zeros(boundary_space) # aerodynamic turbulent fluxes

# init conservation info collector
CS = OnlineConservationCheck([], [],[],[])
# coupling methods
function atmos_push!(atmos_sim, boundary_space, F_A, F_E, F_R, dF_A, parsed_args)
    F_A .= ClimaCore.Fields.zeros(boundary_space)
    dummmy_remap!(F_A, atmos_sim.integrator.p.dif_flux_energy)
    F_E .= ClimaCore.Fields.zeros(boundary_space)
    dummmy_remap!(F_E, atmos_sim.integrator.p.dif_flux_ρq_tot)
    F_R .= ClimaCore.Fields.zeros(boundary_space)
    parsed_args["rad"] == "gray" ? dummmy_remap!(F_R, level(atmos_sim.integrator.p.ᶠradiation_flux, half)) : nothing # TODO: albedo hard coded...
    dF_A .= ClimaCore.Fields.zeros(boundary_space)
    dummmy_remap!(dF_A, atmos_sim.integrator.p.∂F_aero∂T_sfc)
end

function bucket_pull!(bucket_sim, F_A, F_E, F_R, ρ_sfc)
    @. bucket_sim.integrator.p.bucket.ρ_sfc = ρ_sfc
    @. bucket_sim.integrator.p.bucket.SHF = F_A
    @. bucket_sim.integrator.p.bucket.LHF = FT(0.0)
    @. bucket_sim.integrator.p.bucket.E = F_E ./ ρ_cloud_liq(bucket_sim.params.earth_param_set)
    @. bucket_sim.integrator.p.bucket.R_n = F_R
end

function ocean_pull!(slab_ocean_sim, F_A, F_R)
    @. slab_ocean_sim.integrator.p.F_aero = -F_A
    @. slab_ocean_sim.integrator.p.F_rad = -F_R
end

function ice_pull!(slab_ice_sim, F_A, F_R, dF_A)
    @. slab_ice_sim.integrator.p.F_aero = -F_A
    @. slab_ice_sim.integrator.p.F_rad = -F_R
    @. slab_ice_sim.integrator.p.∂F_aero∂T_sfc = dF_A
end

function atmos_pull!(atmos_sim, slab_ice_sim, bucket_sim, slab_ocean_sim, mask, boundary_space, prescribed_sst,  z0m_S,  z0b_S, T_S, ocean_params, SST)
    combined_field = zeros(boundary_space)
    # coupler_get: T_sfc, z_0m, z_0b
    if prescribed_sst == true
        parent(combined_field) .=
            combine_surface.(
                parent(mask) .- parent(slab_ice_sim.integrator.p.ice_mask .* FT(2)),
                parent(bucket_sim.integrator.u.bucket.T_sfc),
                parent(SST),
                parent(slab_ice_sim.integrator.u.T_sfc),
            ) # prescribed SSTs
        dummmy_remap!(T_S, combined_field)
        parent(combined_field) .=
            combine_surface.(
                parent(mask),
                parent(bucket_sim.params.z_0m .* mask),
                parent(ocean_params.z0m .* (abs.(mask .- 1))),
            )
        dummmy_remap!(z0m_S, combined_field)
        parent(combined_field) .=
            combine_surface.(
                parent(mask),
                parent(bucket_sim.params.z_0b .* mask),
                parent(ocean_params.z0b .* (abs.(mask .- 1))),
            )
        dummmy_remap!(z0b_S, combined_field)
    else
        parent(combined_field) .=
            combine_surface.(
                parent(mask) .- parent(slab_ice_sim.integrator.p.ice_mask .* FT(2)),
                parent(bucket_sim.integrator.u.bucket.T_sfc),
                parent(slab_ocean_sim.integrator.u.T_sfc),
                parent(slab_ice_sim.integrator.u.T_sfc),
            )
        dummmy_remap!(T_S, combined_field)
        parent(combined_field) .=
            combine_surface.(
                parent(mask),
                parent(bucket_sim.params.z_0m .* mask),
                parent(slab_ocean_sim.integrator.p.params.z0m .* (abs.(mask .- 1))),
            )
        dummmy_remap!(z0m_S, combined_field)
        parent(combined_field) .=
            combine_surface.(
                parent(mask),
                parent(bucket_sim.params.z_0b .* mask),
                parent(slab_ocean_sim.integrator.p.params.z0b .* (abs.(mask .- 1))),
            )
        dummmy_remap!(z0b_S, combined_field)
        
    end
    # Compute ρ_sfc based atmos properties at lowest level
    set_ρ_sfc!(ρ_sfc, T_S, atmos_sim.integrator)
    # Now compute ocean and sea ice q_sat
    ocean_q_sfc = TD.q_vap_saturation_generic.(atmos_sim.integrator.p.params, slab_ocean_sim.integrator.u.T_sfc, ρ_sfc, TD.Liquid())
    sea_ice_q_sfc = TD.q_vap_saturation_generic.(atmos_sim.integrator.p.params, slab_ice_sim.integrator.u.T_sfc, ρ_sfc, TD.Ice())
    # Pull q_sfc from land, and compute q_sfc on surface with it and the above computed values.
    parent(combined_field) .=
        combine_surface.(
            parent(mask) .- parent(slab_ice_sim.integrator.p.ice_mask .* FT(2)),
            parent(bucket_sim.integrator.p.bucket.q_sfc),
            parent(ocean_q_sfc),
            parent(sea_ice_q_sfc),
        )
    dummmy_remap!(q_sfc, combined_field)

    atmos_sim.integrator.p.rrtmgp_model.surface_temperature .= field2array(T_S) # supplied to atmos for radiation
    # add albedo here
    coords = ClimaCore.Fields.coordinate_field(axes(bucket_sim.integrator.u.bucket.S))
    α_land = surface_albedo.(Ref(bucket_sim.params.albedo), coords, bucket_sim.integrator.u.bucket.S, bucket_sim.params.S_c)

    # calculate turbulent fluxes on atmos grid and save in atmos cache
    info_sfc = (; T_sfc = T_S, ρ_sfc = ρ_sfc, q_sfc = q_sfc, z0m = z0m_S, z0b = z0b_S, ice_mask = slab_ice_sim.integrator.p.ice_mask)
    # This should also need q_sfc - currently it assumes q_sfc = q_sat
    calculate_surface_fluxes_atmos_grid!(atmos_sim.integrator, info_sfc)
end
########
# init coupling
coupler_sim = CouplerSimulation(Δt_cpl, integrator.t, boundary_space, FT, mask)
atmos_pull!(atmos_sim, slab_ice_sim, bucket_sim, slab_ocean_sim, mask, boundary_space, prescribed_sst, z0m_S,  z0b_S, T_S, ocean_params, SST)
atmos_push!(atmos_sim, boundary_space, F_A, F_E, F_R, dF_A, parsed_args)
bucket_pull!(bucket_sim, F_A, F_E, F_R, ρ_sfc)
reinit!(atmos_sim.integrator)
reinit!(bucket_sim.integrator)
if (prescribed_sst !== true) && (prescribed_sic == true)
    ocean_pull!(slab_ocean_sim, F_A, F_R)
    reinit!(slab_ocean_sim.integrator)
end
ice_pull!(slab_ice_sim, F_A, F_R, dF_A)
reinit!(slab_ice_sim.integrator)

if !is_distributed && (@isdefined CS)
        check_conservation(CS, coupler_sim, atmos_sim, bucket_sim, slab_ocean_sim, slab_ice_sim)
end
# At this stage, the integrators all have dY(0) computed based p(0), Y(0), t(0)
@show "Starting coupling loop"
walltime = @elapsed for t in (tspan[1]+Δt_cpl:Δt_cpl:tspan[end])
    @show t
    ## Atmos
    atmos_pull!(atmos_sim, slab_ice_sim, bucket_sim, slab_ocean_sim, mask, boundary_space, prescribed_sst, z0m_S,  z0b_S, T_S, ocean_params, SST);
    # set ρ_sfc(T_sfc(i-1)). Get fluxes based on T_sfc(i-1), q_sfc(i)
    # Compute Y(i) from dY(i-1). comute dY(i) from Y(i) and p(i-1).
    step!(atmos_sim.integrator, t - atmos_sim.integrator.t, true); # NOTE: instead of Δt_cpl, to avoid accumulating roundoff error
 
    #clip TODO: this is bad!! > limiters
    parent(atmos_sim.integrator.u.c.ρq_tot) .= heaviside.(parent(atmos_sim.integrator.u.c.ρq_tot)); # negligible for total energy cons

    # coupler_push!: eventually, atmos.p will have accumulate energy and water fluxes, which will be pushed to F_a, F_E, etc.
    # Right now, these are the same as atmos.p from before the atmos step!
    atmos_push!(atmos_sim, boundary_space, F_A, F_E, F_R, dF_A, parsed_args);

    ## Bucket Land
    bucket_pull!(bucket_sim, F_A, F_E, F_R, ρ_sfc);
    # Compute Y(i) from dY(i-1). q_sfc computed based on T_sfc(i), ρ_sfc(i-1). compute dY(i) from Y(i) and mixed p(i-1) and p(i) 
    step!(bucket_sim.integrator, t - bucket_sim.integrator.t, true);

    ## Slab ocean
    if (prescribed_sst !== true) && (prescribed_sic == true)
        ocean_pull!(slab_ocean_sim, F_A, F_R)
        step!(slab_ocean_sim.integrator, t - slab_ocean_sim.integrator.t, true)
     end
    ## Slab ice
    step!(slab_ice_sim.integrator, t - slab_ice_sim.integrator.t, true)

    if !is_distributed && (@isdefined CS)
        check_conservation(CS, coupler_sim, atmos_sim, bucket_sim, slab_ocean_sim, slab_ice_sim)
    end

end

@show walltime
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
plot2 = Plots.plot(times, (tot .- tot[1]) ./ tot[1], ylabel = "Fractional Energy Error", label = "",xlabel = "time [s]")
plot(plot1,plot2)
savefig("conservation_land_ocean_atmos.png")
@show "Postprocessing"
# collect solutions
sol_atm = atmos_sim.integrator.sol
sol_slab = bucket_sim.integrator.sol
sol_slab_ice = slab_ice_sim.integrator.sol
sol_slab_ocean = prescribed_sst !== true ? slab_ocean_sim.integrator.sol : nothing

include("mpi/mpi_postprocess.jl")


# animations
if debug == false
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