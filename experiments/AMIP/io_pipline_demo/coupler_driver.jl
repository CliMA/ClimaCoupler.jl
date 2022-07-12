# coupler_driver
# don't forget to run with threading: julia --project --threads 8 (MPI not that useful for debugging coarse runs)

import SciMLBase: step!
using OrdinaryDiffEq
using OrdinaryDiffEq: ODEProblem, solve, SSPRK33, savevalues!, Euler
using LinearAlgebra
import Test: @test
using ClimaCore.Utilities: half, PlusHalf
using Dates

include("cli_options.jl")
(s, parsed_args) = parse_commandline()
# Read in some parsed args
prescribed_sst = parsed_args["prescribed_sst"]
energy_check = parsed_args["energy_check"]
const FT = parsed_args["FLOAT_TYPE"] == "Float64" ? Float64 : Float32
land_sim = "bucket"
t_end = FT(time_to_seconds(parsed_args["t_end"]))
tspan = (0, t_end)
Δt_cpl = FT(parsed_args["dt_cpl"])
saveat = time_to_seconds(parsed_args["dt_save_to_sol"])

date0 = Date(2014,1,29) # first date

# overwrite some parsed args :P
parsed_args["coupled"] = true
parsed_args["dt"] = string(Δt_cpl) * "secs"
parsed_args["enable_threading"] = true
parsed_args["microphy"] = "0M"
parsed_args["forcing"] = nothing
parsed_args["idealized_h2o"] = false
parsed_args["vert_diff"] = true
parsed_args["rad"] = "gray"
parsed_args["hyperdiff"] = true
parsed_args["config"] = "sphere"
parsed_args["moist"] = "equil"

# atmos diagnostics
# parsed_args["dt_save_to_disk"] = "1hours" # saves jld2 at this frequency, can be used as a restart also # this crashed upon atmos init (thermo_params)
# ENV["RESTART_FILE"] = ".jld2"

# Get the paths to the necessary data files - land sea mask, sst map, sea ice concentration
include("artifacts.jl")

# import coupler utils
include("coupler_utils/flux_calculator.jl")
include("coupler_utils/conservation_checker.jl")
include("coupler_utils/regridder.jl")
include("coupler_utils/masker.jl")
include("coupler_utils/general_helper.jl")
include("coupler_utils/timer.jl")
include("coupler_utils/bcfile_reader.jl")

# init MPI
include("mpi/mpi_init.jl")

# init atmos model component
include("atmos/atmos_init.jl")
atmos_sim = atmos_init(FT, Y, integrator, params = params);

# init a 2D bounary space at the surface, assuming the same instance (and MPI distribution if applicable) as the atmos domain above
boundary_space = ClimaCore.Fields.level(atmos_sim.domain.face_space, half) # global surface grid

# init land-sea mask
landmask = LandSeaMask(FT, mask_data, "LSMASK", boundary_space)
landmask = swap_space!(landmask, boundary_space) # needed if we are reading from previous run

# init surface (slab) model components
# we need some types that are defined in these files
include("slab/slab_utils.jl")
include("bucket/bucket_init.jl")
include("slab/slab_init.jl")
include("slab_ocean/slab_init.jl")
include("slab_ice/slab_init.jl")

# land
if land_sim == "slab"
    slab_sim = slab_init(FT; tspan, dt = Δt_cpl, space = boundary_space, saveat = saveat, mask = landmask)
end
if land_sim == "bucket"
    slab_sim = bucket_init(FT, FT.(tspan); dt = FT(Δt_cpl), space = boundary_space, saveat = FT(saveat))
end
prescribed_sst = true
if prescribed_sst
    println("No ocean sim - do not expect energy conservation")
    
    # ocean
    SST_info = bcfile_info_init(sst_data, "SST", boundary_space, segment_idx0 = [Int(1729)],interpolate_monthly = true, scaling_function = clean_sst)
    update_midmonth_data!(date0, SST_info)
    SST = interpolate_midmonth_to_daily(date0, SST_info) 
    ocean_params = OceanSlabParameters(FT(20), FT(1500.0), FT(800.0), FT(280.0), FT(1e-3), FT(1e-5), FT(0.06))
    slab_ocean_sim = nothing
    
    # sea ice
    SIC_info = bcfile_info_init(sic_data, "SEAICE", boundary_space, segment_idx0 = [Int(1729)] ,interpolate_monthly = true, scaling_function = clean_sic)
    update_midmonth_data!(date0, SIC_info)
    SIC = interpolate_midmonth_to_daily(date0, SIC_info) 
    ice_mask = get_ice_mask.(SIC .- FT(50), FT) # here 50% and lower is considered ice free
    slab_ice_sim =
    slab_ice_init(FT; tspan = tspan, dt = Δt_cpl, space = boundary_space, saveat = saveat, ice_mask = ice_mask)
else
    slab_ocean_sim =
        slab_ocean_init(FT; tspan = tspan, dt = Δt_cpl, space = boundary_space, saveat = saveat, mask = landmask)
    SST = nothing
    ocean_params = nothing
end

# init coupler
coupler_sim = CouplerSimulation(FT(Δt_cpl), integrator.t, boundary_space, FT, landmask)
include("./push_pull.jl")
# init conservation info collector
atmos_pull!(
    atmos_sim,
    slab_ice_sim,
    slab_sim,
    slab_ocean_sim,
    boundary_space,
    prescribed_sst,
    z0m_S,
    z0b_S,
    T_S,
    ρ_sfc,
    q_sfc,
    ocean_params,
    SST,
    landmask,
)

atmos_push!(atmos_sim, boundary_space, F_A, F_R, F_E, P_liq, parsed_args)
land_pull!(slab_sim, F_A, F_R, F_E, P_liq, ρ_sfc)
ice_pull!(slab_ice_sim, F_A, F_R)
if !prescribed_sst
    ocean_pull!(slab_ocean_sim, F_A, F_R)
    reinit!(slab_ocean_sim.integrator)
end

reinit!(atmos_sim.integrator)
reinit!(slab_sim.integrator)
reinit!(slab_ice_sim.integrator)
if !is_distributed && energy_check && !prescribed_sst
    CS = OnlineConservationCheck([], [], [], [], [], [])
    check_conservation(CS, coupler_sim, atmos_sim, slab_sim, slab_ocean_sim, slab_ice_sim)
end

# coupling loop
@show "Starting coupling loop"
t = 0
walltime = @elapsed for t in ((tspan[1] + Δt_cpl):Δt_cpl:tspan[end])
    #@show t

    global date = @current_date t 

    ## BC load # not a noticeable slowdown with @elapsed

    # load monthly data from files if needed 
    @calendar_callback :(update_midmonth_data!(date, SST_info)) date next_month_date(SST_info)
    @calendar_callback :(update_midmonth_data!(date, SIC_info)) date next_month_date(SIC_info)

    # old approach
    # Dates.days(date - SST_info.all_dates[SST_info.segment_idx[1] + Int(1)]) < FT(0) ? nothing : (update_midmonth_data!(date, SST_info) , @show ("yes:$(date) vs $(SST_info.all_dates[SST_info.segment_idx[1] + Int(1)])")) # TODO: abstract: calendar_callback
    # Dates.days(date - SIC_info.all_dates[SIC_info.segment_idx[1] + Int(1)]) < FT(0) ? nothing : update_midmonth_data!(date, SIC_info) 

    # daily interpolation from monthly data
    SST = interpolate_midmonth_to_daily(date, SST_info)
    SIC = interpolate_midmonth_to_daily(date, SIC_info)
    
    slab_ice_sim.integrator.p.Ya.ice_mask .= get_ice_mask.(SIC .- FT(50), FT) 

    ## Atmos
    atmos_pull!(
        atmos_sim,
        slab_ice_sim,
        slab_sim,
        slab_ocean_sim,
        boundary_space,
        prescribed_sst,
        z0m_S,
        z0b_S,
        T_S,
        ρ_sfc,
        q_sfc,
        ocean_params,
        SST,
        landmask,
    )
    step!(atmos_sim.integrator, t - atmos_sim.integrator.t, true) # NOTE: instead of Δt_cpl, to avoid accumulating roundoff error

    atmos_push!(atmos_sim, boundary_space, F_A, F_R, F_E, P_liq, parsed_args)

    ## Slab land
    land_pull!(slab_sim, F_A, F_R, F_E, P_liq, ρ_sfc)
    step!(slab_sim.integrator, t - slab_sim.integrator.t, true)

    ## Slab ocean
    if !prescribed_sst
        ocean_pull!(slab_ocean_sim, F_A, F_R)
        step!(slab_ocean_sim.integrator, t - slab_ocean_sim.integrator.t, true)
    end

    ## Slab ice
    ice_pull!(slab_ice_sim, F_A, F_R)
    step!(slab_ice_sim.integrator, t - slab_ice_sim.integrator.t, true)

    ## Compute energy
    if !is_distributed && energy_check && !prescribed_sst
        check_conservation(CS, coupler_sim, atmos_sim, slab_sim, slab_ocean_sim, slab_ice_sim)
    end
end

@show walltime

@show "Postprocessing"
if energy_check && !prescribed_sst
    plot_global_energy(CS, coupler_sim, "total_energy_bucket.png", "total_energy_log_bucket.png")
end


@show P_liq
# # animations
if (land_sim == "bucket") && parsed_args["anim"]
    #make it so this works with slab land?
    include("coupler_utils/viz_explorer.jl")
    plot_anim(atmos_sim, slab_sim, slab_ocean_sim, slab_ice_sim, landmask, prescribed_sst, SST)
end

# Cleanup temporary files
# TODO: Where should this live?
rm(REGRID_DIR; recursive = true, force = true)




