include("mpi/mpi_init.jl") # setup MPI context for distributed runs #hide

# # AMIP Driver
# don't forget to run with threading: julia --project --threads 8 (MPI not that useful for debugging coarse runs) #hide

#=
## Overview

AMIP is a standard experimental protocol of the Program for Climate Model Diagnosis & Intercomparison (PCMDI). 
It is used as a model benchmark for the atmospheric and land model components, while sea-surface temperatures (SST) and sea-ice concentration (SIC)
are prescribed using time-interpolations between monthly observed data. We use standard datafiles with original sources:
- SST and SIC: https://gdex.ucar.edu/dataset/158_asphilli.html 
- land-sea mask: https://www.ncl.ucar.edu/Applications/Data/#cdf

For more information, see the PCMDI's specificarions for [AMIP I](https://pcmdi.github.io/mips/amip/) and [AMIP II](https://pcmdi.github.io/mips/amip2/). 

This driver contains two modes. The full `AMIP` mode and a `SlabPlanet` (all surfaces are thermal slabs) mode. Since `AMIP` is not a closed system, the 
`SlabPlanet` mode is useful for checking conservation properties of the coupling. 

=#

#=
## Initialization
=#

import SciMLBase: step!
using OrdinaryDiffEq
using OrdinaryDiffEq: ODEProblem, solve, SSPRK33, savevalues!, Euler
using LinearAlgebra
import Test: @test
using ClimaCore.Utilities: half, PlusHalf
using ClimaCore: InputOutput
using Dates
using UnPack

include("cli_options.jl")
(s, parsed_args) = parse_commandline()

## read in some parsed command line arguments 
mode_name = parsed_args["mode_name"]
run_name = parsed_args["run_name"]
energy_check = parsed_args["energy_check"]
const FT = parsed_args["FLOAT_TYPE"] == "Float64" ? Float64 : Float32
land_sim_name = "bucket"
t_end = FT(time_to_seconds(parsed_args["t_end"]))
tspan = (0, t_end)
Δt_cpl = FT(parsed_args["dt_cpl"])
saveat = time_to_seconds(parsed_args["dt_save_to_sol"])
date0 = date = DateTime(parsed_args["start_date"], dateformat"yyyymmdd")
mono_surface = parsed_args["mono_surface"]

## enforce coupling at every timestep
parsed_args["dt"] = string(Δt_cpl) * "secs"

import ClimaCoupler
pkg_dir = pkgdir(ClimaCoupler)
COUPLER_OUTPUT_DIR = joinpath(pkg_dir, "experiments/AMIP/moist_mpi_earth/output", joinpath(mode_name, run_name))
!isdir(COUPLER_OUTPUT_DIR) && mkpath(COUPLER_OUTPUT_DIR)
REGRID_DIR = joinpath(COUPLER_OUTPUT_DIR, "regrid_tmp/")
mkpath(REGRID_DIR)
@info COUPLER_OUTPUT_DIR
@info parsed_args

# get the paths to the necessary data files - land sea mask, sst map, sea ice concentration
include(joinpath(pkgdir(ClimaCoupler), "artifacts", "artifact_funcs.jl"))
sst_data = joinpath(sst_dataset_path(), "sst.nc")
sic_data = joinpath(sic_dataset_path(), "sic.nc")
mask_data = joinpath(mask_dataset_path(), "seamask.nc")

## import coupler unitilies
include("coupler_utils/flux_calculator.jl")
include("coupler_utils/conservation_checker.jl")
include("coupler_utils/regridder.jl")
include("coupler_utils/masker.jl")
include("coupler_utils/calendar_timer.jl")
include("coupler_utils/general_helper.jl")
include("coupler_utils/bcfile_reader.jl")
include("coupler_utils/variable_definer.jl")
include("coupler_utils/diagnostics_gatherer.jl")
include("coupler_utils/offline_postprocessor.jl")

#=
## Component Model Initialization
Here we set initial and boundary conditions for each component model.
=#

#=
### Atmosphere
This used the `ClimaAtmos.jl` driver, with parameterization options specified in the command line arguments.
=#
## init atmos model component
include("atmos/atmos_init.jl")
atmos_sim = atmos_init(FT, Y, integrator, params = params);

#=
We use a common Space for all global surfaces. This enables the MPI processes to operate on the same columns in both 
the atmospheric and surface components, so exchanges are parallelized. Note this is only possible when the 
atmosphere and surface are of the same horizontal resolution. 
=#
## init a 2D bounary space at the surface
boundary_space = atmos_sim.domain.face_space.horizontal_space

## init land-sea mask
land_mask = LandSeaMask(FT, comms_ctx, mask_data, "LSMASK", boundary_space, mono = mono_surface)

## init surface (slab) model components
include("slab/slab_utils.jl")
include("bucket/bucket_init.jl")
include("slab/slab_init.jl")
include("slab_ocean/slab_init.jl")
include("slab_ice/slab_init.jl")

#=
### Land
We use `ClimaLSM.jl`'s bucket model.
=#
land_sim =
    bucket_init(FT, FT.(tspan), parsed_args["config"]; dt = FT(Δt_cpl), space = boundary_space, saveat = FT(saveat))

#=
### Ocean and Sea Ice
In the `AMIP` mode, all ocean properties are prescribed from a file, while sea-ice temperatures are calculated using observed
SIC and assuming a 2m thickness of the ice. 

In the `SlabPlanet` mode, all ocean and sea ice are dynamical models, namely thermal slabs, with different parameters. 
=#

@info mode_name
if mode_name == "amip"
    @info "AMIP boundary conditions - do not expect energy conservation"

    ## ocean
    SST_info = bcfile_info_init(
        FT,
        comms_ctx,
        sst_data,
        "SST",
        boundary_space,
        interpolate_daily = true,
        scaling_function = clean_sst, ## convert to Kelvin
        land_mask = land_mask,
        date0 = date0,
        mono = mono_surface,
    )

    update_midmonth_data!(date0, SST_info)
    SST_init = interpolate_midmonth_to_daily(date0, SST_info)
    ocean_params = OceanSlabParameters(FT(20), FT(1500.0), FT(800.0), FT(280.0), FT(1e-3), FT(1e-5), FT(0.06))
    ocean_sim = (;
        integrator = (;
            u = (; T_sfc = SST_init),
            p = (; params = ocean_params, ocean_mask = (FT(1) .- land_mask)),
            SST_info = SST_info,
        )
    )
    ## sea ice
    SIC_info = bcfile_info_init(
        FT,
        comms_ctx,
        sic_data,
        "SEAICE",
        boundary_space,
        interpolate_daily = true,
        scaling_function = clean_sic, ## convert to fractions
        land_mask = land_mask,
        date0 = date0,
        mono = mono_surface,
    )
    update_midmonth_data!(date0, SIC_info)
    SIC_init = interpolate_midmonth_to_daily(date0, SIC_info)
    ice_mask = get_ice_mask.(SIC_init, mono_surface)
    ice_sim = ice_init(FT; tspan = tspan, dt = Δt_cpl, space = boundary_space, saveat = saveat, ice_mask = ice_mask)
    mode_specifics = (; name = mode_name, SST_info = SST_info, SIC_info = SIC_info)

elseif mode_name == "slabplanet"
    ## ocean
    ocean_sim = ocean_init(
        FT;
        tspan = tspan,
        dt = Δt_cpl,
        space = boundary_space,
        saveat = saveat,
        ocean_mask = (FT(1) .- land_mask), ## NB: this ocean mask includes areas covered by sea ice (unlike the one contained in the cs)
    )

    ## sea ice
    ice_sim = (;
        integrator = (;
            u = (; T_sfc = ClimaCore.Fields.ones(boundary_space)),
            p = (; params = ocean_sim.params, ice_mask = ClimaCore.Fields.zeros(boundary_space)),
        )
    )
    mode_specifics = (; name = mode_name, SST_info = nothing, SIC_info = nothing)
end

#=
## Coupler Initialization
The coupler needs to contain exchange information, manage the calendar and be able to access all component models. It can also optionally
save online diagnostics. These are all initialized here and saved in a global `CouplerSimulation` struct, `cs`. 
=#

## coupler exchange fields
coupler_field_names = (:T_S, :z0m_S, :z0b_S, :ρ_sfc, :q_sfc, :albedo, :F_A, :F_E, :F_R, :P_liq, :P_snow)
coupler_fields =
    NamedTuple{coupler_field_names}(ntuple(i -> ClimaCore.Fields.zeros(boundary_space), length(coupler_field_names)))

## model simulations
model_sims = (atmos_sim = atmos_sim, ice_sim = ice_sim, land_sim = land_sim, ocean_sim = ocean_sim);

## dates
dates = (; date = [date], date0 = [date0], date1 = [Dates.firstdayofmonth(date0)])

#=
### Online Diagnostics
User can write custom diagnostics in the `coupler_utils/variable_definer.jl`.`
=#
## 3d diagnostics
monthly_3d_diags_names = (:T, :u, :q_tot)
monthly_3d_diags = (;
    fields = NamedTuple{monthly_3d_diags_names}(
        ntuple(i -> ClimaCore.Fields.zeros(atmos_sim.domain.center_space), length(monthly_3d_diags_names)),
    ),
    ct = [0],
)
## 2d diagnostics
monthly_2d_diags_names = (:precipitation, :toa, :T_sfc)
monthly_2d_diags = (;
    fields = NamedTuple{monthly_2d_diags_names}(
        ntuple(i -> ClimaCore.Fields.zeros(boundary_space), length(monthly_2d_diags_names)),
    ),
    ct = [0],
)

## coupler simulation
cs = CouplerSimulation(
    comms_ctx,
    FT(Δt_cpl),
    integrator.t,
    tspan,
    dates,
    boundary_space,
    FT,
    (; land = land_mask, ocean = zeros(boundary_space), ice = zeros(boundary_space)),
    coupler_fields,
    model_sims,
    mode_specifics,
    parsed_args,
    monthly_3d_diags,
    monthly_2d_diags,
);

#=
## Initial States Exchange 
=#
## share states between models
include("./push_pull.jl")
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
## Initialize Conservation Checks
=#
## init conservation info collector
if !is_distributed && energy_check && mode_name == "slabplanet"
    conservation_check = OnlineConservationCheck([], [], [], [], [], [], [])
    check_conservation(conservation_check, cs)
end

#=
## Coupling Loop
=#
function solve_coupler!(cs, energy_check)
    @info "Starting coupling loop"

    @unpack model_sims, Δt_cpl, tspan = cs
    @unpack atmos_sim, land_sim, ocean_sim, ice_sim = model_sims

    ## step in time
    walltime = @elapsed for t in ((tspan[1] + Δt_cpl):Δt_cpl:tspan[end])

        cs.dates.date[1] = current_date(cs, t) # if not global, `date` is not updated.

        ## print date on the first of month 
        @calendar_callback :(@show(cs.dates.date[1])) cs.dates.date[1] cs.dates.date1[1]

        if cs.mode.name == "amip"

            ## monthly read of boundary condition data for sea surface temperature (SST) and sea ice concentration (SIC)
            @calendar_callback :(update_midmonth_data!(cs.dates.date[1], cs.mode.SST_info)) cs.dates.date[1] next_date_in_file(
                cs.mode.SST_info,
            )
            SST = ocean_sim.integrator.u.T_sfc .= interpolate_midmonth_to_daily(cs.dates.date[1], cs.mode.SST_info)
            @calendar_callback :(update_midmonth_data!(cs.dates.date[1], cs.mode.SIC_info)) cs.dates.date[1] next_date_in_file(
                cs.mode.SIC_info,
            )
            SIC = interpolate_midmonth_to_daily(cs.dates.date[1], cs.mode.SIC_info)

            ice_mask = ice_sim.integrator.p.ice_mask .= get_ice_mask.(SIC_init, mono_surface)

            ## accumulate diagnostics at each timestep
            accumulate_diags(collect_diags(cs, propertynames(cs.monthly_3d_diags.fields)), cs.monthly_3d_diags)
            accumulate_diags(collect_diags(cs, propertynames(cs.monthly_2d_diags.fields)), cs.monthly_2d_diags)

            ## save and reset monthly averages
            @calendar_callback :(
                map(x -> x ./= cs.monthly_3d_diags.ct[1], cs.monthly_3d_diags.fields),
                save_hdf5(
                    cs.comms_ctx,
                    cs.monthly_3d_diags.fields,
                    cs.dates.date[1],
                    COUPLER_OUTPUT_DIR,
                    name_tag = "3d_",
                ),
                map(x -> x .= FT(0), cs.monthly_3d_diags.fields),
                cs.monthly_3d_diags.ct .= FT(0),
            ) cs.dates.date[1] cs.dates.date1[1]
            @calendar_callback :(
                map(x -> x ./= cs.monthly_2d_diags.ct[1], cs.monthly_2d_diags.fields),
                save_hdf5(
                    cs.comms_ctx,
                    cs.monthly_2d_diags.fields,
                    cs.dates.date[1],
                    COUPLER_OUTPUT_DIR,
                    name_tag = "2d_",
                ),
                map(x -> x .= FT(0), cs.monthly_2d_diags.fields),
                cs.monthly_2d_diags.ct .= FT(0),
            ) cs.dates.date[1] cs.dates.date1[1]

        end

        ## run component models sequentially for one coupling timestep (Δt_cpl)
        ## 1. atmos
        ClimaComms.barrier(comms_ctx)

        atmos_pull!(cs)
        step!(atmos_sim.integrator, t - atmos_sim.integrator.t, true) # NOTE: instead of Δt_cpl, to avoid accumulating roundoff error
        atmos_push!(cs)

        ## 2. land
        land_pull!(cs)
        step!(land_sim.integrator, t - land_sim.integrator.t, true)

        ## 3. ocean
        if cs.mode.name == "slabplanet"
            ocean_pull!(cs)
            step!(ocean_sim.integrator, t - ocean_sim.integrator.t, true)
        end

        ## 4. sea ice
        if cs.mode.name == "amip"
            ice_pull!(cs)
            step!(ice_sim.integrator, t - ice_sim.integrator.t, true)
        end

        ## compute global energy
        if !simulation.is_distributed && energy_check && cs.mode.name == "slabplanet"
            check_conservation(conservation_check, cs)
        end

        ## step to the next calendar month
        @calendar_callback :(cs.dates.date1[1] += Dates.Month(1)) cs.dates.date[1] cs.dates.date1[1]

    end
    @show walltime

    return cs
end

## run the coupled simulation
solve_coupler!(cs, energy_check);

#=
## Postprocessing 
Currently all postprocessing is performed using the root process only. 
=#

if ClimaComms.iamroot(comms_ctx)
    isdir(COUPLER_OUTPUT_DIR * "_artifacts") ? nothing : mkpath(COUPLER_OUTPUT_DIR * "_artifacts")

    ## energy check plots
    if !is_distributed && energy_check && cs.mode.name == "slabplanet"
        @info "Energy Check"
        plot_global_energy(
            conservation_check,
            cs,
            joinpath(COUPLER_OUTPUT_DIR * "_artifacts", "total_energy_bucket.png"),
            joinpath(COUPLER_OUTPUT_DIR * "_artifacts", "total_energy_log_bucket.png"),
        )
    end

    ## sample animations
    if !is_distributed && parsed_args["anim"]
        @info "Animations"
        include("coupler_utils/viz_explorer.jl")
        plot_anim(cs, COUPLER_OUTPUT_DIR * "_artifacts")
    end

    ## plotting AMIP results
    if cs.mode.name == "amip"
        @info "AMIP plots"

        include("coupler_utils/plotter.jl")

        ## ClimaESM
        include("coupler_utils/amip_visualizer.jl")
        post_spec = (;
            T = (:regridded_3d, :zonal_mean),
            u = (:regridded_3d, :zonal_mean),
            q_tot = (:regridded_3d, :zonal_mean),
            toa = (:regridded_2d, :horizontal_2d),
            precipitation = (:regridded_2d, :horizontal_2d),
            T_sfc = (:regridded_2d, :horizontal_2d),
        )

        plot_spec = (;
            T = (; clims = (190, 320), units = "K"),
            u = (; clims = (-50, 50), units = "m/s"),
            q_tot = (; clims = (0, 50), units = "g/kg"),
            toa = (; clims = (-250, 210), units = "W/m^2"),
            precipitation = (clims = (0, 1e-6), units = "kg/m^2/s"),
            T_sfc = (clims = (225, 310), units = "K"),
        )
        amip_paperplots(
            post_spec,
            plot_spec,
            COUPLER_OUTPUT_DIR,
            files_root = ".monthly",
            output_dir = COUPLER_OUTPUT_DIR * "_artifacts",
        )

        ## NCEP reanalysis
        @info "NCEP plots"
        include("coupler_utils/ncep_visualizer.jl")
        ncep_post_spec = (;
            T = (:zonal_mean,),
            u = (:zonal_mean,),
            q_tot = (:zonal_mean,),
            toa = (:horizontal_2d,),
            precipitation = (:horizontal_2d,),
            T_sfc = (:horizontal_2d,),
        )
        ncep_plot_spec = plot_spec
        ncep_paperplots(
            ncep_post_spec,
            ncep_plot_spec,
            COUPLER_OUTPUT_DIR,
            output_dir = COUPLER_OUTPUT_DIR * "_artifacts",
            month_date = cs.dates.date[1],
        ) ## plot data that correspond to the model's last save_hdf5 call (i.e., last month)
    end

    ## clean up
    rm(COUPLER_OUTPUT_DIR; recursive = true, force = true)
end

#=
## Temporary Unit Tests 
To be moved to `test/`
=#
if !is_distributed && cs.mode.name == "amip"
    @info "Unit Tests"
    include("coupler_utils/unit_tester.jl")
end
