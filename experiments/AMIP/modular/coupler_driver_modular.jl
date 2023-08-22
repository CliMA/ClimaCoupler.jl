include("mpi/mpi_init.jl") # setup MPI context for distributed runs #hide

# # AMIP Driver

#=
## Overview

AMIP is a standard experimental protocol of the Program for Climate Model Diagnosis & Intercomparison (PCMDI).
It is used as a model benchmark for the atmospheric and land model components, while sea-surface temperatures (SST) and sea-ice concentration (SIC)
are prescribed using time-interpolations between monthly observed data. We use standard data files with original sources:
- SST and SIC: https://gdex.ucar.edu/dataset/158_asphilli.html
- land-sea mask: https://www.ncl.ucar.edu/Applications/Data/#cdf

For more information, see the PCMDI's specifications for [AMIP I](https://pcmdi.github.io/mips/amip/) and [AMIP II](https://pcmdi.github.io/mips/amip2/).

This driver contains two modes. The full `AMIP` mode and a `SlabPlanet` (all surfaces are thermal slabs) mode. Since `AMIP` is not a closed system, the
`SlabPlanet` mode is useful for checking conservation properties of the coupling.

=#

#=
## Start Up
Before starting Julia, ensure your environment is properly set up:
```julia
module purge
module load julia/1.8.1 openmpi/4.1.1 hdf5/1.12.1-ompi411 #netcdf-c/4.6.1

export CLIMACORE_DISTRIBUTED="MPI" #include if using MPI, otherwise leave empty
export JUlIA_MPI_BINARY="system"
export JULIA_HDF5_PATH=""
```

Next instantiate/build all packages listed in `Manifest.toml`:
```julia
julia --project -e 'using Pkg; Pkg.instantiate(); Pkg.build()'
julia --project -e 'using Pkg; Pkg.build("MPI"); Pkg.build("HDF5")'
```

The `coupler_driver.jl` is now ready to be run. You can run a SLURM job (e.g., run `sbatch sbatch_job.sh` from the terminal), or
you can run directly from the Julia REPL. The latter is recommended for debugging of lightweight simulations, and should be run
with threading enabled:
```julia
julia --project --threads 8
```
=#

#=
## Initialization
Here we import standard Julia packages, ClimaESM packages, parse in command-line arguments (if none are specified then the defaults in `cli_options.jl` apply).
We then specify the input data file names. If these are not already downloaded, include `artifacts/download_artifacts.jl`.
=#

import SciMLBase: step!, reinit!
using OrdinaryDiffEq
using OrdinaryDiffEq: ODEProblem, solve, SSPRK33, savevalues!, Euler
using LinearAlgebra
import Test: @test
using Dates
using UnPack
using Plots
using Statistics: mean

using ClimaCore.Utilities: half, PlusHalf
using ClimaCore: InputOutput, Fields
import ClimaCore.Spaces as Spaces

if !(@isdefined parsed_args)
    include("cli_options.jl")
    parsed_args = parse_commandline(argparse_settings())
end

## modify parsed args for fast testing from REPL #hide
if isinteractive()
    parsed_args["coupled"] = true #hide
    parsed_args["surface_setup"] = "PrescribedSurface" #hide # necessary to stop Atmos from calculating its own surface fluxes
    parsed_args["moist"] = "equil" #hide
    parsed_args["vert_diff"] = true #hide
    parsed_args["rad"] = "gray" #hide
    parsed_args["energy_check"] = true #hide
    parsed_args["mode_name"] = "slabplanet" #hide
    parsed_args["t_end"] = "10days" #hide
    parsed_args["dt_save_to_sol"] = "3600secs" #hide
    parsed_args["dt_cpl"] = 200 #hide
    parsed_args["dt"] = "200secs" #hide
    parsed_args["mono_surface"] = true #hide
    parsed_args["turb_flux_partition"] = "CombinedStateFluxes" #hide
    parsed_args["h_elem"] = 4 #hide
    # parsed_args["dt_save_restart"] = "5days" #hide
    parsed_args["precip_model"] = "0M" #hide
    parsed_args["job_id"] = "interactive_debug_run" #hide
    parsed_args["monthly_checkpoint"] = true #hide
    parsed_args["uft"] = "Businger" #hide
    parsed_args["anim"] = true #hide
end

## read in some parsed command line arguments
mode_name = parsed_args["mode_name"]
run_name = parsed_args["run_name"]
energy_check = parsed_args["energy_check"]
if !(@isdefined FT)
    const FT = parsed_args["FLOAT_TYPE"] == "Float64" ? Float64 : Float32
end
land_sim_name = "bucket"
t_end = Int(time_to_seconds(parsed_args["t_end"]))
tspan = (Int(0), t_end)
Δt_cpl = Int(parsed_args["dt_cpl"])
saveat = time_to_seconds(parsed_args["dt_save_to_sol"])
date0 = date = DateTime(parsed_args["start_date"], dateformat"yyyymmdd")
mono_surface = parsed_args["mono_surface"]
monthly_checkpoint = parsed_args["monthly_checkpoint"]
restart_dir = parsed_args["restart_dir"]
restart_t = Int(parsed_args["restart_t"])

import ClimaCoupler
import ClimaCoupler.Regridder
import ClimaCoupler.Regridder:
    update_surface_fractions!, combine_surfaces!, combine_surfaces_from_sol!, dummmy_remap!, binary_mask
import ClimaCoupler.ConservationChecker:
    EnergyConservationCheck, WaterConservationCheck, check_conservation!, plot_global_conservation
import ClimaCoupler.Utilities: CoupledSimulation, float_type, swap_space!
import ClimaCoupler.BCReader:
    bcfile_info_init, float_type_bcf, update_midmonth_data!, next_date_in_file, interpolate_midmonth_to_daily
import ClimaCoupler.TimeManager: current_date, datetime_to_strdate, trigger_callback, Monthly, EveryTimestep
import ClimaCoupler.Diagnostics: get_var, init_diagnostics, accumulate_diagnostics!, save_diagnostics, TimeMean
import ClimaCoupler.PostProcessor: postprocess

import ClimaCoupler.Interfacer:
    AtmosModelSimulation,
    SurfaceModelSimulation,
    SurfaceStub,
    SeaIceModelSimulation,
    LandModelSimulation,
    OceanModelSimulation,
    get_field,
    update_field!,
    update_sim!
import ClimaCoupler.FluxCalculator:
    PartitionedStateFluxes,
    CombinedStateFluxes,
    combined_turbulent_fluxes!,
    MoninObukhovScheme,
    partitioned_turbulent_fluxes!
import ClimaCoupler.FieldExchanger:
    import_atmos_fields!,
    import_combined_surface_fields!,
    update_sim!,
    update_model_sims!,
    reinit_model_sims!,
    step_model_sims!
import ClimaCoupler.Checkpointer: checkpoint_model_state, get_model_state_vector, restart_model_state!

pkg_dir = pkgdir(ClimaCoupler)
COUPLER_OUTPUT_DIR = joinpath(pkg_dir, "experiments/AMIP/modular/output", joinpath(mode_name, run_name))
mkpath(COUPLER_OUTPUT_DIR)

REGRID_DIR = joinpath(COUPLER_OUTPUT_DIR, "regrid_tmp/")
mkpath(REGRID_DIR)

COUPLER_ARTIFACTS_DIR = COUPLER_OUTPUT_DIR * "_artifacts"
isdir(COUPLER_ARTIFACTS_DIR) ? nothing : mkpath(COUPLER_ARTIFACTS_DIR)

@info COUPLER_OUTPUT_DIR
@info parsed_args

## get the paths to the necessary data files: land-sea mask, sst map, sea ice concentration
include(joinpath(pkgdir(ClimaCoupler), "artifacts", "artifact_funcs.jl"))
sst_data = joinpath(sst_dataset_path(), "sst.nc")
sic_data = joinpath(sic_dataset_path(), "sic.nc")
land_mask_data = joinpath(mask_dataset_path(), "seamask.nc")

## user-specified diagnostics
include("user_io/user_diagnostics.jl")

#=
## Component Model Initialization
Here we set initial and boundary conditions for each component model.
=#

#=
### Atmosphere
This uses the `ClimaAtmos.jl` driver, with parameterization options specified in the command line arguments.
=#
## init atmos model component
include("components/atmosphere/climaatmos_init.jl")
atmos_sim = atmos_init(FT, parsed_args);
thermo_params = get_thermo_params(atmos_sim) # TODO: this should be shared by all models

#=
We use a common `Space` for all global surfaces. This enables the MPI processes to operate on the same columns in both
the atmospheric and surface components, so exchanges are parallelized. Note this is only possible when the
atmosphere and surface are of the same horizontal resolution.
=#
## init a 2D boundary space at the surface
boundary_space = atmos_sim.domain.face_space.horizontal_space

# init land-sea fraction
land_fraction =
    Regridder.land_fraction(FT, REGRID_DIR, comms_ctx, land_mask_data, "LSMASK", boundary_space, mono = mono_surface)

## init surface (slab) model components
include("components/land/bucket_init.jl")
include("components/land/bucket_utils.jl")
include("components/ocean/slab_ocean_init.jl")
include("components/ocean/slab_seaice_init.jl")

#=
### Land
We use `ClimaLSM.jl`'s bucket model.
=#
land_sim = bucket_init(
    FT,
    FT.(tspan),
    parsed_args["config"],
    parsed_args["albedo_from_file"],
    comms_ctx,
    REGRID_DIR;
    dt = FT(Δt_cpl),
    space = boundary_space,
    saveat = FT(saveat),
    area_fraction = land_fraction,
)

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
        REGRID_DIR,
        sst_data,
        "SST",
        boundary_space,
        comms_ctx,
        interpolate_daily = true,
        scaling_function = clean_sst, ## convert to Kelvin
        land_fraction = land_fraction,
        date0 = date0,
        mono = mono_surface,
    )

    update_midmonth_data!(date0, SST_info)
    SST_init = interpolate_midmonth_to_daily(date0, SST_info)
    ocean_sim = SurfaceStub((;
        T_sfc = SST_init,
        ρ_sfc = ClimaCore.Fields.zeros(boundary_space),
        z0m = FT(1e-3),
        z0b = FT(1e-3),
        beta = FT(1),
        α = FT(0.06),
        area_fraction = (FT(1) .- land_fraction),
        phase = TD.Liquid(),
        thermo_params = thermo_params,
    ))
    ## sea ice
    SIC_info = bcfile_info_init(
        FT,
        REGRID_DIR,
        sic_data,
        "SEAICE",
        boundary_space,
        comms_ctx,
        interpolate_daily = true,
        scaling_function = clean_sic, ## convert to fraction
        land_fraction = land_fraction,
        date0 = date0,
        mono = mono_surface,
    )
    update_midmonth_data!(date0, SIC_info)
    SIC_init = interpolate_midmonth_to_daily(date0, SIC_info)
    ice_fraction = get_ice_fraction.(SIC_init, mono_surface)
    ice_sim = ice_init(
        FT;
        tspan = tspan,
        dt = Δt_cpl,
        space = boundary_space,
        saveat = saveat,
        area_fraction = ice_fraction,
        thermo_params = thermo_params,
    )
    mode_specifics = (; name = mode_name, SST_info = SST_info, SIC_info = SIC_info)

elseif mode_name == "slabplanet"
    ## ocean
    ocean_sim = ocean_init(
        FT;
        tspan = tspan,
        dt = Δt_cpl,
        space = boundary_space,
        saveat = saveat,
        area_fraction = (FT(1) .- land_fraction), ## NB: this ocean fraction includes areas covered by sea ice (unlike the one contained in the cs)
        thermo_params = thermo_params,
    )

    ## sea ice (here set to zero area coverage)
    ice_sim = SurfaceStub((;
        T_sfc = ClimaCore.Fields.ones(boundary_space),
        ρ_sfc = ClimaCore.Fields.zeros(boundary_space),
        z0m = FT(0),
        z0b = FT(0),
        beta = FT(1),
        α = FT(1),
        area_fraction = ClimaCore.Fields.zeros(boundary_space),
        phase = TD.Ice(),
        thermo_params = thermo_params,
    ))

    mode_specifics = (; name = mode_name, SST_info = nothing, SIC_info = nothing)
end

#=
## Coupler Initialization
The coupler needs to contain exchange information, manage the calendar and be able to access all component models. It can also optionally
save online diagnostics. These are all initialized here and saved in a global `CouplerSimulation` struct, `cs`.
=#

## coupler exchange fields
coupler_field_names = (
    :T_S,
    :z0m_S,
    :z0b_S,
    :ρ_sfc,
    :q_sfc,
    :albedo,
    :beta,
    :F_turb_energy,
    :F_turb_moisture,
    :F_turb_ρτxz,
    :F_turb_ρτyz,
    :F_radiative,
    :P_liq,
    :P_snow,
    :F_radiative_TOA,
    :P_net,
)
coupler_fields =
    NamedTuple{coupler_field_names}(ntuple(i -> ClimaCore.Fields.zeros(boundary_space), length(coupler_field_names)))

## model simulations
model_sims = (atmos_sim = atmos_sim, ice_sim = ice_sim, land_sim = land_sim, ocean_sim = ocean_sim);

## dates
dates = (; date = [date], date0 = [date0], date1 = [Dates.firstdayofmonth(date0)], new_month = [false])

#=
### Online Diagnostics
User can write custom diagnostics in the `user_diagnostics.jl`.
=#
monthly_3d_diags = init_diagnostics(
    (:T, :u, :q_tot, :q_liq_ice),
    atmos_sim.domain.center_space;
    save = Monthly(),
    operations = (; accumulate = TimeMean([Int(0)])),
    output_dir = COUPLER_OUTPUT_DIR,
    name_tag = "monthly_mean_3d_",
)

monthly_2d_diags = init_diagnostics(
    (:precipitation_rate, :toa_fluxes, :T_sfc, :tubulent_energy_fluxes),
    boundary_space;
    save = Monthly(),
    operations = (; accumulate = TimeMean([Int(0)])),
    output_dir = COUPLER_OUTPUT_DIR,
    name_tag = "monthly_mean_2d_",
)

diagnostics = (monthly_3d_diags, monthly_2d_diags)

#=
## Initialize Conservation Checks
=#
## init conservation info collector
conservation_checks = nothing
if energy_check
    @assert(
        mode_name == "slabplanet" && !ClimaAtmos.is_distributed(ClimaComms.context(boundary_space)),
        "Only non-distributed slabplanet allowable for energy_check"
    )
    conservation_checks =
        (; energy = EnergyConservationCheck([], [], [], [], [], []), water = WaterConservationCheck([], [], [], []))
end

## coupler simulation
cs = CoupledSimulation{FT}(
    comms_ctx,
    dates,
    boundary_space,
    coupler_fields,
    parsed_args,
    conservation_checks,
    [tspan[1], tspan[2]],
    atmos_sim.integrator.t,
    Δt_cpl,
    (; land = land_fraction, ocean = zeros(boundary_space), ice = zeros(boundary_space)),
    model_sims,
    mode_specifics,
    diagnostics,
);

#=
## Restart component model states if specified
=#
if restart_dir !== "unspecified"
    for sim in cs.model_sims
        if get_model_state_vector(sim) !== nothing
            restart_model_state!(sim, comms_ctx, restart_t; input_dir = restart_dir)
        end
    end
end

#=
## Initialize Component Model Exchange
=#
turbulent_fluxes = nothing
if parsed_args["turb_flux_partition"] == "PartitionedStateFluxes"
    turbulent_fluxes = PartitionedStateFluxes()
elseif parsed_args["turb_flux_partition"] == "CombinedStateFluxes"
    turbulent_fluxes = CombinedStateFluxes()
else
    error("turb_flux_partition must be either PartitionedStateFluxes or CombinedStateFluxes")
end

# 1) coupler combines surface states and calculates rho_sfc using surface and atmos variables
update_surface_fractions!(cs)
import_combined_surface_fields!(cs.fields, cs.model_sims, cs.boundary_space, turbulent_fluxes)
import_atmos_fields!(cs.fields, cs.model_sims, cs.boundary_space, turbulent_fluxes)
update_model_sims!(cs.model_sims, cs.fields, turbulent_fluxes)

# 2) each surface component model calculates its own vapor specific humidity (q_sfc)
# TODO: the q_sfc calculation follows the design of the bucket q_sfc, but it would be neater to abstract this from step!
step!(land_sim, Δt_cpl)
step!(ocean_sim, Δt_cpl)
step!(ice_sim, Δt_cpl)

# 3) coupler re-imports updated surface fields and calculates turbulent fluxes, while updating atmos sfc_conditions
if turbulent_fluxes isa CombinedStateFluxes
    # calculate fluxes using combined surface states on the atmos grid
    import_combined_surface_fields!(cs.fields, cs.model_sims, cs.boundary_space, turbulent_fluxes) # i.e. T_sfc, albedo, z0, beta, q_sfc
    combined_turbulent_fluxes!(cs.model_sims, cs.fields, turbulent_fluxes) # this updates the atmos thermo state, sfc_ts
elseif turbulent_fluxes isa PartitionedStateFluxes
    # calculate turbulent fluxes in surface models and save the weighted average in coupler fields
    partitioned_turbulent_fluxes!(cs.model_sims, cs.fields, cs.boundary_space, MoninObukhovScheme(), thermo_params)

    # update atmos sfc_conditions for surface temperature
    # TODO: this is hard coded and needs to be simplified (need CA modification)
    new_p = get_new_cache(atmos_sim, cs.fields)
    ClimaAtmos.SurfaceConditions.update_surface_conditions!(atmos_sim.integrator.u, new_p, atmos_sim.integrator.t) # sets T_sfc (but SF calculation not necessary - CA)
    atmos_sim.integrator.p.sfc_conditions .= new_p.sfc_conditions
end

# 4) given the new sfc_conditions, atmos calls the radiative flux callback
reinit_model_sims!(cs.model_sims) # NB: for atmos this sets a nonzero radiation flux

# 5) coupler re-imports updated atmos fluxes (radiative fluxes for both `turbulent_fluxes` types
# and also turbulent fluxes if `turbulent_fluxes isa CombinedStateFluxes`,
# and sends them to the surface component models. If `turbulent_fluxes isa PartitionedStateFluxes`
# atmos receives the turbulent fluxes from the coupler.
import_atmos_fields!(cs.fields, cs.model_sims, cs.boundary_space, turbulent_fluxes)
update_model_sims!(cs.model_sims, cs.fields, turbulent_fluxes)

#=
## Coupling Loop
=#
function solve_coupler!(cs)
    @info "Starting coupling loop"

    @unpack model_sims, Δt_cpl, tspan = cs
    @unpack atmos_sim, land_sim, ocean_sim, ice_sim = model_sims

    ## step in time
    walltime = @elapsed for t in ((tspan[1] + Δt_cpl):Δt_cpl:tspan[end])

        cs.dates.date[1] = current_date(cs, t) # if not global, `date` is not updated.

        ## print date on the first of month
        if cs.dates.date[1] >= cs.dates.date1[1]
            @show(cs.dates.date[1])
        end

        if cs.mode.name == "amip"

            ## monthly read of boundary condition data for SST and SIC
            if cs.dates.date[1] >= next_date_in_file(cs.mode.SST_info)
                update_midmonth_data!(cs.dates.date[1], cs.mode.SST_info)
            end
            update_field!(
                ocean_sim,
                Val(:surface_temperature),
                interpolate_midmonth_to_daily(cs.dates.date[1], cs.mode.SST_info),
            )

            if cs.dates.date[1] >= next_date_in_file(cs.mode.SIC_info)
                update_midmonth_data!(cs.dates.date[1], cs.mode.SIC_info)
            end
            interpolate_midmonth_to_daily(cs.dates.date[1], cs.mode.SIC_info)
            update_field!(ice_sim, Val(:area_fraction), get_ice_fraction.(SIC_init, mono_surface))  # TODO: change to SIC

            ## calculate and accumulate diagnostics at each timestep
            ClimaComms.barrier(comms_ctx)
            accumulate_diagnostics!(cs)

            ## save and reset monthly averages
            save_diagnostics(cs)

        end

        ## compute global energy
        !isnothing(cs.conservation_checks) ? check_conservation!(cs, get_slab_energy, get_land_energy) : nothing

        ## run component models sequentially for one coupling timestep (Δt_cpl)
        ClimaComms.barrier(comms_ctx)
        update_model_sims!(cs.model_sims, cs.fields, turbulent_fluxes)

        ## step sims
        step_model_sims!(cs.model_sims, t)

        ## exchange combined fields and (if specified) calculate fluxes using combined states
        update_surface_fractions!(cs)
        import_combined_surface_fields!(cs.fields, cs.model_sims, cs.boundary_space, turbulent_fluxes) # i.e. T_sfc, albedo, z0, beta
        if turbulent_fluxes isa CombinedStateFluxes
            combined_turbulent_fluxes!(cs.model_sims, cs.fields, turbulent_fluxes) # this updates the surface thermo state, sfc_ts, in ClimaAtmos (but also unnecessarily calculates fluxes)
        elseif turbulent_fluxes isa PartitionedStateFluxes
            # calculate turbulent fluxes in surfaces and save the weighted average in coupler fields
            partitioned_turbulent_fluxes!(cs.model_sims, cs.fields, cs.boundary_space, MoninObukhovScheme(), thermo_params)

            # update atmos sfc_conditions for surface temperature - TODO: this needs to be simplified (need CA modification)
            new_p = get_new_cache(atmos_sim, cs.fields)
            ClimaAtmos.SurfaceConditions.update_surface_conditions!(
                atmos_sim.integrator.u,
                new_p,
                atmos_sim.integrator.t,
            ) # to set T_sfc (but SF calculation not necessary - CA modification)
            atmos_sim.integrator.p.sfc_conditions .= new_p.sfc_conditions
        end

        import_atmos_fields!(cs.fields, cs.model_sims, cs.boundary_space, turbulent_fluxes) # radiative and/or turbulent

        ## monthly callbacks
        if trigger_callback(cs, Monthly())
            ## step to the next calendar month
            cs.dates.date1[1] += Dates.Month(1)
            ## checkpoint model state
            if monthly_checkpoint
                for sim in cs.model_sims
                    if get_model_state_vector(sim) !== nothing
                        checkpoint_model_state(sim, comms_ctx, Int(t), output_dir = COUPLER_ARTIFACTS_DIR)
                    end
                end
            end
        end

    end
    @show walltime

    return cs
end

## exit if running performance anaysis #hide
if haskey(ENV, "CI_PERF_SKIP_COUPLED_RUN") #hide
    throw(:exit_profile_init) #hide
end #hide

## run the coupled simulation
solve_coupler!(cs);

#=
## Postprocessing
Currently all postprocessing is performed using the root process only.
=#

if ClimaComms.iamroot(comms_ctx)

    ## energy check plots
    if !isnothing(cs.conservation_checks) && cs.mode.name == "slabplanet"
        @info "Conservation Check Plots"
        plot_global_conservation(
            cs.conservation_checks.energy,
            cs,
            figname1 = joinpath(COUPLER_ARTIFACTS_DIR, "total_energy_bucket.png"),
            figname2 = joinpath(COUPLER_ARTIFACTS_DIR, "total_energy_log_bucket.png"),
        )
        plot_global_conservation(
            cs.conservation_checks.water,
            cs,
            figname1 = joinpath(COUPLER_ARTIFACTS_DIR, "total_water_bucket.png"),
            figname2 = joinpath(COUPLER_ARTIFACTS_DIR, "total_water_log_bucket.png"),
        )
    end

    ## sample animations
    if !is_distributed && parsed_args["anim"]
        @info "Animations"
        include("user_io/viz_explorer.jl")
        plot_anim(cs, COUPLER_ARTIFACTS_DIR)
    end

    ## plotting AMIP results
    if cs.mode.name == "amip"
        @info "AMIP plots"

        ## ClimaESM
        include("user_io/amip_visualizer.jl")
        post_spec = (;
            T = (:regrid, :zonal_mean),
            u = (:regrid, :zonal_mean),
            q_tot = (:regrid, :zonal_mean),
            toa_fluxes = (:regrid, :horizontal_slice),
            precipitation_rate = (:regrid, :horizontal_slice),
            T_sfc = (:regrid, :horizontal_slice),
            tubulent_energy_fluxes = (:regrid, :horizontal_slice),
            q_liq_ice = (:regrid, :zonal_mean),
        )

        plot_spec = (;
            T = (; clims = (190, 320), units = "K"),
            u = (; clims = (-50, 50), units = "m/s"),
            q_tot = (; clims = (0, 30), units = "g/kg"),
            toa_fluxes = (; clims = (-250, 250), units = "W/m^2"),
            precipitation_rate = (clims = (0, 1e-4), units = "kg/m^2/s"),
            T_sfc = (clims = (225, 310), units = "K"),
            tubulent_energy_fluxes = (; clims = (-250, 250), units = "W/m^2"),
            q_liq_ice = (; clims = (0, 10), units = "g/kg"),
        )
        amip_data = amip_paperplots(
            post_spec,
            plot_spec,
            COUPLER_OUTPUT_DIR,
            files_root = ".monthly",
            output_dir = COUPLER_ARTIFACTS_DIR,
        )

        ## NCEP reanalysis
        @info "NCEP plots"
        include("user_io/ncep_visualizer.jl")
        ncep_post_spec = (;
            T = (:zonal_mean,),
            u = (:zonal_mean,),
            q_tot = (:zonal_mean,),
            toa_fluxes = (:horizontal_slice,),
            precipitation_rate = (:horizontal_slice,),
            T_sfc = (:horizontal_slice,),
            tubulent_energy_fluxes = (:horizontal_slice,),
        )
        ncep_plot_spec = plot_spec
        ncep_data = ncep_paperplots(
            ncep_post_spec,
            ncep_plot_spec,
            COUPLER_OUTPUT_DIR,
            output_dir = COUPLER_ARTIFACTS_DIR,
            month_date = cs.dates.date[1],
        ) ## plot data that correspond to the model's last save_hdf5 call (i.e., last month)
    end

    ## clean up
    # rm(COUPLER_OUTPUT_DIR; recursive = true, force = true)
end
