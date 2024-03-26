#=
## Postprocessing
Currently all postprocessing is performed using the root process only.
=#

if ClimaComms.iamroot(comms_ctx)

    ## energy check plots
    if !isnothing(cs.conservation_checks) && cs.mode.name[1:10] == "slabplanet"
        @info "Conservation Check Plots"
        plot_global_conservation(
            cs.conservation_checks.energy,
            cs,
            config_dict["conservation_softfail"],
            figname1 = joinpath(COUPLER_ARTIFACTS_DIR, "total_energy_bucket.png"),
            figname2 = joinpath(COUPLER_ARTIFACTS_DIR, "total_energy_log_bucket.png"),
        )
        plot_global_conservation(
            cs.conservation_checks.water,
            cs,
            config_dict["conservation_softfail"],
            figname1 = joinpath(COUPLER_ARTIFACTS_DIR, "total_water_bucket.png"),
            figname2 = joinpath(COUPLER_ARTIFACTS_DIR, "total_water_log_bucket.png"),
        )
    end

    ## sample animations (not compatible with MPI)
    if !CA.is_distributed(comms_ctx) && config_dict["anim"]
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

        # Compare against observations
        if t_end > 84600
            @info "Error against observations"
            include("user_io/leaderboard.jl")
            compare_vars = ["pr"]
            output_path = joinpath(COUPLER_ARTIFACTS_DIR, "biases.png")
            Leaderboard.plot_biases(atmos_sim.integrator.p.output_dir, compare_vars, cs.dates.date; output_path)
        end
    end

    if isinteractive()
        ## clean up for interactive runs, retain all output otherwise
        rm(COUPLER_OUTPUT_DIR; recursive = true, force = true)

        ## plot all model states and coupler fields (useful for debugging)
        debug(cs)
    end

end