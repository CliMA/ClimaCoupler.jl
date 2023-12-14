using Plots
using ClimaCorePlots
using ClimaCore: Geometry
using Statistics

function plot_anim(cs, out_dir = ".")

    atmos_sim = cs.model_sims.atmos_sim
    slab_land_sim = cs.model_sims.land_sim
    slab_ocean_sim = cs.model_sims.ocean_sim
    slab_ice_sim = cs.model_sims.ice_sim
    mode_name = cs.mode.name
    SST = cs.fields.T_S

    # plot the lowest atmos (center) levels of key variables
    sol_atm = atmos_sim.integrator.sol

    anim = Plots.@animate for u in sol_atm.u
        Plots.plot(
            Fields.level(Geometry.UVVector.(u.c.uₕ).components.data.:1, 5),
            title = "Surface zonal wind",
            xlabel = "",
            ylabel = "cubed sphere panels",
        )
    end
    Plots.mp4(anim, joinpath(out_dir, "anim_u.mp4"), fps = 20)

    anim = Plots.@animate for u in sol_atm.u
        Plots.plot(
            Fields.level(u.c.ρe_tot, 1) .- Fields.level(sol_atm.u[1].c.ρe_tot, 1),
            clims = (-5000, 50000),
            title = "Surface total energy",
            xlabel = "",
            ylabel = "cubed sphere panels",
        )
    end
    Plots.mp4(anim, joinpath(out_dir, "anim_rhoe_anom.mp4"), fps = 20)

    anim = Plots.@animate for u in sol_atm.u
        Plots.plot(
            Fields.level(u.c.ρq_tot ./ u.c.ρ, 1),
            title = "Surface total humidity",
            xlabel = "",
            ylabel = "cubed sphere panels",
        )
    end
    Plots.mp4(anim, joinpath(out_dir, "anim_qt.mp4"), fps = 20)

    # plot combined surfaces
    combined_field = zeros(boundary_space)
    sol_slab = slab_land_sim.integrator.sol

    if mode_name == "slabplanet"
        sol_slab_ocean = slab_ocean_sim.integrator.sol
        anim = Plots.@animate for (bucketu, oceanu) in zip(sol_slab.u, sol_slab_ocean.u)
            land_T_sfc = get_land_temp_from_state(cs.model_sims.land_sim, bucketu)
            combine_surfaces_from_sol!(
                combined_field,
                cs.surface_fractions,
                (; land = land_T_sfc, ocean = oceanu.T_sfc, ice = FT(0)),
            )
            # TODO do anomaly from zonal mean temp
            # global_mean = mean(combined_field)
            # anomaly = combined_field .- global_mean
            Plots.plot(combined_field, title = "Surface temperature [K]", xlabel = "", ylabel = "cubed sphere panels")
        end
    elseif mode_name == "slabplanet_eisenman"
        slab_ice_sim = slab_ice_sim.integrator.sol
        anim = Plots.@animate for (bucketu, iceu) in zip(sol_slab.u, slab_ice_sim.u)
            land_T_sfc = get_land_temp_from_state(cs.model_sims.land_sim, bucketu)
            combine_surfaces_from_sol!(
                combined_field,
                cs.surface_fractions,
                (; land = land_T_sfc, ocean = FT(0), ice = iceu.T_sfc),
            )
            Plots.plot(combined_field, title = "Surface temperature [K])", xlabel = "", ylabel = "cubed sphere panels")
        end

    elseif mode_name == "amip"
        sol_slab_ice = slab_ice_sim.integrator.sol
        anim = Plots.@animate for (bucketu, iceu) in zip(sol_slab.u, sol_slab_ice.u)
            land_T_sfc = get_land_temp_from_state(cs.model_sims.land_sim, bucketu)
            combine_surfaces_from_sol!(
                combined_field,
                cs.surface_fractions,
                (; land = land_T_sfc, ocean = SST, ice = iceu.T_sfc),
            )
            Plots.plot(combined_field, title = "Surface temperature [K])", xlabel = "", ylabel = "cubed sphere panels")
        end
    end
    Plots.mp4(anim, joinpath(out_dir, "earth_T.mp4"), fps = 20)

    combined_field = zeros(boundary_space)
    anim = Plots.@animate for bucketu in sol_slab.u
        combine_surfaces_from_sol!(
            combined_field,
            cs.surface_fractions,
            (; land = bucketu.bucket.W, ocean = 0.0, ice = 0.0),
        )
        Plots.plot(combined_field, title = "Sub-surface water storage [m]", xlabel = "", ylabel = "cubed sphere panels")
    end
    Plots.mp4(anim, joinpath(out_dir, "bucket_W.mp4"), fps = 20)

    combined_field = zeros(boundary_space)
    anim = Plots.@animate for bucketu in sol_slab.u
        combine_surfaces_from_sol!(
            combined_field,
            cs.surface_fractions,
            (; land = bucketu.bucket.σS, ocean = 0.0, ice = 0.0),
        )
        Plots.plot(combined_field, title = "Snow cover fraction over land", xlabel = "", ylabel = "cubed sphere panels")
    end
    Plots.mp4(anim, joinpath(out_dir, "bucket_snow.mp4"), fps = 20)

    if mode_name == "slabplanet_eisenman"
        sol_ice = cs.model_sims.ice_sim.integrator.sol
        combined_field = zeros(boundary_space)
        anim = Plots.@animate for sol_iceu in sol_ice.u
            combine_surfaces_from_sol!(
                combined_field,
                cs.surface_fractions,
                (; land = 0.0, ocean = 0.0, ice = sol_iceu.h_ice),
            )
            Plots.plot(combined_field)
        end
        Plots.mp4(anim, joinpath(out_dir, "eisenman_seaice.mp4"), fps = 20)
    end
    # plot surface fluxes
    # TODO as part of the flux accumulation PR
end
