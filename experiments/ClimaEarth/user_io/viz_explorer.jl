import Plots
import ClimaCorePlots
import ClimaCore as CC
import ClimaCoupler: Regridder

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
        Plots.plot(CC.Fields.level(CC.Geometry.UVVector.(u.c.uₕ).components.data.:1, 5))
    end
    Plots.mp4(anim, joinpath(out_dir, "anim_u.mp4"), fps = 10)

    anim = Plots.@animate for u in sol_atm.u
        Plots.plot(CC.Fields.level(u.c.ρe_tot, 1) .- CC.Fields.level(sol_atm.u[1].c.ρe_tot, 1), clims = (-5000, 50000))
    end
    Plots.mp4(anim, joinpath(out_dir, "anim_rhoe_anom.mp4"), fps = 10)

    anim = Plots.@animate for u in sol_atm.u
        Plots.plot(CC.Fields.level(u.c.ρq_tot ./ u.c.ρ, 1))
    end
    Plots.mp4(anim, joinpath(out_dir, "anim_qt.mp4"), fps = 10)

    # get all surface fractions and package into a named tuple
    land_fraction = Interfacer.get_field(land_sim, Val(:area_fraction))
    ocean_fraction = Interfacer.get_field(ocean_sim, Val(:area_fraction))
    ice_fraction = Interfacer.get_field(ice_sim, Val(:area_fraction))
    surface_fractions = (; land = land_fraction, ocean = ocean_fraction, ice = ice_fraction)

    # plot combined surfaces
    combined_field = zeros(boundary_space)
    sol_slab = slab_land_sim.integrator.sol

    if mode_name == "slabplanet"
        sol_slab_ocean = slab_ocean_sim.integrator.sol
        anim = Plots.@animate for (bucketu, oceanu) in zip(sol_slab.u, sol_slab_ocean.u)
            land_T_sfc = get_land_temp_from_state(cs.model_sims.land_sim, bucketu)
            Regridder.combine_surfaces_from_sol!(
                combined_field,
                surface_fractions,
                (; land = land_T_sfc, ocean = oceanu.T_sfc, ice = FT(0)),
            )
            Plots.plot(combined_field)
        end
    elseif mode_name == "slabplanet_eisenman"
        slab_ice_sim = slab_ice_sim.integrator.sol
        anim = Plots.@animate for (bucketu, iceu) in zip(sol_slab.u, slab_ice_sim.u)
            land_T_sfc = get_land_temp_from_state(cs.model_sims.land_sim, bucketu)
            Regridder.combine_surfaces_from_sol!(
                combined_field,
                surface_fractions,
                (; land = land_T_sfc, ocean = FT(0), ice = iceu.T_sfc),
            )
            Plots.plot(combined_field)
        end

    elseif mode_name == "amip"
        sol_slab_ice = slab_ice_sim.integrator.sol
        anim = Plots.@animate for (bucketu, iceu) in zip(sol_slab.u, sol_slab_ice.u)
            land_T_sfc = get_land_temp_from_state(cs.model_sims.land_sim, bucketu)
            Regridder.combine_surfaces_from_sol!(
                combined_field,
                surface_fractions,
                (; land = land_T_sfc, ocean = SST, ice = iceu.T_sfc),
            )
            Plots.plot(combined_field)
        end
    end
    Plots.mp4(anim, joinpath(out_dir, "earth_T.mp4"), fps = 10)

    combined_field = zeros(boundary_space)
    anim = Plots.@animate for bucketu in sol_slab.u
        Regridder.combine_surfaces_from_sol!(
            combined_field,
            surface_fractions,
            (; land = bucketu.bucket.W, ocean = 0.0, ice = 0.0),
        )
        Plots.plot(combined_field)
    end
    Plots.mp4(anim, joinpath(out_dir, "bucket_W.mp4"), fps = 10)

    combined_field = zeros(boundary_space)
    anim = Plots.@animate for bucketu in sol_slab.u
        Regridder.combine_surfaces_from_sol!(
            combined_field,
            surface_fractions,
            (; land = bucketu.bucket.σS, ocean = 0.0, ice = 0.0),
        )
        Plots.plot(combined_field)
    end
    Plots.mp4(anim, joinpath(out_dir, "bucket_snow.mp4"), fps = 10)

    if mode_name == "slabplanet_eisenman"
        sol_ice = cs.model_sims.ice_sim.integrator.sol
        combined_field = zeros(boundary_space)
        anim = Plots.@animate for sol_iceu in sol_ice.u
            Regridder.combine_surfaces_from_sol!(
                combined_field,
                surface_fractions,
                (; land = 0.0, ocean = 0.0, ice = sol_iceu.h_ice),
            )
            Plots.plot(combined_field)
        end
        Plots.mp4(anim, joinpath(out_dir, "eisenman_seaice.mp4"), fps = 10)
    end
    # plot surface fluxes
    # TODO as part of the flux accumulation PR
end
