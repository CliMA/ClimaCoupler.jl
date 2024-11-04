import Plots
import ClimaCorePlots
import ClimaCore as CC
import ClimaCoupler: Regridder

function plot_anim(cs, out_dir = ".")

    atmos_sim = cs.model_sims.atmos_sim
    land_sim = cs.model_sims.land_sim
    ocean_sim = cs.model_sims.ocean_sim
    ice_sim = cs.model_sims.ice_sim
    mode_name = cs.mode.name

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
    FT = CC.Spaces.undertype(boundary_space)

    land_T_sfc = Interfacer.get_field(cs.model_sims.land_sim, Val(:surface_temperature))
    SST = cs.fields.T_S

    # Surface temperature plots
    if mode_name == "slabplanet"
        # Ocean surface temperature
        sol_ocean = ocean_sim.integrator.sol
        anim_T = Plots.@animate for oceanu in sol_ocean.u
            Regridder.combine_surfaces_from_sol!(
                combined_field,
                surface_fractions,
                (; land = land_T_sfc, ocean = oceanu.T_sfc, ice = FT(0)),
            )
            Plots.plot(combined_field)
        end
    elseif mode_name == "slabplanet_eisenman"
        # Ice surface temperature
        sol_ice = ice_sim.integrator.sol
        anim_T = Plots.@animate for iceu in sol_ice.u
            Regridder.combine_surfaces_from_sol!(
                combined_field,
                surface_fractions,
                (; land = land_T_sfc, ocean = FT(0), ice = iceu.T_sfc),
            )
            Plots.plot(combined_field)
        end

        # Ice height
        sol_ice = ice_sim.integrator.sol
        anim = Plots.@animate for sol_iceu in sol_ice.u
            Regridder.combine_surfaces_from_sol!(
                combined_field,
                surface_fractions,
                (; land = FT(0), ocean = FT(0), ice = sol_iceu.h_ice),
            )
            Plots.plot(combined_field)
        end
        Plots.mp4(anim, joinpath(out_dir, "eisenman_seaice.mp4"), fps = 10)
    elseif mode_name == "amip"
        # Ice surface temperature
        sol_ice = ice_sim.integrator.sol
        anim_T = Plots.@animate for iceu in sol_ice.u
            Regridder.combine_surfaces_from_sol!(
                combined_field,
                surface_fractions,
                (; land = land_T_sfc, ocean = SST, ice = iceu.T_sfc),
            )
            Plots.plot(combined_field)
        end
    end
    Plots.mp4(anim_T, joinpath(out_dir, "earth_T.mp4"), fps = 10)

    # Land surface plots
    sol_land = land_sim.integrator.sol

    # Water content
    anim = Plots.@animate for bucketu in sol_land.u
        Regridder.combine_surfaces_from_sol!(
            combined_field,
            surface_fractions,
            (; land = bucketu.bucket.W, ocean = FT(0), ice = FT(0)),
        )
        Plots.plot(combined_field)
    end
    Plots.mp4(anim, joinpath(out_dir, "bucket_W.mp4"), fps = 10)

    # Snow cover fraction
    anim = Plots.@animate for bucketu in sol_land.u
        Regridder.combine_surfaces_from_sol!(
            combined_field,
            surface_fractions,
            (; land = bucketu.bucket.σS, ocean = FT(0), ice = FT(0)),
        )
        Plots.plot(combined_field)
    end
    Plots.mp4(anim, joinpath(out_dir, "bucket_snow.mp4"), fps = 10)
end
