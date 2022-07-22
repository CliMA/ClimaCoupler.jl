using ClimaCorePlots

function plot_anim(atmos_sim, slab_sim, slab_ocean_sim, slab_ice_sim, land_sea_mask, mode_name, SST)
    sol_atm = atmos_sim.integrator.sol

    anim = Plots.@animate for u in sol_atm.u
        Plots.plot(Fields.level(Geometry.UVVector.(u.c.uₕ).components.data.:1, 1))
    end
    Plots.mp4(anim, "anim_u.mp4", fps = 10)

    anim = Plots.@animate for u in sol_atm.u
        Plots.plot(Fields.level(Geometry.UVVector.(u.c.uₕ).components.data.:1, 5))
    end
    Plots.mp4(anim, "anim_u_7km.mp4", fps = 10)

    anim = Plots.@animate for u in sol_atm.u
        Plots.plot(Fields.level(u.c.ρe_tot, 1))
    end
    Plots.mp4(anim, "anim_rhoe.mp4", fps = 10)

    anim = Plots.@animate for u in sol_atm.u
        Plots.plot(Fields.level(u.c.ρe_tot, 1) .- Fields.level(sol_atm.u[1].c.ρe_tot, 1), clims = (-5000, 50000))
        println(parent(Fields.level(u.c.ρe_tot, 1) .- Fields.level(sol_atm.u[1].c.ρe_tot, 1))[1])
    end
    Plots.mp4(anim, "anim_rhoe_anom.mp4", fps = 10)

    anim = Plots.@animate for u in sol_atm.u
        Plots.plot(Fields.level(u.c.ρe_tot, 5) .- Fields.level(sol_atm.u[1].c.ρe_tot, 5), clims = (-1000, 3000))
    end
    Plots.mp4(anim, "anim_rhoe_anom_7km.mp4", fps = 10)

    anim = Plots.@animate for u in sol_atm.u
        Plots.plot(Fields.level(u.c.ρq_tot ./ u.c.ρ, 1))
    end
    Plots.mp4(anim, "anim_qt.mp4", fps = 10)

    anim = Plots.@animate for u in sol_atm.u
        Plots.plot(Fields.level(u.c.ρq_tot ./ u.c.ρ, 5))
    end
    Plots.mp4(anim, "anim_qt_7km.mp4", fps = 10)
    combined_field = zeros(boundary_space)
    sol_slab = slab_sim.integrator.sol
    FT = eltype(land_sea_mask)
    univ_mask = parent(land_sea_mask) .- parent(slab_ice_sim.integrator.p.Ya.ice_mask .* FT(2))

    if mode_name == "aquaplanet"
        
        sol_slab_ocean = slab_ocean_sim.integrator.sol
        T_ice = ice_sim.integrator.u.T_sfc

        anim = Plots.@animate for (bucketu, oceanu) in zip(sol_slab.u, sol_slab_ocean.u)
            parent(combined_field) .=
                combine_surface.(FT, univ_mask, parent(bucketu.bucket.T_sfc), parent(oceanu.T_sfc), parent(T_ice))
           

            Plots.plot(combined_field)
        end
    elseif mode_name == "amip"
         sol_slab_ice = slab_ice_sim.integrator.sol
        anim = Plots.@animate for (bucketu, iceu) in zip(sol_slab.u, sol_slab_ice.u)
            parent(combined_field) .=
                combine_surface.(FT, univ_mask, parent(bucketu.bucket.T_sfc), parent(SST), parent(iceu.T_sfc))
            Plots.plot(combined_field)
        end
    end

    Plots.mp4(anim, "earth_T.mp4", fps = 10)

    combined_field = zeros(boundary_space)
    anim = Plots.@animate for bucketu in sol_slab.u
        parent(combined_field) .= combine_surface.(FT, univ_mask, parent(bucketu.bucket.W), 0.0, 0.0)
        Plots.plot(combined_field)
    end
    Plots.mp4(anim, "bucket_W.mp4", fps = 10)

    combined_field = zeros(boundary_space)
    anim = Plots.@animate for bucketu in sol_slab.u
        parent(combined_field) .= combine_surface.(FT, univ_mask, parent(bucketu.bucket.σS), 0.0, 0.0)
        Plots.plot(combined_field)
    end
    Plots.mp4(anim, "bucket_σS.mp4", fps = 10)

end
