using ClimaCorePlots

function plot_anim() # TODO: uses global defs



    anim = Plots.@animate for u in sol_atm.u
        Plots.plot(Fields.level(Geometry.UVVector.(u.c.uₕ).components.data.:1, 1))
    end
    Plots.mp4(anim, "anim_u.mp4", fps = 10)

    anim = Plots.@animate for u in sol_atm.u
        Plots.plot(Fields.level(Geometry.UVVector.(u.c.uₕ).components.data.:1, 5))
    end
    Plots.mp4(anim, "anim_u_7km.mp4", fps = 10)

    anim = Plots.@animate for u in sol_atm.u
        Plots.plot(Fields.level(u.c.ρe, 1))
    end
    Plots.mp4(anim, "anim_rhoe.mp4", fps = 10)
    combined_field = zeros(boundary_space)
    
    anim = Plots.@animate for (bucketu,oceanu, iceu) in zip(sol_slab.u,sol_slab_ocean.u, sol_slab_ice.u)
        parent(combined_field) .=
            combine_surface.(
                parent(mask) .- parent(slab_ice_sim.integrator.p.ice_mask .* FT(2)),
                parent(bucketu.bucket.T_sfc),
                parent(oceanu.T_sfc),
                parent(iceu.T_sfc),
            )
        dummmy_remap!(T_S, combined_field)
        
        Plots.plot(T_S, clims = (240, 330))
    end
    Plots.mp4(anim, "earth_T.mp4", fps = 10)

    anim = Plots.@animate for u in sol_atm.u
        Plots.plot(Fields.level(u.c.ρe, 1) .- Fields.level(sol_atm.u[1].c.ρe, 1), clims = (-5000, 50000))
        println(parent(Fields.level(u.c.ρe, 1) .- Fields.level(sol_atm.u[1].c.ρe, 1))[1])
    end
    Plots.mp4(anim, "anim_rhoe_anom.mp4", fps = 10)

    anim = Plots.@animate for u in sol_atm.u
        Plots.plot(Fields.level(u.c.ρe, 5) .- Fields.level(sol_atm.u[1].c.ρe, 5), clims = (-1000, 3000))
    end
    Plots.mp4(anim, "anim_rhoe_anom_7km.mp4", fps = 10)

    anim = Plots.@animate for u in sol_atm.u
        Plots.plot(Fields.level(u.c.ρq_tot, 1))
    end
    Plots.mp4(anim, "anim_rhoqt.mp4", fps = 10)

    anim = Plots.@animate for u in sol_atm.u
        Plots.plot(Fields.level(u.c.ρq_tot, 2), clims = (0, 0.005))
    end
    Plots.mp4(anim, "anim_rhoqt_1km_v2.mp4", fps = 10)

    anim = Plots.@animate for u in sol_atm.u
        Plots.plot(Fields.level(u.c.ρq_tot, 5), clims = (0, 0.001))
    end
    Plots.mp4(anim, "anim_rhoqt_7km_v2.mp4", fps = 10)
end
