using ClimaCorePlots
function plot_anim(atmos_sim, slab_sim, slab_ocean_sim, slab_ice_sim) # TODO: uses global defs

    sol_atm = atmos_sim.integrator.sol
    sol_slab = slab_sim.integrator.sol
    sol_slab_ocean = slab_ocean_sim.integrator.sol
    sol_slab_ice = slab_ice_sim.integrator.sol

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

    anim = Plots.@animate for u in sol_slab.u
        Plots.plot(u.T_sfc)#,  clims = (240, 330))
    end
    Plots.mp4(anim, "slab_T.mp4", fps = 10)

    anim = Plots.@animate for u in sol_atm.u
        Plots.plot(Fields.level(u.c.ρe, 1) .- Fields.level(sol_atm.u[1].c.ρe, 1), clims = (-5000, 50000))
    end
    Plots.mp4(anim, "anim_rhoe_anom.mp4", fps = 10)

    anim = Plots.@animate for u in sol_atm.u
        Plots.plot(Fields.level(u.c.ρe, 5) .- Fields.level(sol_atm.u[1].c.ρe, 5), clims = (-1000, 3000))
    end
    Plots.mp4(anim, "anim_rhoe_anom_7km.mp4", fps = 10)

    anim = Plots.@animate for u in sol_atm.u
        Plots.plot(Fields.level(u.c.ρq_tot, 1))#.- Fields.level(sol_atm.u[1].c.ρt_tot,1),  clims = (-5000, 50000) )
    end
    Plots.mp4(anim, "anim_rhoqt.mp4", fps = 10)

    anim = Plots.@animate for u in sol_atm.u
        Plots.plot(Fields.level(u.c.ρq_tot, 2), clims = (0, 0.005))#.- Fields.level(sol_atm.u[1].c.ρt_tot,1),  clims = (-5000, 50000) )
    end
    Plots.mp4(anim, "anim_rhoqt_1km_v2.mp4", fps = 10)

    anim = Plots.@animate for u in sol_atm.u
        Plots.plot(Fields.level(u.c.ρq_tot, 5), clims = (0, 0.001))#.- Fields.level(sol_atm.u[1].c.ρt_tot,1),  clims = (-5000, 50000) )
    end
    Plots.mp4(anim, "anim_rhoqt_7km_v2.mp4", fps = 10)

    anim = Plots.@animate for u in sol_atm.u
        Plots.plot(mask)
    end
    Plots.mp4(anim, "mask.mp4", fps = 10)

    times = 0:saveat:t_end
    anim = Plots.@animate for t_i in 1:1:length(sol_slab.u)
        t = t_i / 24 / 60 / 60
        u = sol_slab.u[t_i]
        u_o = sol_slab_ocean.u[1]
        u_i = sol_slab_ice.u[t_i]
        combined_field = similar(u.T_sfc)
        #parent(combined_field) .= combine_surface.(parent(mask), parent(u.T_sfc), parent(SST) )
        parent(combined_field) .=
            combine_surface.(
                parent(mask) .- parent(slab_ice_sim.integrator.p.Ya.ice_mask .* FT(2)),
                parent(u.T_sfc),
                parent(u_o.T_sfc),
                parent(u_i.T_sfc),
            )
        Plots.plot(combined_field, clims = (265, 310), title = ("day: $t"))
    end
    Plots.mp4(anim, "slab_T_combo.mp4", fps = 10)

    if :h_ice in propertynames(sol_slab_ice.u[1])
        anim = Plots.@animate for t_i in 1:1:length(sol_slab.u)
            u = slab_ice_sim.integrator.sol.u[t_i]
            t = t_i / 24 / 60 / 60
            Plots.plot(
                u.h_ice .* swap_space!(abs.(mask .- FT(1)), axes(u.h_ice)),
                clims = (0, 0.35),
                title = ("day: $t"),
            )#.- Fields.level(sol_atm.u[1].c.ρt_tot,1),  clims = (-5000, 50000) )
        end
        Plots.mp4(anim, "h_ice.mp4", fps = 10)
    end

end
