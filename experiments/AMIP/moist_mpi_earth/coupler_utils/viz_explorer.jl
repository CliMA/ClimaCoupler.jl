function plot_anim() # TODO: uses global defs

    using ClimaCorePlots

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
        println(parent(Fields.level(u.c.ρe, 1) .- Fields.level(sol_atm.u[1].c.ρe, 1))[1])
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
        Plots.plot(mask)#.- Fields.level(sol_atm.u[1].c.ρt_tot,1),  clims = (-5000, 50000) )
    end
    Plots.mp4(anim, "mask.mp4", fps = 10)

    times = 0:saveat:t_end
    anim = Plots.@animate for t_i in 1:1:length(sol_slab.u)
        t = t_i / 24 / 60 / 60
        u = sol_slab.u[t_i]
        u_i = sol_slab_ice.u[t_i]
        combined_field = similar(u.T_sfc)
        #parent(combined_field) .= combine_surface.(parent(mask), parent(u.T_sfc), parent(SST) )
        parent(combined_field) .=
            combine_surface.(
                parent(mask) .- parent(slab_ice_sim.integrator.p.ice_mask .* FT(2)),
                parent(u.T_sfc),
                parent(SST),
                parent(u_i.T_sfc),
            )
        Plots.plot(combined_field, clims = (265, 310), title = ("day: $t"))
    end
    Plots.mp4(anim, "slab_T_combo.mp4", fps = 10)

    u1 = sol_slab.u[1]
    u_i1 = sol_slab_ice.u[1]
    combined_field1 = similar(sol_slab.u[1].T_sfc)
    parent(combined_field1) .=
        combine_surface.(
            parent(mask) .- parent(slab_ice_sim.integrator.p.ice_mask .* FT(2)),
            parent(u1.T_sfc),
            parent(SST),
            parent(u_i1.T_sfc),
        )
    anim = Plots.@animate for t_i in 1:1:length(sol_slab.u)
        t = t_i / 24 / 60 / 60
        u = sol_slab.u[t_i]
        u_i = sol_slab_ice.u[t_i]
        combined_field = similar(u.T_sfc)
        parent(combined_field) .=
            combine_surface.(
                parent(mask) .- parent(slab_ice_sim.integrator.p.ice_mask .* FT(2)),
                parent(u.T_sfc),
                parent(SST),
                parent(u_i.T_sfc),
            )
        Plots.plot(combined_field .- swap_space!(combined_field1, axes(combined_field)), title = ("day: $t"))
    end
    Plots.mp4(anim, "slab_T_combo_anom.mp4", fps = 10)
end
