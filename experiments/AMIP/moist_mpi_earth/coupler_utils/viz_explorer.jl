using ClimaCorePlots

function plot_anim(cs, out_dir = ".")

    atmos_sim = cs.model_sims.atmos_sim
    slab_land_sim = cs.model_sims.land_sim
    slab_ocean_sim = cs.model_sims.ocean_sim
    slab_ice_sim = cs.model_sims.ice_sim
    land_sea_mask = cs.land_mask
    mode_name = cs.mode.name
    SST = cs.fields.T_S

    # plot the lowest atmos (center) levels of key variables
    sol_atm = atmos_sim.integrator.sol

    anim = Plots.@animate for u in sol_atm.u
        Plots.plot(Fields.level(Geometry.UVVector.(u.c.uₕ).components.data.:1, 1))
    end
    Plots.mp4(anim, joinpath(out_dir, "anim_u.mp4"), fps = 10)

    anim = Plots.@animate for u in sol_atm.u
        Plots.plot(Fields.level(u.c.ρe_tot, 1) .- Fields.level(sol_atm.u[1].c.ρe_tot, 1), clims = (-5000, 50000))
        println(parent(Fields.level(u.c.ρe_tot, 1) .- Fields.level(sol_atm.u[1].c.ρe_tot, 1))[1])
    end
    Plots.mp4(anim, joinpath(out_dir, "anim_rhoe_anom.mp4"), fps = 10)

    anim = Plots.@animate for u in sol_atm.u
        Plots.plot(Fields.level(u.c.ρq_tot ./ u.c.ρ, 1))
    end
    Plots.mp4(anim, joinpath(out_dir, "anim_qt.mp4"), fps = 10)

    # plot combined surfaces
    combined_field = zeros(boundary_space)
    R_F = similar(combined_field)
    sol_slab = slab_land_sim.integrator.sol
    sol_slab_ice = slab_ice_sim.integrator.sol
    FT = eltype(land_sea_mask)
    univ_mask = parent(land_sea_mask) .- parent(slab_ice_sim.integrator.p.Ya.ice_mask .* FT(2))

    if mode_name == "slabplanet"
        sol_slab_ocean = slab_ocean_sim.integrator.sol

        anim = Plots.@animate for (bucketu, oceanu, iceu) in zip(sol_slab.u, sol_slab_ocean.u, sol_slab_ice.u)
            parent(combined_field) .=
                combine_surface.(FT, univ_mask, parent(bucketu.bucket.T_sfc), parent(oceanu.T_sfc), parent(iceu.T_sfc))

            Plots.plot(combined_field)
        end
    elseif mode_name == "amip"
        anim = Plots.@animate for (bucketu, iceu) in zip(sol_slab.u, sol_slab_ice.u)
            parent(combined_field) .=
                combine_surface.(FT, univ_mask, parent(bucketu.bucket.T_sfc), parent(SST), parent(iceu.T_sfc))

            Plots.plot(combined_field)
        end
    end
    Plots.mp4(anim, joinpath(out_dir, "earth_T.mp4"), fps = 10)

    combined_field = zeros(boundary_space)
    anim = Plots.@animate for bucketu in sol_slab.u
        parent(combined_field) .= combine_surface.(FT, univ_mask, parent(bucketu.bucket.W), 0.0, 0.0)

        Plots.plot(combined_field)
    end
    Plots.mp4(anim, joinpath(out_dir, "bucket_W.mp4"), fps = 10)

    # plot surface fluxes
    # TODO as part of the flux accumulation PR
end
