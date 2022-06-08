# generalise this into coupler-specific function

abstract type AbstractCheck end

struct OnlineConservationCheck{A} <: AbstractCheck
    ρe_tot_atmos::A
    ρe_tot_land::A
    ρe_tot_ocean::A
    ρe_tot_seaice::A
end
function check_conservation(
    cs::OnlineConservationCheck,
    coupler_sim,
    atmos_sim = nothing,
    land_sim = nothing,
    ocean_sim = nothing,
    seaice_sim = nothing,
    radiation = true,
)

    times = (coupler_sim.Δt):(coupler_sim.Δt):(coupler_sim.t)
    z = parent(Fields.coordinate_field(face_space).z)
    Δz_bot = FT(0.5) * (z[2, 1, 1, 1, 1] - z[1, 1, 1, 1, 1])

    u_atm = atmos_sim !== nothing ? atmos_sim.integrator.u.c.ρe : nothing
    u_lnd = land_sim !== nothing ? swap_space!(land_sim.integrator.u.bucket.T_sfc, coupler_sim.boundary_space) : nothing
    u_ocn = ocean_sim !== nothing ? swap_space!(ocean_sim.integrator.u.T_sfc, coupler_sim.boundary_space) : nothing
    u_ice = seaice_sim !== nothing ? swap_space!(seaice_sim.integrator.u.T_sfc, coupler_sim.boundary_space) : nothing
    # global sums
    if atmos_sim !== nothing
        atmos_e = sum(u_atm)

        if radiation

            LWu_TOA = Fields.level(
                array2field(FT.(atmos_sim.integrator.p.rrtmgp_model.face_lw_flux_up), face_space),
                10 + half,
            )
            SWd_TOA = Fields.level(
                array2field(FT.(atmos_sim.integrator.p.rrtmgp_model.face_sw_flux_dn), face_space),
                10 + half,
            )
            SWu_TOA = Fields.level(
                array2field(FT.(atmos_sim.integrator.p.rrtmgp_model.face_sw_flux_up), face_space),
                10 + half,
            )

            Δz_top = FT(0.5) * (z[end, 1, 1, 1, 1] - z[end - 1, 1, 1, 1, 1])
            radiation_sources = -sum(SWd_TOA .- LWu_TOA .- SWu_TOA) ./ Δz_top
            radiation_sources_tsrs = radiation_sources .* atmos_sim.integrator.t # * coupler_sim.Δt # accumulate flux over coupling timestep (TODO: accumulate in atmos)

            atmos_e = atmos_e .+ radiation_sources_tsrs # J 
        end
    end
    univ_mask = coupler_sim.mask
    if (prescribed_sst != true) && (prescribed_sic == true)
        # zero out the ocean where there's sea ice
        # zero out the land where there is ocean & ice

        mask = coupler_sim.mask
        univ_mask = parent(mask) .- parent(slab_ice_sim.integrator.p.ice_mask .* FT(2))
        ice_mask(u_ice_1, univ_mask) = (univ_mask ≈ FT(-2) ? u_ice_1 : FT(0))
        parent(u_ice) .= ice_mask.(parent(u_ice), univ_mask)
    end
    if (prescribed_sst != true) 
        ocean_mask(u_ocn_1, univ_mask) = (univ_mask ≈ FT(0) ? u_ocn_1 : FT(0))
        parent(u_ocn) .= ocean_mask.(parent(u_ocn), univ_mask)
    end
    land_mask(u_lnd_1, univ_mask) = (univ_mask ≈ FT(1) ? u_lnd_1 : FT(0))
    parent(u_lnd) .= land_mask.(parent(u_lnd), univ_mask)


    land_e = land_sim !== nothing ? sum(get_bucket_energy(land_sim, u_lnd))./ Δz_bot : FT(0)
    ocean_e = ocean_sim !== nothing ? sum(get_slab_energy(ocean_sim, u_ocn))./ Δz_bot : FT(0)
    seaice_e = seaice_sim !== nothing ? sum(get_slab_energy(seaice_sim, u_ice))./ Δz_bot : FT(0)
    
    # save to coupler cache
    atmos_sim !== nothing ? push!(cs.ρe_tot_atmos, atmos_e) : nothing
    land_sim !== nothing ? push!(cs.ρe_tot_land, land_e) : nothing
    ocean_sim !== nothing ? push!(cs.ρe_tot_ocean, ocean_e) : nothing
    seaice_sim !== nothing ? push!(cs.ρe_tot_seaice, seaice_e) : nothing

end


struct OfflineConservationCheck <: AbstractCheck end
function check_conservation(
    ::OfflineConservationCheck,
    coupler_sim,
    atmos_sim = nothing,
    land_sim = nothing,
    ocean_sim = nothing,
    seaice_sim = nothing,
    radiation = true,
    figname = "test.png",
)

    times = (0.0):(coupler_sim.Δt):(coupler_sim.t)
    z = parent(Fields.coordinate_field(face_space).z)
    Δz_bot = FT(0.5) * (z[2, 1, 1, 1, 1] - z[1, 1, 1, 1, 1])

    solu_atm = atmos_sim !== nothing ? atmos_sim.integrator.sol.u : nothing
    solu_lnd =
        land_sim !== nothing ?
        Fields.FieldVector(
            T_sfc = [swap_space!(u.bucket.T_sfc, coupler_sim.boundary_space) for u in land_sim.integrator.sol.u],
        ) : nothing
    solu_ocn =
        ocean_sim !== nothing ?
        Fields.FieldVector(
            T_sfc = [swap_space!(u.T_sfc, coupler_sim.boundary_space) for u in ocean_sim.integrator.sol.u],
        ) : nothing
    solu_ice =
        seaice_sim !== nothing ?
        Fields.FieldVector(
            T_sfc = [swap_space!(u.T_sfc, coupler_sim.boundary_space) for u in seaice_sim.integrator.sol.u],
        ) : nothing

    # global sums
    if atmos_sim !== nothing

        atmos_e = [sum(u.c.ρe) for u in solu_atm]

        if radiation

            LWu_TOA = Fields.level(
                array2field(FT.(atmos_sim.integrator.p.rrtmgp_model.face_lw_flux_up), face_space),
                10 + half,
            )
            SWd_TOA = Fields.level(
                array2field(FT.(atmos_sim.integrator.p.rrtmgp_model.face_sw_flux_dn), face_space),
                10 + half,
            )
            SWu_TOA = Fields.level(
                array2field(FT.(atmos_sim.integrator.p.rrtmgp_model.face_sw_flux_up), face_space),
                10 + half,
            )

            Δz_top = FT(0.5) * (z[end, 1, 1, 1, 1] - z[end - 1, 1, 1, 1, 1])
            radiation_sources = -sum(SWd_TOA .- LWu_TOA .- SWu_TOA) ./ Δz_top
            radiation_sources_tsrs = radiation_sources * times

            atmos_e = atmos_e .+ radiation_sources_tsrs # J 
        end
    end

    land_e = land_sim !== nothing ? [sum(get_bucket_energy(land_sim, u)) for u in solu_lnd] ./ Δz_bot : [FT(0)]
    ocean_e = ocean_sim !== nothing ? [sum(get_slab_energy(ocean_sim, u)) for u in solu_ocn] ./ Δz_bot : [FT(0)]
    seaice_e = seaice_sim !== nothing ? [sum(get_slab_energy(seaice_sim, u)) for u in solu_ice] ./ Δz_bot : [FT(0)]

    diff_ρe_tot_atmos = atmos_e .- atmos_e[1]
    diff_ρe_tot_slab = (land_e .- land_e[1])
    diff_ρe_tot_slab_ocean = (ocean_e .- ocean_e[1])
    diff_ρe_tot_slab_seaice = (seaice_e .- seaice_e[1])

    Plots.plot(diff_ρe_tot_atmos, label = "atmos")
    Plots.plot!(diff_ρe_tot_slab, label = "land")
    Plots.plot!(diff_ρe_tot_slab_ocean, label = "ocean")
    Plots.plot!(diff_ρe_tot_slab_seaice, label = "seaice")
    tot = atmos_e .+ ocean_e .+ land_e .+ seaice_e
    times_days = floor.(times ./ (24 * 60 * 60))
    Plots.plot!(
        tot .- tot[1],
        label = "tot",
        xlabel = "time [days]",
        ylabel = "energy(t) - energy(t=0) [J]",
        xticks = (collect(1:length(times))[1:50:end], times_days[1:50:end]),
    )
    Plots.savefig(figname)

end
#=
struct ConservationCheck{A} <: AbstractCheck
    ρe_tot_atmos::A
    ρe_tot_bucket::A
    dE_expected::A
    dW_expected::A
    water_tot_bucket::A
end

function check_conservation_callback(cs, atmos_sim, bucket_sim, energy_flux, water_flux)

    Δz_1 =
        (parent(ClimaCore.Fields.coordinate_field(atmos_sim.domain.face_space).z)[2] -
        parent(ClimaCore.Fields.coordinate_field(atmos_sim.domain.face_space).z)[1]) ./ FT(2.0)
    dt = bucket_sim.integrator.dt
    atmos_field = atmos_sim.integrator.u.c.ρe  # J 
    bucket_field = get_bucket_energy(bucket_sim, bucket_sim.integrator.u.bucket.T_sfc)  
    #[NB: sum of the boundary field inherits the depth from the first atmospheric layer, which ≂̸ bucket depth]
    ρe_tot_atmos = sum(atmos_field) # ∫ ρe dV
    ρe_tot_bucket = sum(bucket_field)./ Δz_1 # ∫ ρc T*d_soil dz / Δz_1 
    dE_expected = sum(energy_flux .* dt ./ Δz_1) 
    dW_expected = sum(water_flux .* dt ./ Δz_1) 
    water_tot_bucket = sum(bucket_sim.integrator.u.bucket.W .+bucket_sim.integrator.u.bucket.Ws) ./ Δz_1
    push!(cs.ρe_tot_atmos, ρe_tot_atmos)
    push!(cs.ρe_tot_bucket, ρe_tot_bucket)
    push!(cs.dE_expected, dE_expected)
    push!(cs.water_tot_bucket, water_tot_bucket)
    push!(cs.dW_expected, dW_expected)
end

using ClimaCorePlots
function conservation_plot(CS::ConservationCheck, t, figname)
    diff_ρe_tot_atmos = CS.ρe_tot_atmos[2:end] .- CS.ρe_tot_atmos[1] # Eatmos - Eatmos(0), negative = energy leaving atmos to bucket
    diff_ρe_tot_bucket = CS.ρe_tot_bucket[2:end] .- CS.ρe_tot_bucket[1] # Eland - Eland(0), positive = energy entering from atmos
    Plots.plot(t[2:end],diff_ρe_tot_atmos, label = "atmos")
    Plots.plot!(t[2:end],diff_ρe_tot_bucket, label = "bucket") #tot = Eland+Eatmos - Eland(0)- Eatmos(0) = E_earth - E_earth(0)
    tot = diff_ρe_tot_atmos .+ diff_ρe_tot_bucket
    dE = CS.dE_expected
    Plots.plot!(t[2:end],tot, label = "tot")
    Plots.plot!(t[2:end],cumsum(dE[1:end-1]), label= "Cumulative sum of ∫F*dA*dt")
    diff_by_step_atmos = CS.ρe_tot_atmos[2:end] - CS.ρe_tot_atmos[1:end-1]
    diff_by_step_bucket = CS.ρe_tot_bucket[2:end] - CS.ρe_tot_bucket[1:end-1]
    Plots.savefig(string(figname,"_energy.png"))

    Plots.plot(t[2:end], abs.(diff_by_step_atmos .- dE[1:end-1]), yaxis = :log, label = "d E_atmos - dE step")
    Plots.plot!(t[2:end], abs.(diff_by_step_bucket.+ dE[1:end-1]), yaxis = :log, label = "d E_bucket+ dE_step")
    Plots.plot!(ylabel = "Change in energy in a step (J)")
    Plots.plot!(xlabel = "time(s)")
    diff_by_step_bucket_water = CS.water_tot_bucket[2:end] - CS.water_tot_bucket[1:end-1]
    Plots.savefig(string(figname,"_denergy.png"))
    dW = CS.dW_expected
    Plots.plot(t[2:end], abs.(diff_by_step_bucket_water.+ dW[1:end-1]), yaxis = :log, label = "d W_bucket+ dW_step")
    Plots.plot!(ylabel = "Change in water in a step (m)")
    Plots.plot!(xlabel = "time(s)")
    Plots.savefig(string(figname,"_water.png"))
    
end

=#
#=
# for land-sea-atmos
times = 0:saveat:t_end
solu_atm = sol_atm.u
h_space = make_horizontal_space(horizontal_mesh, quad, nothing) #TODO move this to the beginning (once same the instance error sorted)
solu_slab = Fields.FieldVector(T_sfc = [Fields.Field(Fields.field_values(u.bucket.T_sfc), h_space) for u in sol_slab.u])
solu_slab_ocean = Fields.FieldVector(T_sfc = [Fields.Field(Fields.field_values(u.T_sfc), h_space) for u in sol_slab_ocean.u])


atmos_e = [sum(u.c.ρe) for u in solu_atm] # J 
z = parent(ClimaCore.Fields.coordinate_field(atmos_sim.domain.face_space).z)
Δz_1 = z[2] - z[1]
slab_e = [sum(get_bucket_energy(bucket_sim, u)) for u in solu_slab] 
slab_ocean_e = [sum(get_slab_energy(slab_ocean_sim, u)) for u in solu_slab_ocean] 



diff_ρe_tot_atmos = atmos_e .- atmos_e[3]
diff_ρe_tot_slab = (slab_e .- slab_e[3])
diff_ρe_tot_slab_ocean = (slab_ocean_e .- slab_ocean_e[3])

Plots.plot(diff_ρe_tot_atmos, label = "atmos")
Plots.plot!(diff_ρe_tot_slab, label = "slab")
Plots.plot!(diff_ρe_tot_slab_ocean, label = "slab_ocean")
tot = atmos_e .+ slab_ocean_e .+ slab_e
times_days = floor.(times ./ (24*60*60))
Plots.plot!(tot .- tot[1], label = "tot", xlabel = "time [days]", ylabel = "energy(t) - energy(t=0) [J]", xticks = ( collect(1:length(times))[1:50:end], times_days[1:50:end]) )
Plots.savefig(figname)
=#
