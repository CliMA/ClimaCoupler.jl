# generalise this into coupler-specific function

abstract type AbstractCheck end

struct OnlineConservationCheck{A} <: AbstractCheck
    ρe_tot_atmos::A
    ρe_tot_land::A
    ρe_tot_ocean::A
    ρe_tot_seaice::A
end
function check_conservation(cs::OnlineConservationCheck, coupler_sim, atmos_sim = nothing, land_sim = nothing, ocean_sim = nothing, seaice_sim = nothing, radiation = true)

    u_atm = atmos_sim !== nothing ? atmos_sim.integrator.u.c.ρe : nothing
    u_lnd = land_sim  !== nothing ? swap_space!(land_sim.integrator.u.T_sfc, coupler_sim.boundary_space) : nothing
    u_ocn = ocean_sim !== nothing ? swap_space!(ocean_sim.integrator.u.T_sfc, coupler_sim.boundary_space) : nothing
    u_ice = seaice_sim !== nothing ? swap_space!(seaice_sim.integrator.u.T_sfc, coupler_sim.boundary_space) : nothing

    # global sums
    if atmos_sim !== nothing
        atmos_e = sum(u_atm)

        if radiation
            
            z = parent(Fields.coordinate_field(face_space).z)
            Δz_top = FT(0.5) * (z[end,1,1,1,1] - z[end-1,1,1,1,1])
            nz = length(z[:,1,1,1,1])

            LWu_TOA = Fields.level(array2field(FT.(atmos_sim.integrator.p.rrtmgp_model.face_lw_flux_up), face_space), nz-half)
            SWd_TOA = Fields.level(array2field(FT.(atmos_sim.integrator.p.rrtmgp_model.face_sw_flux_dn), face_space), nz-half)
            SWu_TOA = Fields.level(array2field(FT.(atmos_sim.integrator.p.rrtmgp_model.face_sw_flux_up), face_space), nz-half)

            radiation_sources = - sum(SWd_TOA .- LWu_TOA .- SWu_TOA) ./ Δz_top
            radiation_sources_tsrs = radiation_sources .* atmos_sim.integrator.t # * coupler_sim.Δt # accumulate flux over coupling timestep (TODO: accumulate in atmos)

            atmos_e = atmos_e .+ radiation_sources_tsrs # J 
        end
    end

    if (prescribed_sst !== true) && (prescribed_sic == true)
        # zero out the ocean where there's sea ice
        mask = coupler_sim.mask
        univ_mask = parent(mask) .- parent(slab_ice_sim.integrator.p.ice_mask .* FT(2))
        value1, value2 = (FT(-0.5), FT(0.5))
        ocean_mask(u_ocn_1, univ_mask) = ((univ_mask >= FT(value1) && (univ_mask <= FT(value2))) ? u_ocn_1 : FT(0)) 
        u_ocn_ = similar(u_ocn)
        parent(u_ocn_) .=     ocean_mask.(parent(u_ocn), parent(univ_mask))
        u_ocn = u_ocn_
    end

    land_e = land_sim !== nothing ? sum(get_slab_energy(land_sim, u_lnd)) : FT(0)
    ocean_e = ocean_sim !== nothing ? sum(get_slab_energy(ocean_sim, u_ocn)) : FT(0)
    seaice_e = seaice_sim !== nothing ? sum(get_slab_energy(seaice_sim, u_ice)) : FT(0)

    # save to coupler cache
    atmos_sim !== nothing ? push!(cs.ρe_tot_atmos, atmos_e) : nothing 
    land_sim !== nothing ? push!(cs.ρe_tot_land, land_e) : nothing 
    ocean_sim !== nothing ? push!(cs.ρe_tot_ocean, ocean_e) : nothing 
    seaice_sim !== nothing ? push!(cs.ρe_tot_seaice, seaice_e) : nothing 
    
end

struct OfflineConservationCheck{A} <: AbstractCheck end
function check_conservation(::OfflineConservationCheck, coupler_sim, atmos_sim = nothing, land_sim = nothing, ocean_sim = nothing, seaice_sim = nothing, radiation = true, figname = "test.png")

    times = coupler_sim.Δt:coupler_sim.Δt:coupler_sim.t
    z = parent(Fields.coordinate_field(face_space).z)
    Δz_bot = FT(0.5) * (z[2,1,1,1,1] - z[1,1,1,1,1])

    solu_atm = atmos_sim !== nothing ? atmos_sim.integrator.sol.u : nothing
    solu_lnd = land_sim !== nothing ? Fields.FieldVector(T_sfc = [swap_space!(u.T_sfc, coupler_sim.boundary_space) for u in land_sim.integrator.sol.u]) : nothing
    solu_ocn = ocean_sim !== nothing ? Fields.FieldVector(T_sfc = [swap_space!(u.T_sfc, coupler_sim.boundary_space) for u in ocean_sim.integrator.sol.u]) : nothing    
    solu_ice = seaice_sim !== nothing ? Fields.FieldVector(T_sfc = [swap_space!(u.T_sfc, coupler_sim.boundary_space) for u in seaice_sim.integrator.sol.u]) : nothing    
    
    # global sums
    if atmos_sim !== nothing
        
        atmos_e = [sum(u.c.ρe) for u in solu_atm]

        if radiation
            
            LWu_TOA = Fields.level(array2field(FT.(atmos_sim.integrator.p.rrtmgp_model.face_lw_flux_up), face_space), 10+half)
            SWd_TOA = Fields.level(array2field(FT.(atmos_sim.integrator.p.rrtmgp_model.face_sw_flux_dn), face_space), 10+half)
            SWu_TOA = Fields.level(array2field(FT.(atmos_sim.integrator.p.rrtmgp_model.face_sw_flux_up), face_space), 10+half)

            Δz_top = FT(0.5) * (z[end,1,1,1,1] - z[end-1,1,1,1,1])
            radiation_sources = - sum(SWd_TOA .- LWu_TOA .- SWu_TOA) ./ Δz_top
            radiation_sources_tsrs = radiation_sources * times

            atmos_e = atmos_e .+ radiation_sources_tsrs # J 
        end
    end

    land_e = land_sim !== nothing ? [sum(get_slab_energy(land_sim, u)) for u in solu_lnd]  ./ Δz_bot : [FT(0)]
    ocean_e = ocean_sim !== nothing ? [sum(get_slab_energy(ocean_sim, u)) for u in solu_ocn]  ./ Δz_bot : [FT(0)]
    seaice_e = seaice_sim !== nothing ? [sum(get_slab_energy(seaice_sim, u)) for u in solu_ice] ./ Δz_bot : [FT(0)]

    diff_ρe_tot_atmos = atmos_e .- atmos_e[1]
    diff_ρe_tot_slab = (land_e .- land_e[1]) 
    diff_ρe_tot_slab_ocean = (ocean_e .- ocean_e[1]) 
    diff_ρe_tot_slab_seaice = (seaice_e .- seaice_e[1]) 
    
    Plots.plot(diff_ρe_tot_atmos, label = "atmos")
    Plots.plot!(diff_ρe_tot_slab, label = "land")
    Plots.plot!(diff_ρe_tot_slab_ocean, label = "ocean")
    Plots.plot!(diff_ρe_tot_slab_ocean, label = "seaice")
    tot = atmos_e .+ ocean_e .+ land_e .+ seaice_e
    times_days = floor.(times ./ (24*60*60))
    Plots.plot!(tot .- tot[1], label = "tot", xlabel = "time [days]", ylabel = "energy(t) - energy(t=0) [J]", xticks = ( collect(1:length(times))[1:50:end], times_days[1:50:end]) )
    Plots.savefig(figname)

end

#=
times = coupler_sim.Δt:coupler_sim.Δt:coupler_sim.t
figname = "Z_lets_see"
diff_ρe_tot_atmos = (CS.ρe_tot_atmos .- CS.ρe_tot_atmos[1])
diff_ρe_tot_slab = (CS.ρe_tot_land .- CS.ρe_tot_land[1]) ./ 250
diff_ρe_tot_slab_seaice = (CS.ρe_tot_seaice .- CS.ρe_tot_seaice[1]) ./ 250
diff_ρe_tot_slab_ocean = (CS.ρe_tot_ocean .- CS.ρe_tot_ocean[1]) ./ 250

Plots.plot(diff_ρe_tot_atmos, label = "atmos")
Plots.plot!(diff_ρe_tot_slab, label = "land")
Plots.plot!(diff_ρe_tot_slab_ocean, label = "ocean")
Plots.plot!(diff_ρe_tot_slab_seaice, label = "seaice")
tot = CS.ρe_tot_atmos  .+ (CS.ρe_tot_ocean .+ CS.ρe_tot_land .+ CS.ρe_tot_seaice) ./ 250
times_days = floor.(times ./ (24*60*60))
Plots.plot!(tot .- tot[1], label = "tot", xlabel = "time [days]", ylabel = "energy(t) - energy(t=0) [J]", xticks = ( collect(1:length(times))[1:50:end], times_days[1:50:end]) )


Plots.savefig(figname)
=#

# # NB
# sp= axes(SWu_TOA)
# 4π*6.371229f6^2 == sum(ones(sp)) / 2.5e3
