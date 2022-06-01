# generalise this into coupler-specific function

abstract type AbstractCheck end

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
