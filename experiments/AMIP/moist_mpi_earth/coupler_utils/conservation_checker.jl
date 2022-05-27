# generalise this into coupler-specific function

abstract type AbstractCheck end

struct ConservationCheck{A} <: AbstractCheck
    ρe_tot_atmos::A
    ρe_tot_bucket::A
end

function check_conservation_callback(cs, atmos_sim, bucket_sim)

    atmos_sim.domain.center_space.vertical_topology.mesh.domain
    Δz_1 =
        parent(ClimaCore.Fields.coordinate_field(atmos_sim.domain.face_space).z)[2] -
        parent(ClimaCore.Fields.coordinate_field(atmos_sim.domain.face_space).z)[1]

    atmos_field = atmos_sim.integrator.u.c.ρe  # J 
    bucket_field = get_bucket_energy(bucket_sim, bucket_sim.integrator.u.bucket.T_sfc) ./ Δz_1  # J  [NB: sum of the boundary field inherits the depth from the first atmospheric layer, which ≂̸ bucket depth]

    ρe_tot_atmos = sum(atmos_field)
    ρe_tot_bucket = sum(bucket_field)

    push!(cs.ρe_tot_atmos, ρe_tot_atmos)
    push!(cs.ρe_tot_bucket, ρe_tot_bucket)
end

using ClimaCorePlots
function plot(CS::ConservationCheck)
    diff_ρe_tot_atmos = CS.ρe_tot_atmos .- CS.ρe_tot_atmos[3]
    diff_ρe_tot_bucket = (CS.ρe_tot_bucket .- CS.ρe_tot_bucket[3])
    Plots.plot(diff_ρe_tot_atmos, label = "atmos")
    Plots.plot!(diff_ρe_tot_bucket, label = "bucket")
    tot = diff_ρe_tot_atmos .+ diff_ρe_tot_bucket
    Plots.plot!(tot .- tot[1], label = "tot")
end

function conservation_plot(atmos_sim, bucket_sim, solu_atm, solu_bucket, figname = "tst_c.png")
    z = parent(ClimaCore.Fields.coordinate_field(atmos_sim.domain.face_space).z);
    Δz_1 = (z[2] - z[1])./2.0
    
    atmos_e = [sum(u.c.ρe) for u in solu_atm] # J 


    bucket_e = [sum(get_bucket_energy(bucket_sim, u)) for u in solu_bucket]  # J  [NB: sum of the boundary field inherits the depth from the first atmospheric layer, which ≂̸ bucket depth]

    diff_ρe_tot_atmos = atmos_e .- atmos_e[3]
    diff_ρe_tot_bucket = (bucket_e .- bucket_e[3])
    Plots.plot(diff_ρe_tot_atmos, label = "atmos")
    Plots.plot!(diff_ρe_tot_bucket, label = "bucket")
    tot = diff_ρe_tot_atmos .+ diff_ρe_tot_bucket
    Plots.plot!(tot .- tot[1], label = "tot", xlabel = "time [s]", ylabel = "energy(t) - energy(t=0) [s]")
    Plots.savefig(figname)

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
