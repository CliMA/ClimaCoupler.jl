# generalise this into coupler-specific function

struct ConservationCheck{A}
    ρe_tot_atmos::A
    ρe_tot_slab::A
end

function check_conservation(cs, atmos_sim, slab_sim)
    atmos_field = atmos_sim.integrator.u.thermodynamics.ρe_tot
    slab_field = get_total_energy(slab_sim)

    ρe_tot_atmos = sum(atmos_field)
    ρe_tot_slab = sum(slab_field)

    push!(cs.ρe_tot_atmos, ρe_tot_atmos)
    push!(cs.ρe_tot_slab, ρe_tot_slab)
end

function get_total_energy(slab_sim)
    ρe_tot = slab_sim.params.ρ .* slab_sim.params.c .* slab_sim.integrator.u.T_sfc .* slab_sim.params.h
end

# struct ConservationCheck{F, E, I, S}
#     fields::F
#     exception::E
#     interval::I
#     show::S
# end
