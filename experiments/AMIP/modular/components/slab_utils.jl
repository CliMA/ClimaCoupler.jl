"""
    SlabSimulation{P, Y, D, I}

The simulation structure for slab models.
"""
struct SlabSimulation{F, P, Y, D, I}
    FT::F
    params::P
    Y_init::Y
    domain::D
    integrator::I
end

"""
    get_slab_energy(slab_sim, T_sfc)

Returns the internal energy per unit area of the slab.
"""
get_slab_energy(slab_sim, T_sfc) =
    slab_sim.integrator.p.params.ρ .* slab_sim.integrator.p.params.c .* T_sfc .* slab_sim.integrator.p.params.h

get_slab_temperature(slab_sim, colidx) = slab_sim.integrator.u.T[colidx]
get_slab_humidity(slab_sim, colidx) = nothing
