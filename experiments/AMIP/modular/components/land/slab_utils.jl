"""
    SlabSimulation{P, Y, D, I}

The simulation structure for slab models.
"""
struct SlabSimulation{P, Y, D, I}
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
