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

Returns the volumetric internal energy of the slab.
"""
get_slab_energy(slab_sim, T_sfc) = slab_sim.params.ρ .* slab_sim.params.c .* T_sfc .* slab_sim.params.h
