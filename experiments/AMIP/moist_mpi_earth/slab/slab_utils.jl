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
    get_slab_energy(slab_sim, boundary_space)

Returns the internal energy per unit area of the slab.
"""
function get_slab_energy(slab_sim::SlabSimulation, boundary_space)
    T_sfc = swap_space!(slab_sim.integrator.u.T_sfc, boundary_space)
    return slab_sim.params.ρ .* slab_sim.params.c .* T_sfc.* slab_sim.params.h
end
