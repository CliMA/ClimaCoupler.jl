# Common to slab land, slab ocean, slab ice
struct SlabSimulation{P, Y, D, I}
    params::P
    Y_init::Y
    domain::D
    integrator::I
end
