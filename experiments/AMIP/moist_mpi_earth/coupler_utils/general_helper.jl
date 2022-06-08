struct CouplerSimulation{I, F, B, T, M}
    Δt::I
    t::F
    boundary_space::B
    FT::T
    mask::M
end

get_u(sim, t) = Geometry.UVVector.(sim.integrator.sol.u[t].c.uₕ).components.data.:1

function swap_space!(field, new_space)
    field_out = zeros(new_space)
    parent(field_out) .= parent(field)
    return field_out
end

heaviside(var) = var < 0 ? 0 : var
