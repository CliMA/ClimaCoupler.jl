# most of these functions are temporary helpers until upstream issues are resolved

# TODO: unify with coupler interface
struct CouplerSimulation{FT, B, M}
    boundary_space::B
    mask::M
    t::FT
    Δt::FT
end

CouplerSimulation{FT}(args...) where {FT} = CouplerSimulation{FT, typeof.(args[1:2])...}(args...)
float_type(::CouplerSimulation{FT}) where {FT} = FT

get_u(sim, t) = Geometry.UVVector.(sim.integrator.sol.u[t].c.uₕ).components.data.:1

function swap_space!(field_out, field_in)
    parent(field_out) .= parent(field_in)
    return field_out
end

heaviside(var) = var < 0 ? 0 : var
