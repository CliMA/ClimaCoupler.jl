# most of these functions are temporary helpers until upstream issues are resolved

# TODO: unify with coupler interface
struct CouplerSimulation{I, F, S, D, B, T, P}
    Δt_cpl::I
    t::F
    tspan::S
    dates::D
    boundary_space::B
    FT::T
    surface_masks::NamedTuple
    fields::NamedTuple
    model_sims::NamedTuple
    mode::NamedTuple
    parsed_args::P
    monthly_state_diags::NamedTuple
end

function swap_space!(field, new_space)
    field_out = zeros(new_space)
    parent(field_out) .= parent(field)
    return field_out
end

heaviside(var) = var < 0 ? 0 : var
