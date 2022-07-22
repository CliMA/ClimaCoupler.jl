# most of these functions are temporary helpers until upstream issues are resolved

# TODO: unify with coupler interface
struct CouplerSimulation{I, F, D, B, T, M, P}
    Δt::I
    t::F
    dates::D
    boundary_space::B
    FT::T
    mask::M
    fields::NamedTuple
    model_sims::NamedTuple
    mode::NamedTuple
    parsed_args::P
end

function swap_space!(field, new_space)
    field_out = zeros(new_space)
    parent(field_out) .= parent(field)
    return field_out
end

heaviside(var) = var < 0 ? 0 : var
