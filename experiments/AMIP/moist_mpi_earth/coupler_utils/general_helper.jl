# most of these functions are temporary helpers until upstream issues are resolved

# TODO: unify with coupler interface
struct CouplerSimulation{C, I, F, S, D, B, T, P}
    comms_ctx::C
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
    monthly_3d_diags::NamedTuple
    monthly_2d_diags::NamedTuple
end

function swap_space!(field, new_space)
    field_out = zeros(new_space)
    parent(field_out) .= parent(field)
    return field_out
end

heaviside(var, FT) = var < FT(0) ? FT(0) : var
