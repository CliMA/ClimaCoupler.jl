# most of these functions are temporary helpers until upstream issues are resolved

# TODO: unify with coupler interface
struct CouplerSimulation{FT, C, S, D, B, P}
    comms_ctx::C
    tspan::S
    dates::D
    boundary_space::B
    parsed_args::P
    t::FT
    Δt_cpl::FT
    surface_masks::NamedTuple
    fields::NamedTuple
    model_sims::NamedTuple
    mode::NamedTuple
    monthly_3d_diags::NamedTuple
    monthly_2d_diags::NamedTuple
end

CouplerSimulation{FT}(args...) where {FT} = CouplerSimulation{FT, typeof.(args[1:5])...}(args...)
float_type_cs(::CouplerSimulation{FT}) where {FT} = FT

function swap_space!(field, new_space)
    field_out = zeros(new_space)
    parent(field_out) .= parent(field)
    return field_out
end

heaviside(var, FT) = var < FT(0) ? FT(0) : var
