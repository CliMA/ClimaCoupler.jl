# most of these functions are temporary helpers until upstream issues are resolved

get_u(sim, t) = Geometry.UVVector.(sim.integrator.sol.u[t].c.uₕ).components.data.:1

function swap_space!(field, new_space)
    field_out = zeros(new_space)
    parent(field_out) .= parent(field)
    return field_out
end

heaviside(var) = var < 0 ? 0 : var

# Plots.plot(get_u(atmos_sim, 20) .- swap_space!(get_u(atmos_sim_old, 20), axes(get_u(atmos_sim, 20))) )
