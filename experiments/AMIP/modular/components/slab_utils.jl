"""
    get_slab_energy(slab_sim, T_sfc)

Returns the internal energy per unit area of the slab.
"""
get_slab_energy(slab_sim, T_sfc) =
    slab_sim.integrator.p.params.ρ .* slab_sim.integrator.p.params.c .* T_sfc .* slab_sim.integrator.p.params.h

"""
      weighted_dss_slab!(Y::ClimaCore.Fields.FieldVector, p::NamedTuple, _)

Computes the weighted direct stiffness summation and updates `Y` in place.
In the case of a column domain, no dss operations are performed.
"""
function weighted_dss_slab!(Y::ClimaCore.Fields.FieldVector, p::NamedTuple, _)
    for key in propertynames(Y)
        field = getproperty(Y, key)
        buffer = get_dss_buffer(axes(field), p)
        ClimaCore.Spaces.weighted_dss!(field, buffer)
    end
end

get_dss_buffer(::ClimaCore.Spaces.SpectralElementSpace2D, p) = p.dss_buffer
