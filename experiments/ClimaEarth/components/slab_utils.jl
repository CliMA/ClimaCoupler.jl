"""
      weighted_dss_slab!(Y::CC.Fields.FieldVector, p::NamedTuple, _)

Computes the weighted direct stiffness summation and updates `Y` in place.
In the case of a column domain, no dss operations are performed.
"""
function weighted_dss_slab!(Y::CC.Fields.FieldVector, p::NamedTuple, _)
    for key in propertynames(Y)
        field = getproperty(Y, key)
        buffer = get_dss_buffer(axes(field), p)
        CC.Spaces.weighted_dss!(field, buffer)
    end
end

get_dss_buffer(::CC.Spaces.SpectralElementSpace2D, p) = p.dss_buffer
