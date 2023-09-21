"""
    Base.show(io::IO, dict::Dict)

This prints the keys and values of a Dict in sorted order.
"""
function Base.show(io::IO, dict::Dict)
    for k in sort!(collect(keys(dict)))
        println(io, " $k => $(dict[k])")
    end
end
