function Base.show(io::IO, dict::Dict)
    for k in sort!(collect(keys(dict)))
        println(io, " $k => $(dict[k])")
    end
 end
