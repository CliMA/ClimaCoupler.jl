using Literate
Literate.markdown("foo.jl", "."; documenter=true)
Literate.notebook("foo.jl", "."; documenter=true)
