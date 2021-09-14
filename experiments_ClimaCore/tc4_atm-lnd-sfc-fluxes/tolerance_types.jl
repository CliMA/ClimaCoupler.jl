#### Tolerances

export ResidualTolerance, SolutionTolerance

# abstract type AbstractTolerance{FT <: AbstractFloat} end

"""
    ResidualTolerance

A tolerance type based on the residual of the equation ``f(x) = 0``
"""
struct ResidualTolerance{FT} <: AbstractTolerance{FT}
    tol::FT
end

"""
    (tol::ResidualTolerance)(x1, x2, y)

Evaluates residual tolerance, based on ``|f(x)|``
"""
(tol::ResidualTolerance)(x1, x2, y) = all(abs.(y) .< tol.tol)

"""
    SolutionTolerance

A tolerance type based on the solution ``x`` of the equation ``f(x) = 0``
"""
struct SolutionTolerance{FT} <: AbstractTolerance{FT}
    tol::FT
end

"""
    (tol::SolutionTolerance)(x1, x2, y)

Evaluates solution tolerance, based on ``|x2-x1|``
"""
(tol::SolutionTolerance)(x1, x2, y) = all(abs.(x2 .- x1) .< tol.tol)
