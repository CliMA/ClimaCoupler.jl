#### Solution types

export CompactSolution, VerboseSolution

abstract type SolutionType end

"""
    VerboseSolution <: SolutionType

Used to return a [`VerboseSolutionResults`](@ref)

!!! note
    VerboseSolution is designed for debugging purposes
    only, and is not GPU-friendly.
"""
struct VerboseSolution <: SolutionType end

abstract type AbstractSolutionResults{AbstractFloat} end

"""
    VerboseSolutionResults{AT} <: AbstractSolutionResults{AT}

Result returned from `find_zero` when
`VerboseSolution` is passed as the `soltype`.
"""
struct VerboseSolutionResults{AT} <: AbstractSolutionResults{AT}
    "solution ``x^*`` of the root of the equation ``f(x^*) = 0``"
    root::AT
    "indicates convergence"
    converged::Bool
    "error of the root of the equation ``f(x^*) = 0``"
    err::AT
    "number of iterations performed"
    iter_performed::Int
    "solution per iteration"
    root_history::Vector{AT}
    "error of the root of the equation ``f(x^*) = 0`` per iteration"
    err_history::Vector{AT}
end
SolutionResults(soltype::VerboseSolution, args...) =
    VerboseSolutionResults(args...)

"""
    CompactSolution <: SolutionType

Used to return a [`CompactSolutionResults`](@ref)
"""
struct CompactSolution <: SolutionType end

"""
    CompactSolutionResults{AT} <: AbstractSolutionResults{AT}

Result returned from `find_zero` when
`CompactSolution` is passed as the `soltype`.
"""
struct CompactSolutionResults{AT} <: AbstractSolutionResults{AT}
    "solution ``x^*`` of the root of the equation ``f(x^*) = 0``"
    root::AT
    "indicates convergence"
    converged::Bool
end
SolutionResults(soltype::CompactSolution, root, converged, args...) =
    CompactSolutionResults(root, converged)

init_history(::VerboseSolution, ::Type{AT}) where {AT <: AbstractArray} = AT[]
init_history(::CompactSolution, ::Type{AT}) where {AT <: AbstractArray} =
    nothing

function push_history!(
    history::Vector{AT},
    x::AT,
    ::VerboseSolution,
) where {AT <: AbstractArray}
    push!(history, x)
end
function push_history!(
    history::Nothing,
    x::AT,
    ::CompactSolution,
) where {AT <: AbstractArray}
    nothing
end
