#### NonlinearSolvers.jl
include("tolerance_types.jl")
include("solution_types.jl")

using ForwardDiff
export solve!

# Main entry point: Dispatch to specific method
"""
    solve!(
        method::AbstractNonlinearSolverMethod{FT},
        soltype::SolutionType = CompactSolution(),
        tol::Union{Nothing, AbstractTolerance} = nothing,
        maxiters::Union{Nothing, Int} = 10_000,
    )

Solve the non-linear system given
 - `method` the numerical method
 - `soltype` the solution type (`CompactSolution` or `VerboseSolution`)
 - `tol` the stopping tolerance
 - `maxiters` the maximum number of iterations to perform
"""
function solve!(
    method::AbstractNonlinearSolverMethod{FT},
    soltype::SolutionType = CompactSolution(),
    tol::Union{Nothing, AbstractTolerance} = nothing,
    maxiters::Union{Nothing, Int} = 1,
) where {FT <: FTypes, F <: Function}
    if tol === nothing
        tol = ResidualTolerance{eltype(FT)}(sqrt(eps(FT)))
    end
    return solve!(method, method_args(method)..., soltype, tol, maxiters)
end

##### NewtonsMethodAD

export NewtonsMethodAD

"""
    NewtonsMethodAD(f!::F!, x_init::A) where {F!, A <: AbstractArray}

A non-linear system of equations type.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct NewtonsMethodAD{FT, F!, A, JA} <: AbstractNonlinearSolverMethod{FT}
    "Function to find the root of"
    f!::F!
    "Initial guess"
    x_init::A
    "Storage"
    x1::A
    "Storage"
    J::JA
    "Storage"
    J⁻¹::JA
    "Storage"
    F::A
    function NewtonsMethodAD(f!::F!, x_init::A) where {F!, A }
        x1 = similar(x_init)
        J = similar(x_init, (length(x_init), length(x_init)))
        J⁻¹ = similar(J)
        F = similar(x_init)
        JA = typeof(J)
        FT = eltype(x_init)
        return new{FT, F!, A, JA}(f!, x_init, x1, J, J⁻¹, F)
    end
end

method_args(m::NewtonsMethodAD) = (m.x_init, m.x1, m.f!, m.F, m.J, m.J⁻¹)
function solve!(
    ::NewtonsMethodAD,
    x0::AT,
    x1::AT,
    f!::F!,
    F::FA,
    J::JA,
    J⁻¹::J⁻¹A,
    soltype::SolutionType,
    tol::AbstractTolerance{FT},
    maxiters::Int,
) where {FA, J⁻¹A, JA, F! <: Function, AT, FT}

    x_history = init_history(soltype, AT)
    F_history = init_history(soltype, AT)
    if soltype isa VerboseSolution
        f!(F, x0)
        ForwardDiff.jacobian!(J, f!, F, x0)
        push_history!(x_history, x0, soltype)
        push_history!(F_history, F, soltype)
    end
    for i in 1:maxiters
        @show i
        @show x0
        f!(F, x0)
        @show "pre jacobian"
        @show J
        @show f!
        @show F
        @show x0
        ForwardDiff.jacobian!(J, f!, F, x0)
        x1 .= x0 .- J \ F
        @show "post jacobian"
        push_history!(x_history, x1, soltype)
        push_history!(F_history, F, soltype)
        if tol(x0, x1, F)
            return SolutionResults(
                soltype,
                x1,
                true,
                F,
                i,
                x_history,
                F_history,
            )
        end
        x0 = x1
    end
    return SolutionResults(
        soltype,
        x0,
        false,
        F,
        maxiters,
        x_history,
        F_history,
    )
end


function newton_simple(f,Df,x0,epsilon,max_iter)
    """Approximate solution of f(x)=0 by Newton's method.

    Parameters
    ----------
    f : function
        Function for which we are searching for a solution f(x)=0.
    Df : function
        Derivative of f(x).
    x0 : number
        Initial guess for a solution f(x)=0.
    epsilon : number
        Stopping criteria is abs(f(x)) < epsilon.
    max_iter : integer
        Maximum number of iterations of Newton's method.

    Returns
    -------
    xn : number
        Implement Newton's method: compute the linear approximation
        of f(x) at xn and find x intercept by the formula
            x = xn - f(xn)/Df(xn)
        Continue until abs(f(xn)) < epsilon and return xn.
        If Df(xn) == 0, return None. If the number of iterations
        exceeds max_iter, then return None.

    Examples
    --------
    >>> f = lambda x: x**2 - x - 1
    >>> Df = lambda x: 2*x - 1
    >>> newton(f,Df,1,1e-8,10)
    Found solution after 5 iterations.
    1.618033988749989
    """
    xn = x0
    for n in 1:max_iter
        fxn = f(xn)
        if abs(fxn) < epsilon
            println("Found solution after",n,"iterations.")
            return xn
        end
        Dfxn = Df(xn)
        if Dfxn == 0
            print("Zero derivative. No solution found.")
            return nothing
        end
        xn = xn - fxn/Dfxn
        println("Exceeded maximum iterations. No solution found.")
        return nothing
    end
end



# f!(F, x0)
# ForwardDiff.jacobian!(J, f!, F, x0)
# x1 .= x0 .- J \ F