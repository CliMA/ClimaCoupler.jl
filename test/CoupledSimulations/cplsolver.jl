using Test
using ClimaCoupler
using ClimateMachine.GenericCallbacks, ClimateMachine.ODESolvers
using Dates, Unitful, LinearAlgebra

#= Coupled Oscillators
https://scholar.harvard.edu/files/schwartz/files/lecture3-coupled-oscillators.pdf

Two coupled oscillators start from rest but with at least
one perturbed from equilibrium (hence, motion). In our simulation set up,
each individual mass is only aware of its own position, not that of the
other. Through the coupler, the force each mass exerts on the other is passed
to the other mass. A coupler is certainly not necessary to solve the problem,
we could set up a simple matrix system and also have an analytical solution,
however this is a good demonstration of the spectrum of coupling set ups:
from tight monolithic coupling (matrix system) to the black box (our coupler).

| k   κ   k |
|~~~◯~~~◯~~~|
|   m   m   |
|   1   2   |

F11 = -kx1 - κx1
F21 = κx2
F11 + F21 = m x1''

F22 = -kx2 - κx2
F12 = κx1 
F22 + F12 = m x2''

dx/dt = v; dv/dt = a = F/m
=#
@testset "Coupled Timestepping" begin
    FT = Float64

    function exact_sol(t, k, κ, m, x1_0, x2_0)
        ωs, ωf = sqrt(k / m), sqrt((k + 2κ) / m)
        As, Af = x1_0 + x2_0, x1_0 - x2_0
        x1 = 0.5 * (As * cos.(ωs .* t) + Af * cos.(ωf .* t))
        x2 = 0.5 * (As * cos.(ωs .* t) - Af * cos.(ωf .* t))
        return x1, x2
    end

    mutable struct RecordState
        u::Vector{Any}
        counter::Int
    end
    RecordState() = RecordState(Any[], 0)

    rs1, rs2 = RecordState(), RecordState()

    mutable struct SpringRHS{T}
        v::T
        F_in::T
    end

    function (f!::SpringRHS)(du, u, p, t, args...; increment = false)
        f!.v .= u[1]
        du[1] = u[2] .+ increment * du[1] # u[1] is position
        du[2] = (-(p.k + p.κ) * f!.v[1] + f!.F_in[1]) / p.m .+ increment * du[2] # u[2] is velocity
    end

    p = (k = FT(4), κ = FT(2), m = FT(1))
    cdt = FT(0.01)
    endtime = 20
    t = 0.0:(cdt * 10):endtime
    nsteps1, nsteps2 = 2, 3

    timestepper(f, state, nsteps) = LSRK54CarpenterKennedy(f, state, dt = cdt / nsteps, t0 = 0.0)

    function main(::Type{FT}) where {FT}
        x1, x2 = [FT(1), FT(0)], [FT(0), FT(0)] # initial position & velocity
        f1! = SpringRHS([FT(0)], [FT(0)])
        f2! = SpringRHS([FT(0)], [FT(0)])

        # cb1, cb2: record total iterations for each model
        cb1 = EveryXSimulationSteps(1) do
            rs1.counter += 1
        end
        cb2 = EveryXSimulationSteps(1) do
            rs2.counter += 1
        end

        mass1 = CplModel(nothing, f1!, nothing, x1, timestepper(f1!, x1, nsteps1), nsteps1, cb1)
        mass2 = CplModel(nothing, f2!, nothing, x2, timestepper(f2!, x2, nsteps2), nsteps2, cb2)

        coupler = CplState()
        coupler_register!(coupler, :F12, [p.κ * x1[1]], nothing, DateTime(0), u"N")
        coupler_register!(coupler, :F21, [p.κ * x2[1]], nothing, DateTime(0), u"N")

        comp1 = (pre_step = prestep1, component_model = mass1, post_step = poststep1)
        comp2 = (pre_step = prestep2, component_model = mass2, post_step = poststep2)
        components = (mass1 = comp1, mass2 = comp2)
        cpl_solver = CplSolver(component_list = components, coupler = coupler, coupling_dt = cdt, t0 = 0.0)

        # set up coupler callbacks
        # cb3, cb4: record state info with time.
        push!(rs1.u, deepcopy(mass1.state[1]))
        cb3 = EveryXSimulationSteps(10) do
            push!(rs1.u, deepcopy(mass1.state[1]))
        end
        push!(rs2.u, deepcopy(mass2.state[1]))
        cb4 = EveryXSimulationSteps(10) do
            push!(rs2.u, deepcopy(mass2.state[1]))
        end
        cbvector = (cb3, cb4)

        return cpl_solver, cbvector
    end

    # presteps: recv force from other mass
    function prestep1(csolver)
        mass1 = csolver.component_list.mass1.component_model
        mass1.discretization.F_in .= coupler_get(csolver.coupler, :F21, nothing, DateTime(0), u"N")
    end

    function prestep2(csolver)
        mass2 = csolver.component_list.mass2.component_model
        mass2.discretization.F_in .= coupler_get(csolver.coupler, :F12, nothing, DateTime(0), u"N")
    end

    # poststeps: send force upon other mass
    function poststep1(csolver)
        mass1 = csolver.component_list.mass1.component_model
        coupler_put!(csolver.coupler, :F12, p.κ * mass1.state[1], nothing, DateTime(0), u"N")
    end

    function poststep2(csolver)
        mass2 = csolver.component_list.mass2.component_model
        coupler_put!(csolver.coupler, :F21, p.κ * mass2.state[1], nothing, DateTime(0), u"N")
    end

    function run(cpl_solver, params, numberofsteps, cbvector)
        solve!(nothing, cpl_solver, params; numberofsteps = numberofsteps, callbacks = cbvector)
    end

    simulation, cbvector = main(FT)
    nsteps = Int(endtime / cdt)
    run(simulation, p, nsteps, cbvector)

    # subiteration execution
    @test rs2.counter == (nsteps2 / nsteps1) * rs1.counter
    # we solve the problem well enough (toggle plots on to check)
    @test norm(exact_sol(t, p..., 1.0, 0.0)[1] .- rs1.u) <= 1e-1
    @test norm(exact_sol(t, p..., 1.0, 0.0)[2] .- rs2.u) <= 1e-1
end

# Quality Assurance Plotting
# using Plots
# p1 = plot(t, exact_sol(t, p..., 1.0, 0.0)[1], label="exact",
#         xlabel="t", ylabel="x1", title="mass 1")
# plot!(t, rs1.u, label="sim", ls=:dash)
# p2 = plot(t, exact_sol(t, p..., 1.0, 0.0)[2], label="exact",
#         xlabel="t", ylabel="x2", title="mass 2")
# plot!(t, rs2.u, label="sim", ls=:dash)
# plot(p1, p2, layout = (2, 1), legend = true, lw=3)
