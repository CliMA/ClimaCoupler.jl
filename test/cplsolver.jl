using Test
using CouplerMachine
using ClimateMachine.GenericCallbacks, ClimateMachine.ODESolvers
using Dates, Unitful

#= Coupled Oscillators
https://scholar.harvard.edu/files/schwartz/files/lecture3-coupled-oscillators.pdf

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
=#
FT = Float64

function exact_sol(t, k, κ, m, x1_0, x2_0)
    ωs, ωf = sqrt(k/m), sqrt((k + 2κ)/m)
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
# function (f!::SpringRHS)(du, u, p, t, α=1.0, β=0.0)
#     du[1] = α * 1.0 + β * du[1]
#     du[2] = α * 1.0 + β * du[2]
# end

# function (f!::SpringRHS)(du, u, p, t; increment=false)
#     du[1] = 1.0 + increment * du[1]
#     du[2] = 1.0 + increment * du[2]
# end


p = (
    k = FT(4),
    κ = FT(2),
    m = FT(1),
)
cdt = FT(.02)
nsteps1, nsteps2 = 1, 1

timestepper(f, state, nsteps) = LSRK54CarpenterKennedy(f, state, dt = cdt / nsteps, t0 = 0.0)

function main(::Type{FT}) where {FT}
    x1, x2 = [FT(1), FT(0)], [FT(0), FT(0)] # initial position & velocity
    f1! = SpringRHS([FT(0)], [FT(0)])
    f2! = SpringRHS([FT(0)], [FT(0)])

    mass1 = CplModel(nothing, f1!, nothing, x1, timestepper(f1!, x1, nsteps1), nsteps1)
    mass2 = CplModel(nothing, f2!, nothing, x2, timestepper(f2!, x2, nsteps2), nsteps2)

    coupler = CplState()
    coupler_register!(coupler, :F12, [p.κ * x1[1]], nothing, DateTime(0), u"N")
    coupler_register!(coupler, :F21, [p.κ * x2[1]], nothing, DateTime(0), u"N")

    comp1 = (pre_step = prestep1, component_model = mass1, post_step = poststep1)
    comp2 = (pre_step = prestep2, component_model = mass2, post_step = poststep2)
    components = (mass1 = comp1, mass2 = comp2)
    cpl_solver = CplSolver(
        component_list = components,
        coupler = coupler,
        coupling_dt = cdt,
        t0 = 0.0,
    )

    push!(rs1.u, deepcopy(mass1.state[1]))
    cb1 = EveryXSimulationSteps(5) do
        push!(rs1.u, deepcopy(mass1.state[1]))
    end
    push!(rs2.u, deepcopy(mass2.state[1]))
    cb2 = EveryXSimulationSteps(5) do
        push!(rs2.u, deepcopy(mass2.state[1]))
    end
    cb3 = EveryXSimulationSteps(1) do
        rs1.counter += 1
    end
    cb4 = EveryXSimulationSteps(1) do
        rs2.counter += 1
    end
    cbvector = (cb1, cb2, cb3, cb4)

    return cpl_solver, cbvector
end

function run(cpl_solver, params, numberofsteps, cbvector)
    solve!(
        nothing,
        cpl_solver,
        params;
        numberofsteps = numberofsteps,
        callbacks = cbvector,
    )
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

simulation, cbvector = main(FT)
endtime = 50
nsteps = Int(endtime / cdt)
run(simulation, p, nsteps, cbvector)

t = 0.0:cdt*5:endtime
# using Plots
# plot(t, exact_sol(t, p..., 1.0, 0.0)[1], label="exact")
# plot!(t, rs1.u, label="sim")