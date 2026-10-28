#=
# Cell exchange

Each process holds the exchange-grid polygons of its own SE elements and its
own block of the ocean cells, and a polygon's cell may be stored on another
process. A gather from the cells first pulls the values of the cells the
polygons of this process use; a scatter to the cells pushes the partial sums of
these polygons back to the process storing each cell, which adds them in
process order. On one process there is nothing to exchange.
=#

"""
    CellExchange

The ocean cells one process exchanges with the others, numbered by *slot*:
the `n_local` cells it stores, in the storage order of its ocean grid, then the
cells its polygons use from each other process, in process order.

A pull sends the values of the stored slots `send_slots`, grouped by the
receiving process, and receives the other cells into `remote`. A push sends
the partial sums `remote_partial` of the other cells and receives the partial
sums of the stored cells into `returned`, laid out as `send_slots`;
`returned_ptr`/`returned_pos` are the CSR of the positions in `returned` of
each stored cell. `plan` holds the processes and message lengths of both.
"""
struct CellExchange{C, P, VI, VF, G}
    context::C
    plan::P
    n_local::Int
    send_slots::VI
    returned_ptr::VI
    returned_pos::VI
    sent::VF
    remote::VF
    remote_partial::VF
    returned::VF
    pull::G
    push::G
end

function CellExchange(
    context,
    plan,
    n_local,
    send_slots,
    returned_ptr,
    returned_pos,
    sent,
    remote,
    remote_partial,
    returned,
)
    (; send_pids, send_lengths, recv_pids, recv_lengths) = plan
    pull = ClimaComms.graph_context(
        context,
        sent,
        send_lengths,
        send_pids,
        remote,
        recv_lengths,
        recv_pids,
    )
    push = ClimaComms.graph_context(
        context,
        remote_partial,
        recv_lengths,
        recv_pids,
        returned,
        send_lengths,
        send_pids,
    )
    return CellExchange(
        context,
        plan,
        n_local,
        send_slots,
        returned_ptr,
        returned_pos,
        sent,
        remote,
        remote_partial,
        returned,
        pull,
        push,
    )
end

"""
    CellExchange(FT, context, plan, n_local)

A zeroed `CellExchange` on the CPU with `FT` buffers, for the `plan` of
[`plan_cell_exchange`](@ref).
"""
function CellExchange(::Type{FT}, context, plan, n_local) where {FT}
    send_slots = Int32.(plan.send_slots)
    returned_ptr, returned_pos = _build_csr(plan.send_slots, n_local)
    n_send, n_recv = length(send_slots), sum(plan.recv_lengths; init = 0)
    return CellExchange(
        context,
        plan,
        n_local,
        send_slots,
        returned_ptr,
        returned_pos,
        zeros(FT, n_send),
        zeros(FT, n_recv),
        zeros(FT, n_recv),
        zeros(FT, n_send),
    )
end

# The graph contexts are bound to the buffers, so they are rebuilt on the new ones.
Adapt.adapt_structure(to, cx::CellExchange) = CellExchange(
    cx.context,
    cx.plan,
    cx.n_local,
    Adapt.adapt(to, cx.send_slots),
    Adapt.adapt(to, cx.returned_ptr),
    Adapt.adapt(to, cx.returned_pos),
    Adapt.adapt(to, cx.sent),
    Adapt.adapt(to, cx.remote),
    Adapt.adapt(to, cx.remote_partial),
    Adapt.adapt(to, cx.returned),
)

"""
    plan_cell_exchange(context, grid_oc, n_cells, used_cells) -> (slot, owned, plan)

The slot of each of the `n_cells` global ocean cells (numbered `i` fastest) that
this process stores or whose number is in `used_cells`, the global numbers of the
cells it stores, in slot order, and the processes and message lengths of the
exchange. Collective over `context`.
"""
function plan_cell_exchange(context, grid_oc, n_cells, used_cells)
    pid, nprocs = ClimaComms.mypid(context), ClimaComms.nprocs(context)
    owned, owner = CR.local_indices(grid_oc, n_cells, pid)
    needed = [Int[] for _ in 1:nprocs]
    for c in sort!(unique(used_cells))
        push!(needed[owner(c)], c)
    end
    incoming = get_ConservativeRegriddingCCExt().incoming_requests(context, needed)
    plan = CR.exchange_plan(needed, incoming, owned, pid)

    slot = Dict(c => s for (s, c) in enumerate(owned))
    s = length(owned)
    for p in plan.recv_pids, c in needed[p]
        slot[c] = s += 1
    end
    return slot, owned, plan
end

@kernel function _copy_slots_kernel!(dst, src, slots)
    i = @index(Global)
    @inbounds dst[i] = src[slots[i]]
end

"""
    pull_cells!(cx::CellExchange, cell_values)

Receive into `cx.remote` the values of the cells this process uses from the
others, sending them the values in `cell_values` of the cells they use.
"""
function pull_cells!(cx::CellExchange, cell_values)
    backend = KernelAbstractions.get_backend(cx.sent)
    isempty(cx.sent) || launch_kernel!(
        _copy_slots_kernel!,
        backend,
        length(cx.sent),
        cx.sent,
        cell_values,
        cx.send_slots,
    )
    ClimaComms.start(cx.pull)
    ClimaComms.finish(cx.pull)
    return cx.remote
end

"""
    push_partials!(cx::CellExchange)

Send `cx.remote_partial` to the processes storing those cells and receive the
partial sums of the stored cells into `cx.returned`.
"""
function push_partials!(cx::CellExchange)
    ClimaComms.start(cx.push)
    ClimaComms.finish(cx.push)
    return cx.returned
end
