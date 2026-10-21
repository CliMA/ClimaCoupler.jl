#=
# Exchange space

The exchange grid, the polygon fluxes and the ocean and sea-ice models live on
one process. On a single process the exchange space is the boundary space. On
several, the root rank builds a single-process replica of the boundary space's
mesh, element order and quadrature; the coupler fields the exchange reads are
gathered into it and the fields it produces are scattered back. ClimaCore gives
each rank a contiguous, rank-ordered block of the element order, so the
element blocks of a field gathered in rank order are in the replica's order.
=#

import MPI

"""
    exchange_space(boundary_space)

The space the exchange grid is built on: `boundary_space` on one process; on
several, a single-process replica of it on the root rank and `nothing` elsewhere.
"""
function exchange_space(boundary_space)
    context = ClimaComms.context(boundary_space)
    context isa ClimaComms.SingletonCommsContext && return boundary_space
    ClimaComms.iamroot(context) || return nothing
    grid = CC.Spaces.grid(boundary_space)
    topology = grid.topology
    replica_topology = CC.Topologies.Topology2D(
        ClimaComms.SingletonCommsContext(ClimaComms.device(context)),
        topology.mesh,
        topology.elemorder,
        nothing,
        topology.orderindex,
    )
    replica_grid = CC.Grids.SpectralElementGrid2D(
        replica_topology,
        grid.quadrature_style;
        enable_bubble = grid.enable_bubble,
        discretization = grid.discretization,
    )
    return CC.Spaces.SpectralElementSpace2D(replica_grid)
end

# The coupler fields the exchange reads from the atmosphere, the fluxes it
# returns, and one scratch field.
exchange_field_names() = (
    :T_atmos,
    :q_tot_atmos,
    :q_liq_atmos,
    :q_ice_atmos,
    :ρ_atmos,
    :u_int,
    :v_int,
    :height_int,
    :height_sfc,
    :SW_d,
    :LW_d,
    :P_liq,
    :P_snow,
    :F_turb_ρτxz,
    :F_turb_ρτyz,
    :F_sh,
    :F_lh,
    :F_turb_moisture,
    :scalar_temp,
)

"""
    gather_to_root!(remapping, csf, name::Symbol)
    gather_to_root!(remapping, csf, names::Tuple)

Copy the coupler field `csf.name` onto the exchange space, into
`remapping.exchange_fields.name`.
"""
function gather_to_root!(remapping, csf, name::Symbol)
    field = getproperty(csf, name)
    context = ClimaComms.context(axes(field))
    if context isa ClimaComms.SingletonCommsContext
        getproperty(remapping.exchange_fields, name) .= field
        return nothing
    end
    staging = remapping.local_scratch
    staging .= field
    recvbuf = if ClimaComms.iamroot(context)
        exchange_field = getproperty(remapping.exchange_fields, name)
        MPI.VBuffer(parent(exchange_field), remapping.rank_counts)
    else
        nothing
    end
    ClimaComms.device_synchronize(ClimaComms.device(context))
    MPI.Gatherv!(parent(staging), recvbuf, context.mpicomm; root = 0)
    return nothing
end
gather_to_root!(remapping, csf, names::Tuple) =
    foreach(name -> gather_to_root!(remapping, csf, name), names)

"""
    scatter_from_root!(field, remapping, name::Symbol)

Copy `remapping.exchange_fields.name` from the exchange space into the
boundary-space `field`.
"""
function scatter_from_root!(field, remapping, name::Symbol)
    context = ClimaComms.context(axes(field))
    if context isa ClimaComms.SingletonCommsContext
        field .= getproperty(remapping.exchange_fields, name)
        return nothing
    end
    staging = remapping.local_scratch
    sendbuf = if ClimaComms.iamroot(context)
        exchange_field = getproperty(remapping.exchange_fields, name)
        MPI.VBuffer(parent(exchange_field), remapping.rank_counts)
    else
        nothing
    end
    ClimaComms.device_synchronize(ClimaComms.device(context))
    MPI.Scatterv!(sendbuf, parent(staging), context.mpicomm; root = 0)
    field .= staging
    return nothing
end
