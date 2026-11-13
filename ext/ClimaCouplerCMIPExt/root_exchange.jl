#=
# Exchange on the root rank

With `exchange_on_root`, the ocean, the sea ice, the exchange grid and the
polygon fluxes live on the root rank, on a single-process replica of the
boundary space's mesh, element order and quadrature. The coupler fields the
exchange reads are gathered into the replica and the fields it produces are
scattered back. ClimaCore gives each rank a contiguous, rank-ordered block of
the element order, so the element blocks of a field gathered in rank order are
in the replica's order.
=#

import MPI

# The coupler fields the exchange reads and one scratch field
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
    :SW_d_diffuse,
    :LW_d,
    :P_liq,
    :P_snow,
    :scalar_temp,
)

"""
    gather_to_root!(remapping, csf, name::Symbol)
    gather_to_root!(remapping, csf, names::Tuple)

Gather the coupler field `csf.name` into `remapping.exchange_fields.name` on the
root rank.
"""
function gather_to_root!(remapping, csf, name::Symbol)
    staging = remapping.local_scratch
    staging .= getproperty(csf, name)
    context = ClimaComms.context(axes(staging))
    recvbuf = if ClimaComms.iamroot(context)
        exchange_field = getproperty(remapping.exchange_fields, name)
        MPI.VBuffer(parent(exchange_field), remapping.rank_counts)
    else
        nothing
    end
    KernelAbstractions.synchronize(KernelAbstractions.get_backend(parent(staging)))
    MPI.Gatherv!(parent(staging), recvbuf, context.mpicomm; root = 0)
    return nothing
end
gather_to_root!(remapping, csf, names::Tuple) =
    foreach(name -> gather_to_root!(remapping, csf, name), names)

"""
    scatter_from_root!(field, remapping, name::Symbol)

Scatter `remapping.exchange_fields.name` from the root rank into the
boundary-space `field`.
"""
function scatter_from_root!(field, remapping, name::Symbol)
    staging = remapping.local_scratch
    context = ClimaComms.context(axes(staging))
    sendbuf = if ClimaComms.iamroot(context)
        exchange_field = getproperty(remapping.exchange_fields, name)
        MPI.VBuffer(parent(exchange_field), remapping.rank_counts)
    else
        nothing
    end
    KernelAbstractions.synchronize(KernelAbstractions.get_backend(parent(staging)))
    MPI.Scatterv!(sendbuf, parent(staging), context.mpicomm; root = 0)
    field .= staging
    return nothing
end
