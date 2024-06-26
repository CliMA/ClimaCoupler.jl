import ClimaCoupler: Checkpointer, Interfacer

"""
    Base.show(io::IO, dict::Dict)

This prints the keys and values of a Dict in sorted order.
"""
function Base.show(io::IO, dict::Dict)
    for k in sort!(collect(keys(dict)))
        println(io, " $k => $(dict[k])")
    end
end

# user callbacks
"""
    checkpoint_sims(cs::CoupledSimulation, _)

This is a callback function that checkpoints all simulations defined in the current coupled simulation.
"""
function checkpoint_sims(cs::Interfacer.CoupledSimulation, _)
    for sim in cs.model_sims
        if Checkpointer.get_model_prog_state(sim) !== nothing
            t = Dates.datetime2epochms(cs.dates.date[1])
            t0 = Dates.datetime2epochms(cs.dates.date0[1])
            Checkpointer.checkpoint_model_state(sim, cs.comms_ctx, Int((t - t0) / 1e3), output_dir = cs.dirs.artifacts)
        end
    end
end
