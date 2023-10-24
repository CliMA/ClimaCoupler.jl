# general callbacks pilot

abstract type CouplerCallback end

@kwdef struct HourlyCallback{FT} <: CouplerCallback
    dt::FT = FT(1) # hours
    func::Function = nothing
    ref_date::Array = [Dates.DateTime(0)]
end

function checkpoint_func(cs, ::CouplerCallback)
    for sim in cs.model_sims
        if get_model_state_vector(sim) !== nothing
            t = Dates.datetime2epochms(cs.dates.date[1])
            t0 = Dates.datetime2epochms(cs.dates.date0[1])
            checkpoint_model_state(sim, cs.comms_ctx, Int( (t - t0)/ 1e3), output_dir = cs.dirs.artifacts)
        end
    end
end

function trigger_callback!(cs::CoupledSimulation, cb::CouplerCallback)
    current_date = cs.dates.date[1]
    if current_date >= cb.ref_date[1]
        ClimaComms.iamroot(cs.comms_ctx) ? @show(cs.dates.date[1]) : nothing
        cb.func(cs, cb)
        cb.ref_date[1] = cs.dates.date[1] + Dates.Hour(cb.dt)
    end
end

