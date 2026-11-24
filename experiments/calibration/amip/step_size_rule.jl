import EnsembleKalmanProcesses as EKP
import LinearAlgebra
import Statistics

"""
    EffectiveEntriesController(effective_entries)

`EKP.DataMisfitController` with the number of loss entries replaced by `effective_entries`.
Use it with observations whose noise covariance has been scaled by
`entries / effective_entries` (`make_observation_vector` in `generate_observations.jl`).

DataMisfitController sets the step to `Δt = max(M / (2 mean(Φ)), sqrt(M / (2 var(Φ))))`, with
`Φ` half of each member's misfit and `M` the number of entries, so the update weighs the data
as `M` independent entries. Entries closer together than the decorrelation length repeat each 
other's error, and the update is then too confident by about the number of entries in one 
correlated patch. Scaling the noise covariance cannot fix this, because `Φ` scales with it 
and the step cancels it. This controller uses `effective_entries` in place of `M`, so the step
no longer cancels the tempering.

Changing only the step is not enough. Each UKI iteration inflates the ensemble covariance by
`1 + Δt` before the update, so the posterior it converges to is set by the information `Γ`
claims, whatever the step. On an untempered native-grid loss the step-only rule gives a
first update of the right width and then collapses within about five iterations.

It never terminates, like `DataMisfitController(terminate_at = 1e6)`.
"""
struct EffectiveEntriesController{FT <: AbstractFloat} <: EKP.LearningRateScheduler
    effective_entries::FT
end

function EKP.calculate_timestep!(
    ekp::EKP.EnsembleKalmanProcess,
    g::AbstractMatrix,
    Δt_new::Union{Nothing, AbstractFloat},
    scheduler::EffectiveEntriesController,
)
    isnothing(Δt_new) ||
        @info "EffectiveEntriesController chooses its own step, ignoring Δt_new = $Δt_new"
    diff = g .- reshape(EKP.get_obs(ekp), :, 1)
    X = EKP.lmul_obs_noise_cov_inv(ekp, diff)
    Φ = [LinearAlgebra.dot(diff[:, j], X[:, j]) / 2 for j in 1:size(g, 2)]
    m = scheduler.effective_entries
    Δt = max(m / (2 * Statistics.mean(Φ)), sqrt(m / (2 * Statistics.var(Φ))))
    push!(EKP.get_Δt(ekp), Δt)
    return nothing
end

"""
    effective_entries(observation, observables)

Total number of independent values `observation` counts as across `observables`: each
observable contributes `o.effective_entries` when set, or `1.0` for an untempered 1-entry
scalar (`isnothing(o.effective_entries)` with `length(s) == 1`). Errors if a multi-entry
observable (`length(s) > 1`) leaves `isnothing(o.effective_entries)`, which would mix an
untempered multi-cell map into a tempered loss.
"""
function effective_entries(observation::EKP.Observation, observables::AbstractVector)
    samples = EKP.get_samples(observation)
    length(observables) == length(samples) || error(
        "$(length(observables)) observables for $(length(samples)) covariance blocks",
    )
    total = 0.0
    for (s, o) in zip(samples, observables)
        if isnothing(o.effective_entries)
            length(s) == 1 || error(
                "Multi-entry observable $(o.name) ($(length(s)) entries) must set effective_entries when used with EffectiveEntriesController",
            )
            total += 1.0
        else
            total += Float64(o.effective_entries)
        end
    end
    return total
end
