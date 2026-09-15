#=
# Flux snapshots

Host-side copies of the coupler's fluxes on every grid they pass through, written
to JLD2 for `Plotting.plot_flux_snapshot`. Components opt in by extending
`capture_flux_snapshot` and `flux_snapshot_geometry`.
=#

"""
    capture_flux_snapshot(sim, csf)

Copy to the host every flux `sim` exchanges, on each grid it passes through, as
a nested `NamedTuple` of arrays grouped by component (`atmos`, `ocean`, ...).
`csf` holds the coupler fields. Returns `nothing` if `sim` has nothing to capture
(the default).
"""
capture_flux_snapshot(sim, csf) = nothing

"""
    flux_snapshot_geometry(sim, csf)

Plain-array geometry of the grids the arrays from
[`capture_flux_snapshot`](@ref) live on, or `nothing` (the default).
"""
flux_snapshot_geometry(sim, csf) = nothing

import Dates
import JLD2
import ClimaCoupler: Interfacer

"""
    flux_snapshot_dir(coupler_output_dir)

Directory flux snapshots are written to.
"""
flux_snapshot_dir(coupler_output_dir) = joinpath(coupler_output_dir, "fluxes")

# Write a nested NamedTuple of arrays as JLD2 paths: `(; a = (; b = x))` → "a/b".
function _write_groups!(file, prefix, value::NamedTuple)
    for (name, child) in pairs(value)
        _write_groups!(file, isempty(prefix) ? string(name) : "$prefix/$name", child)
    end
    return nothing
end
_write_groups!(file, path, value) = (file[path] = value; nothing)

"""
    write_flux_snapshot(model_sims, csf, dir, date; suffix = "")
    write_flux_snapshot(cs::Interfacer.CoupledSimulation; suffix = "")

Capture every component's flux snapshot and write them to one file,
`<dir>/fluxes_<yyyy-mm-ddTHHMMSS><suffix>.jld2`, grouped by component. The first
call in `dir` also writes `exchange_grid_geometry.jld2`. Returns the snapshot
path, or `nothing` if no component had anything to capture.
"""
function write_flux_snapshot(model_sims, csf, dir, date; suffix = "")
    groups = (;)
    for sim in values(model_sims)
        captured = capture_flux_snapshot(sim, csf)
        isnothing(captured) || (groups = merge(groups, captured))
    end
    isempty(groups) && return nothing

    mkpath(dir)
    geometry_path = joinpath(dir, "exchange_grid_geometry.jld2")
    if !isfile(geometry_path)
        for sim in values(model_sims)
            geometry = flux_snapshot_geometry(sim, csf)
            isnothing(geometry) && continue
            JLD2.jldopen(file -> _write_groups!(file, "", geometry), geometry_path, "w")
            break
        end
    end

    stamp = Dates.format(date, "yyyy-mm-ddTHHMMSS")
    path = joinpath(dir, "fluxes_$(stamp)$(suffix).jld2")
    JLD2.jldopen(path, "w") do file
        file["date"] = string(date)
        _write_groups!(file, "", groups)
    end
    return path
end

write_flux_snapshot(cs::Interfacer.CoupledSimulation; suffix = "") = write_flux_snapshot(
    cs.model_sims,
    cs.fields,
    flux_snapshot_dir(cs.dir_paths.coupler_output_dir),
    Interfacer.current_date(cs, cs.t[]);
    suffix,
)

const _COMBINED_TURBULENT_FLUXES =
    (:F_sh, :F_lh, :F_turb_moisture, :F_turb_ρτxz, :F_turb_ρτyz)

"""
    turbulent_fluxes_have_nan(csf)

Whether any combined turbulent flux the atmosphere receives contains a NaN. It
reduces on the device, so it is cheap enough to run every coupling step. A NaN in
any surface's fluxes reaches these fields within the step.
"""
turbulent_fluxes_have_nan(csf) =
    any(name -> any(isnan, parent(getproperty(csf, name))), _COMBINED_TURBULENT_FLUXES)

"""
    NaNFluxSnapshot()

Callback that writes one flux snapshot, suffixed `_nan`, on the first coupling
step where [`turbulent_fluxes_have_nan`](@ref) is true, and does nothing after.
"""
struct NaNFluxSnapshot
    fired::Base.RefValue{Bool}
end
NaNFluxSnapshot() = NaNFluxSnapshot(Ref(false))

function maybe_write_nan_snapshot!(callback::NaNFluxSnapshot, model_sims, csf, dir, date)
    callback.fired[] && return nothing
    turbulent_fluxes_have_nan(csf) || return nothing
    callback.fired[] = true
    path = write_flux_snapshot(model_sims, csf, dir, date; suffix = "_nan")
    if isnothing(path)
        @warn "NaN in the coupler turbulent fluxes at $date"
    else
        @warn "NaN in the coupler turbulent fluxes at $date; wrote a flux snapshot" path
    end
    return path
end

(callback::NaNFluxSnapshot)(cs) = maybe_write_nan_snapshot!(
    callback,
    cs.model_sims,
    cs.fields,
    flux_snapshot_dir(cs.dir_paths.coupler_output_dir),
    Interfacer.current_date(cs, cs.t[]),
)
