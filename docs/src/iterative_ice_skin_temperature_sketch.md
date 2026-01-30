# Sketch: Iterative Ice Skin Temperature Using SurfaceFluxes.jl

This document sketches how to add an **iterative ice skin temperature** for sea ice in ClimaCoupler while keeping **SurfaceFluxes.jl** as the sole source of turbulent fluxes. The ice surface temperature is updated using the same **flux balance** logic as in ClimaOcean's `SkinTemperature` (conductive flux through ice vs. turbulent + radiative fluxes at the surface).

---

## 1. Goal

- **SurfaceFluxes.jl**: continues to compute turbulent fluxes (sensible, latent, momentum) from atmospheric state and surface state; no change to its API.
- **Ice skin temperature**: solved iteratively so that at convergence, the surface temperature satisfies the flux balance at the ice–atmosphere interface (turbulent + radiative fluxes balanced by conductive flux through the ice).
- **ClimaOcean**: reuse its flux-balance formula for the ice side (e.g. `SkinTemperature{ClimaSeaIce.ConductiveFlux}`) so we don’t reimplement physics; call into a small, stateless helper if possible.

---

## 2. When to Use Iteration

- **Only for sea ice** when the chosen formulation is “iterative skin”.
- In practice: only when `sim isa ClimaSeaIceSimulation` (or a future sea-ice type that opts in) and a flag/config requests iterative skin (e.g. `use_iterative_ice_skin_temperature = true`).
- All other surfaces (ocean, land, prescribed ice) keep the current single-call path: one `get_field(:surface_temperature)` and one `get_surface_fluxes` call.

---

## 3. High-Level Flow (Single Column)

For each column (or each point on the coupler boundary space) where sea ice is present:

1. **Initial guess** for ice surface temperature:  
   `T_sfc = ice_model.top_surface_temperature` (or previous timestep / previous iterate).
2. **Loop** until convergence (or max iterations):
   - Build surface thermodynamic state from current `T_sfc` (surface density, saturation humidity over ice, etc.).
   - Call **SurfaceFluxes.jl**: `get_surface_fluxes(..., thermo_state_sfc, ...)` → get `F_sh`, `F_lh` (and momentum, etc.). All fluxes positive upward.
   - Build **radiative terms** at the surface (same convention as ClimaOcean):
     - Outgoing LW: `Q_LW_out = ε σ T_sfc^4`
     - Net absorbed SW (e.g. `(1 - α) * SW_d`); net absorbed LW (e.g. `LW_d - Q_LW_out`). So total “flux out of surface” from radiation: `Q_rad_out = Q_LW_out - (1-α)*SW_d - LW_d` or equivalently structure it as in ClimaOcean (Qu, Qd, etc.).
   - **Flux balance**:  
     Given turbulent heat flux `Q_turb = F_sh + F_lh` (positive upward = out of surface) and radiative flux at surface, and ice state (thickness `h`, critical thickness `hc`, conductivity `k`, melting temperature at base `T_base`), compute the **new** `T_sfc` such that conductive flux through the ice balances the total flux at the surface.  
     Use the same formula as ClimaOcean’s `flux_balance_temperature(::SkinTemperature{ClimaSeaIce.ConductiveFlux}, ...)` (see ClimaOcean `interface_states.jl`).
   - **Check convergence**: e.g. `|T_sfc_new - T_sfc| < tol` (and optionally cap iterations).
   - Set `T_sfc = T_sfc_new` and repeat.
3. **Output**: converged `T_sfc` and the corresponding turbulent fluxes from the last SurfaceFluxes call. Use these to update coupler fields and the ice model (same as today).

---

## 4. Alignment with ClimaCoupler and Sea Ice Patterns

The sketch follows the same patterns used elsewhere in ClimaCoupler and its sea ice modules.

### 4.1 Where component-specific flux logic lives

- **FluxCalculator** (in `src/FluxCalculator.jl`) defines the **base** `compute_surface_fluxes!(csf, sim::AbstractSurfaceSimulation, ...)`, which uses `get_field(sim, Val(:surface_temperature))`, calls SurfaceFluxes once, then `update_turbulent_fluxes!(sim, fields)` and adds to `csf`.
- **Component modules extend** `FluxCalculator.compute_surface_fluxes!` for their simulation type when they need different flux logic: `BucketSimulation` in **climaland_bucket.jl**, `ClimaLandSimulation` in **climaland_integrated.jl**. Sea ice currently has no override and uses the base path.
- The iterative ice skin path should be added as a **component extension** in **experiments/ClimaEarth/components/ocean/clima_seaice.jl** (same file that already extends `update_turbulent_fluxes!` and `FieldExchanger.update_sim!`), **not** in `src/FluxCalculator.jl`.

### 4.2 Component provides state via get_field; declares behavior

- Surface state is provided via `Interfacer.get_field(sim, Val(:field_name))`. Ice state for the flux balance (e.g. thickness `h`, critical thickness `hc`, conductivity `k`, melting temperature at base) should be provided the same way: new `get_field` methods in **clima_seaice.jl** for `ClimaSeaIceSimulation`, with remapping from the ice grid done inside those getters.
- Components declare options via `get_field` (e.g. `:roughness_model`). For iterative skin, use e.g. `get_field(sim, Val(:use_iterative_ice_skin))` → `true`/`false`. The component’s `compute_surface_fluxes!` runs the iteration only when `true`; otherwise it can `invoke` the base implementation.

### 4.3 Writing back and tail of compute_surface_fluxes!

- The base builds `fields = (; F_turb_ρτxz, F_turb_ρτyz, F_lh, F_sh, F_turb_moisture)` and calls `update_turbulent_fluxes!(sim, fields)`. The component extends `update_turbulent_fluxes!(sim::ClimaSeaIceSimulation, fields)` in **clima_seaice.jl** to write those fluxes into the ice model. For iterative ice, include converged **T_sfc** in `fields` (e.g. `T_sfc_converged`) and extend the same method to write it into the ice model’s `top_surface_temperature`.
- The base ends with: zero out where `area_fraction ≈ 0`, `update_turbulent_fluxes!(sim, fields)`, then add to `csf` (fluxes, L_MO, ustar, buoyancy_flux). The iterative sea ice override should use the **same tail** so coupler fields stay consistent.

---

## 5. Where It Lives in the Codebase

### 5.1 Sea ice component (clima_seaice.jl)

- **New**: extend `FluxCalculator.compute_surface_fluxes!` for `ClimaSeaIceSimulation` in **experiments/ClimaEarth/components/ocean/clima_seaice.jl** (not in FluxCalculator):

  ```julia
  function compute_surface_fluxes!(
      csf,
      sim::ClimaSeaIceSimulation,
      atmos_sim::AbstractAtmosSimulation,
      thermo_params,
  )
      # ... if use_iterative_ice_skin(sim) ...
      # 1. Get initial T_sfc from sim (e.g. get_field(sim, Val(:surface_temperature)))
      # 2. Get ice state on boundary space: h, hc, k, T_base (new get_field or remap from ice model)
      # 3. Loop:
      #    - Build thermo_state_sfc from T_sfc
      #    - fluxes = get_surface_fluxes(...)  # SurfaceFluxes.jl
      #    - Q_turb = F_sh + F_lh; build Q_rad from csf.SW_d, csf.LW_d, T_sfc, α, ε, σ
      #    - T_sfc_new = ice_flux_balance_temperature(T_sfc, Q_turb, Q_rad, h, hc, k, T_base, ...)
      #    - break if converged
      #    - T_sfc = T_sfc_new
      # 4. Update sim’s stored top_surface_temperature (if desired) and call update_turbulent_fluxes!(sim, fields)
      # 5. Add fluxes to csf (area-weighted) same as now
  end
  ```

- When `get_field(sim, Val(:use_iterative_ice_skin))` is `false` or unimplemented: call the base via `invoke(...)` so current behavior is unchanged. When `true`: run the iteration, then the same tail as the base. The **base** `compute_surface_fluxes!(csf, sim::AbstractSurfaceSimulation, ...)` in **FluxCalculator** stays unchanged.

### 5.2 Ice-Side Flux Balance

Two options:

- **A. Thin wrapper in ClimaOcean**  
  In ClimaOcean’s InterfaceComputations, add a small **stateless** function that takes:
  - Current guess `T_sfc`,
  - Turbulent heat flux (sensible + latent, positive upward),
  - Radiative terms (e.g. Qu, Qd or equivalent),
  - Ice properties: `h`, `hc`, `k`, `T_base` (melting T at base),
  - Thermodynamic parameters (e.g. for saturation over ice),

  and returns the **new** `T_sfc` from one step of the flux balance (or the full fixed-point iteration internally). ClimaCoupler then calls this from `compute_surface_fluxes!` for sea ice and passes in the fluxes from SurfaceFluxes and the radiative terms computed in the coupler.

- **B. Copy formula into ClimaCoupler**  
  Implement in ClimaCoupler a function that mirrors ClimaOcean’s `flux_balance_temperature(::SkinTemperature{ClimaSeaIce.ConductiveFlux}, ...)` (single step: given fluxes and ice state, return new T_sfc). No dependency on ClimaOcean’s internal state types; only the algebraic formula is shared (and documented as such).

Recommendation: **A** if ClimaOcean is willing to expose a small API; **B** if we want to avoid adding a ClimaOcean API and accept one-time formula sync.

---

## 6. Data Flow (Summary)

```
Atmosphere state (from csf)     Ice state (from ClimaSeaIceSimulation)
T_atmos, q_atmos, ρ_atmos,      h, hc, k, T_base (melting),
u_int, v_int, height_int       α, ε, σ
         │                              │
         ▼                              ▼
    ┌─────────────────────────────────────────┐
    │  Iteration (per column / boundary point)│
    │  ┌───────────────────────────────────┐  │
    │  │ T_sfc_guess                       │  │
    │  │   → surface state → SurfaceFluxes │  │  ← SurfaceFluxes.jl
    │  │   → F_sh, F_lh, ...               │  │
    │  │   → Q_turb, Q_rad                 │  │
    │  │   → flux_balance → T_sfc_new      │  │  ← ClimaOcean-style balance
    │  │   → repeat until converged        │  │
    │  └───────────────────────────────────┘  │
    └─────────────────────────────────────────┘
         │
         ▼
    Final T_sfc, F_sh, F_lh, F_turb_ρτxz, F_turb_ρτyz, F_turb_moisture
         │
         ├── update_turbulent_fluxes!(sim, fields)  → ice model
         └── csf.F_* += fluxes * area_fraction      → coupler
```

---

## 7. Interfaces with SurfaceFluxes.jl

- **No change** to SurfaceFluxes.jl’s public API.
- ClimaCoupler still calls the same entry point: e.g. `FluxCalculator.get_surface_fluxes(surface_fluxes_params, u_int, thermo_state_atmos, h_int, u_sfc, thermo_state_sfc, h_sfc, d, config)`.
- The only difference is that for iterative sea ice, this call is made **inside a loop** with an updated `thermo_state_sfc` (and hence updated `T_sfc`) each time; the last call’s fluxes are the ones used.

So: **SurfaceFluxes.jl continues to define the turbulent-flux interface**; the coupler only adds an outer iteration and an ice-side flux balance when the surface is sea ice with iterative skin enabled.

---

## 8. ClimaSeaIceSimulation Extensions

- **Optional**: `get_field(sim, Val(:ice_flux_balance_state))` or similar that returns, on the boundary space (or remapped from the ice grid), the fields needed for the flux balance: `h`, `hc`, conductivity `k`, and base (melting) temperature. If not provided, these can be derived from existing ice model fields and remapping in `compute_surface_fluxes!`.
- **Optional**: a way to write back the converged `T_sfc` into the ice model’s `top_surface_temperature` so the next timestep (or other physics) sees a consistent value; otherwise the iteration is self-contained within the flux computation.

---

## 9. Implementation Order

1. **Flux balance helper**  
   Implement (in ClimaCoupler or ClimaOcean) the single-step flux balance: inputs = (T_sfc, Q_turb, radiative terms, h, hc, k, T_base, …), output = T_sfc_new. Add unit tests against ClimaOcean’s `flux_balance_temperature` for a few (T_sfc, flux) pairs.

2. **Per-point iteration in clima_seaice.jl**  
   Add `FluxCalculator.compute_surface_fluxes!(csf, sim::ClimaSeaIceSimulation, ...)` in **clima_seaice.jl** that: when iterative skin is enabled, gets initial T_sfc and ice state via get_field; runs the loop (SurfaceFluxes → fluxes → radiative terms → flux balance → T_sfc_new); breaks on convergence; then uses the same tail as the base (update_turbulent_fluxes!, add to csf). When iterative skin is disabled, invoke the base implementation.

3. **Remapping / get_field**  
   In **clima_seaice.jl**, add get_field(sim::ClimaSeaIceSimulation, Val(:ice_flux_balance_state)) (or per-field getters) so ice state (h, hc, k, T_base) is available on the coupler boundary space, with remapping from the ice grid done inside the component.

4. **Config / dispatch**  
   Gate the iterative path on `get_field(sim, Val(:use_iterative_ice_skin))` (defined in clima_seaice.jl for ClimaSeaIceSimulation) so that existing runs keep the current single-call behavior unless they opt in.

5. **Integration tests**  
   Run a short coupled case with sea ice and iterative skin, and compare (e.g. T_sfc and surface fluxes) against a single-call run or against ClimaOcean’s own atmosphere–sea-ice coupling if available.

---

## 10. References

- ClimaOcean `InterfaceComputations`:  
  - `compute_interface_state.jl`, `interface_states.jl` (SkinTemperature, flux_balance_temperature),  
  - `atmosphere_sea_ice_fluxes.jl` (how atmosphere–ice fluxes and interface state are computed).
- ClimaCoupler:  
  - `FluxCalculator.compute_surface_fluxes!`,  
  - `ClimaSeaIceSimulation` in `experiments/ClimaEarth/components/ocean/clima_seaice.jl`.
