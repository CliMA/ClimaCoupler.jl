## rlut_pigroups (2026-08-19/20, 6 iterations)

Config `experiments/calibration/amip/config/rlut_pigroups.jl`,
run dir `amip_calibration_rlut_pigroups`.

# Parameters

1. entr_param_vec_E2 — N(0, 0.5) on [−5, 5]
2. entr_param_vec_E3 — N(0, 0.5) on [−5, 5]
3. entr_param_vec_E6 — N(0.3, 0.1) on [0, 1]

# Observational data

1. Upwelling longwave radiation at TOA (rlut, CERES `toa_lw_all_mon`)

- October 2010, the same month every iteration (no minibatching over years);
  7-day spinup from the 2010-09-24 initial condition, subseasonal mode
- Covariance estimated from October 2000, ..., 2025 (26 samples, rank 25)

# Preprocessing

1. Regrid to the model's output grid — 144 × 72 (2.5°).
2. No latitude window
3. Coarsen data to 5 degrees by computing cos(latitude)-weighted average
4. Land and ocean (no mask); n = 2592

# Covariance matrix

1. SVDplusD covariance matrix
2. Model error scale 0.05 for diagonal matrix (floor ≈ 11.2 W/m²)
3. Regularization

# Outcome

Loss didn't change. Suspected that the model error scale is too high.

---

## rlut_pigroups_ocean (2026-08-20, 4 iterations)

Config `experiments/calibration/amip/config/rlut_pigroups_ocean.jl`,
run dir `amip_calibration_rlut_pigroups_ocean`.

# Parameters

Identical to rlut_pigroups (E2, E3 ~ N(0, 0.5); E6 ~ N(0.3, 0.1)).

# Observational data

1. Upwelling longwave radiation at TOA (rlut, CERES `toa_lw_all_mon`)

- October 2010, the same month every iteration; 7-day spinup from the
  2010-09-24 initial condition, subseasonal mode
- Covariance estimated from October 2000, ..., 2025

# Preprocessing

1. Regrid to the model's output grid — 144 × 72 (2.5°).
2. No latitude window
3. Coarsen data to 5 degrees by computing cos(latitude)-weighted average
4. **Ocean only**: `ClimaAnalysis.apply_landmask` (threshold 0.5) as the last
   step before the sample collection; n = 1696 of 2592

# Covariance matrix

1. SVDplusD covariance matrix
2. **Model error scale 0.02** for diagonal matrix (floor ≈ 4.5 W/m², ≈ 1×
   the interannual spread — the well-specified point)
3. Regularization

# Outcome

Global bias −0.99 increase to −0.07 W/m² with RMSE 10.24 decrease to 10.07.
NaNs in iteration 1.

---

## lwcre_pigroups_ocean (2026-08-21, 3 iterations)

Config `experiments/calibration/amip/config/lwcre_pigroups_ocean.jl`,
run dir `amip_calibration_lwcre_pigroups_ocean`.

# Parameters

Identical to rlut_pigroups (E2, E3 ~ N(0, 0.5); E6 ~ N(0.3, 0.1)).

# Observational data

1. Longwave cloud radiative effect (lwcre = rlutcs − rlut, CERES)

- October 2010, the same month every iteration; 7-day spinup from the
  2010-09-24 initial condition, subseasonal mode
- Covariance estimated from October 2000, ..., 2025

# Preprocessing

1. Regrid to the model's output grid — 144 × 72 (2.5°).
2. No latitude window
3. Coarsen data to 5 degrees by computing cos(latitude)-weighted average
4. Ocean only (`apply_landmask`, threshold 0.5); n = 1696

# Covariance matrix

1. SVDplusD covariance matrix
2. Model error scale 0.02
3. Regularization

# Outcome

Model error scale is too low, so the ensemble collapse in the first update.
The loss increases.

---

## rlut_son_ocean (2026-08-21/23, 4 iterations)

Config `experiments/calibration/amip/config/rlut_son_ocean.jl`,
run dir `amip_calibration_rlut_son_ocean`. First AMIP-mode calibration
(evolving observed SST/sea ice) and first seasonal-mean target.
Run in two launches at clean iteration boundaries (1–2, then 3–4) because
122-day members exceed the 12 h worker walltime.

# Parameters

1. entr_param_vec_E2 — **N(0, 0.3)** on [−5, 5]
2. entr_param_vec_E3 — **N(0, 0.3)** on [−5, 5]
3. entr_param_vec_E6 — N(0.3, 0.1) on [0, 1]

Change std dev from 0.3 to 0.5

# Observational data

1. Upwelling longwave radiation at TOA (rlut, CERES `toa_lw_all_mon`)

- **SON 2010 seasonal mean** (time average of the September, October and
  November monthly means), the same season every iteration; 1-month spinup
  from a 2010-08-01 start, so the simulation covers Aug–Nov and the spinup
  month is discarded by the season window
- Covariance estimated from SON 2000, ..., 2025 (26 seasonal means, rank 25)

# Preprocessing

1. Regrid to the model's output grid — 144 × 72 (2.5°).
2. No latitude window
3. Coarsen data to 5 degrees by computing cos(latitude)-weighted average
4. Ocean only (`apply_landmask`, threshold 0.5); n = 1696
5. Seasonal mean

# Covariance matrix

1. SVDplusD covariance matrix
2. Model error scale 0.02
3. Regularization

# Outcome

One ensemble member ends with NaNs in the first iteration.
RMSE decrease from 10.3 to 9.87 and bias decrease from 2.06 to 1.61.
