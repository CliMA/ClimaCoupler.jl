# WeatherBench 10d experiment log

Local run log for 10-day ProgEDMF 0M WeatherBench forecasts on Derecho.
Newest first. Package `NEWS.md` is the public changelog; this is just the
lab notebook.

Shared unless a row says otherwise:

- Coupler: `copies2/ClimaCoupler.jl`, branch `weatherbench2_pr_derecho2`
- Atmos path: `/glade/u/home/cchristo/clima/copies2/ClimaAtmos.jl`
- IC: `…/initial_conditions_0p25deg_dev/no_hz_stretch`, `era5_ic_full_pressure: true`
- Grid: `h_elem: 30`, `z_elem: 130`, `z_max: 70000`, `dt`/`dt_cpl`: `40secs`
- Network: `col_w3h_ctx_julia.nc`
- Start date in the scored plots: `20200101`

## 2026-10-04 — ML stronger (gain 0.75 + tq)

YAML: `weatherbench_progedmf_0m_10d_mlcorr.yml`  
Output: `…/weatherbench_10d_mlcorr_strong/20200101/`  
Atmos: `ml_internal_corr` @ `0c27a48ce`  
Job: `7709593` (`cc_wx_prog_20200101`)  
PBS/logs: `generated/batch_10d_mlcorr_strong/`

Paired control (same Atmos, no ML): job `7708949` → `…/weatherbench_10d_control_main/20200101/`

| Key | Value |
|---|---|
| `ml_correction_variables` | `tq` (was `t`) |
| `ml_correction_gain` | `0.75` (was `0.5`; 1.0 is the offline-matched next step) |
| `ml_correction_t_cap` | `1.5` K/h (was `0.5`) |
| taper | 50 → 20 hPa (`p_full: 5000`, `p_zero: 2000`; was 100 → 50) |
| `ml_correction_q_cap` | default `2e-4` (leave unless this run looks clipped) |
| `ml_correction_dt` | `1hours` |

Gain and humidity applied together (not split).

## 2026-10-04 — next ML (in progress)

YAML: `weatherbench_progedmf_0m_10d_mlcorr.yml`  
Output: `…/weatherbench_10d_mlcorr_strong/`  
Atmos: `ml_internal_corr` @ `0c27a48ce`

Already set in the YAML (this commit of the notebook):

- `ml_correction_t_cap: 1.5` (was 0.5 K/h)
- taper `p_full: 5000`, `p_zero: 2000` (50 → 20 hPa; was 100 → 50 hPa)

Still to edit before submit (user): **gain** and **`ml_correction_variables`**
(`t` vs `tq`). Leave a new output-dir suffix if you change those.

## 2026-10-03 — isolation controls (single-run A/B vs the “good” machine)

Not ML. Same 10d control physics, with Atmos main pieces reverted one at a time.
Scored against the first ML/control pair; late-lead differences were noise.

| Label | Atmos branch | What changed vs `ml_internal_corr` | Output |
|---|---|---|---|
| control (latest) | `ml_internal_corr` | none | `…/weatherbench_10d_control/` |
| control cf041 | `test/cf4869_cm041` | revert CF #4869 + CM 0.43 → 0.41 / Params 1.1.17; `dt` 50s | `…/weatherbench_10d_control_cf041/` |
| control noTKE | `test/cf4869_cm041_nolTKE` | cf041 + revert `l_TKE` | (as run) |

Conclusion: keep **latest** Atmos (`ml_internal_corr`), not the reverts.

## 2026-10-03 — first online ML (baseline)

YAML: first `weatherbench_progedmf_0m_10d_mlcorr.yml`  
Output: `…/weatherbench_10d_mlcorr/20200101/`  
Atmos: `ml_internal_corr` (then ~`c72de942a` / later `0c27a48ce` after rebase)

| Key | Value |
|---|---|
| `ml_correction_variables` | `t` |
| `ml_correction_gain` | `0.5` |
| `ml_correction_t_cap` | `0.5` K/h (Atmos default) |
| taper | 100 → 50 hPa (Atmos default) |
| `ml_correction_dt` | `1hours` |

Result (one start date): stable. T850 and q850 better than control from day 1;
Z500 mainly after day ~7. q850 improved even with T-only (thermo coupling).
Weaker than offline, as expected at gain 0.5 and no q tendency.

## 2026-10-03 — first 10d control on this stack

YAML: `weatherbench_progedmf_0m_10d_control.yml`  
Output: `…/weatherbench_10d_control/20200101/`  
Atmos: `ml_internal_corr`, `ml_correction` unset.

Paired with the first ML run above.
