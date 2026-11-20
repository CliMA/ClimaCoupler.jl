# Startup-time measurement harness

Phase-resolved accounting of launch -> first coupled step, splitting wall time into Julia
inference, host codegen, GPU kernel compilation, GC and I/O.

```bash
CLIMACOMMS_DEVICE=CUDA COUPLER_PERF_OUT=<dir> \
  julia --project=experiments/CMIP experiments/CMIP/startup/cmip_startup.jl \
  --config_file <config>
```

* `phase_probe.jl` -- records, per phase, wall time, `Base.cumulative_compile_time_ns()`
  (exact Julia compile time) and a sampled category breakdown. Sampling with C frames is the
  only mechanism that sees *both* host codegen and GPU kernel compilation; SnoopCompile sees
  inference only, which on this workload is the smaller half.
* `cmip_startup.jl` -- full harness, samples the setup and step phases.
* `cmip_startup_nosample.jl` -- timing only. Prefer this unless a breakdown is needed:
  sampling a run this size adds substantial post-processing time.
* `analyze_flat.py` -- renders the flat profile dumps into an attribution table. Use the
  main-thread root frame as the denominator; the profiler samples all threads, so total
  sample count exceeds wall time.

Load phases are deliberately not sampled: the profiler interacts badly with package loading.

Measurements from this harness are collected in `experiments/warmup/precompile_info.md`.
