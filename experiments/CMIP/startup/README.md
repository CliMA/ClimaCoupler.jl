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
  inference only, which on this workload is the smaller half (inference 18.8 min vs LLVM
  codegen 34.5 min in a 61-minute setup).
* `cmip_startup.jl` -- full harness, samples the setup and step phases.
* `cmip_startup_nosample.jl` -- timing only. Prefer this unless a breakdown is needed: on a
  run this size `Profile.print` adds ~1 h of post-processing, and sampling correlated with
  hangs when two heavy jobs ran concurrently.
* `analyze_flat.py` -- renders the flat profile dumps into an attribution table. Use the
  main-thread root frame as the denominator; the profiler samples all threads, so total
  sample count exceeds wall time.

Load phases are deliberately not sampled: the profiler interacts badly with package loading.
