# Coupler stepping

There are two modes in which we can run the coupler:
- sequential - models are called in series (e.g. ECMWF IFS)
![](figures/coupling_leapfrog.png)
- concurrent - models are called in parallel (e.g. OASIS-MCT)
    - parallelisation tools
        - multi-node: via MPI (e.g., each model on a different node)
        - single-node: via julia multithreading and shared arrays
    - layout - workload optimisation optimisation
    - lagging - to avoid race conditions
![](figures/coupling_explicit.png)
