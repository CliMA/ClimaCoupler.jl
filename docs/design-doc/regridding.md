# **Regridding & redistribution**

- we use linear remapping, so that:
 $\phi^{target} = R (\phi^{source}) $

# Properties
### conservative
$$
    \sum_{j=1}^{f^s}\psi_j^sJ_j^s = \sum_{j=1}^{f^t}\psi_i^tJ_i^t 
$$
### consistent
$$
\mathbf{1}^t=\mathbf{R}\mathbf{1}^s
$$ 
### monotonicity-preserving
$$
\forall (i,j) \in [1,...,f^t] \times [1,...,f^s] R_{ij} \geq 0
$$
Note:
- if R is conservative and consistent, the source and target meshes must have the same area. 
- monotone linear maps can be at most 2nd order if face size of target ~ face size of source. Otherwise 1st order. 

# General steps
1. obtain source and target meshes: $\Omega_j^s$ and $\Omega_i^t$
2. generate the overlap mesh: $\Omega_{ij}^{ov} = \Omega_j^s \cap \Omega_i^t$
3. define continuous function of target and source dofs ($\tilde{J}^s_j(x)$ and $\tilde{J}^t_i(x)$)
    - FV: $\tilde{J}$ = 1 within face and 0 elsewhere
    - SE: continuous function associated with each dof are computed using a tensor product of polynomials
4. define $R_{ij} = \frac{1}{J^t_i}\int \tilde{J}^s_j(x) \tilde{J}^t_i(x) dx$ 
    - if numerical integration non-conservative, project coefficients into the space of conservative maps to recover conservation)
    
# ClimaCore usage
- in the intermediate term, we develop our own regridding for cartesian domains, and for regridding on a sphere we use the published TempestRemap (which does not fully support Cartesian regridding):
    - Cartesian: [ClimaCore regridding](https://github.com/CliMA/ClimaCore.jl/blob/main/src/Operators/remapping.jl) 
        - remap operator generation
            - `local_weights = space.local_geometry.WJ`
            - `linear_remap_op = overlap_weights / local_weights(target)`
        -  remap! operator application
            - `mul!(vec(parent(target_field)), R.map, vec(parent(source_field)))`
    - Spherical: [ClimaCoreTempestRemap](https://github.com/CliMA/ClimaCore.jl/tree/ln/exodus/lib/ClimaCoreTempestRemap)
        - this helps write and read files formatted in TempestRemap format
        - we use the TempestRemap_jll helper to call TempestRemap directly from Julia to generate the mapping weights from source and target meshes, and the map application can be done within the ClimaCoupler
        - Links:
            - [TempestRemap](https://github.com/ClimateGlobalChange/tempestremap) repo
            - [TempestRemap_jll](https://github.com/JuliaPackaging/Yggdrasil/tree/master/T/TempestRemap) - package that automatically configures TempestRemap and allows calls directly from Julia
            - [our notes](tempestremap_notes.md) on Tempest summarizing relevant info

# Plan
- Cartesian [Ben,Valeria,Lenka-Jan/Feb] 
    - finish off the current test cases
    - 2D SE cases
- Spherical [Lenka,Valeria,Ben-Jan/Feb]
    - read TempestRemap documentation + familiarise with code [2d]
    - get familiar with mesh differences [3d]
    - adapter design + implementation [5d]
    - prototype julia version [long-term]
