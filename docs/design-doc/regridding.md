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
- Example: 1D FV > SE
    - 1st order monotone FV map: $R_{ij} = J^{ov}_{ij} / J^t_i$, where $\Omega_{ij}^{ov} = \Omega_i^t \cap \Omega_i^t$ is a local overlapping region which belongs to a set of $\Omega^{ov}$ so that $\Omega_{j,i}^{s,t} = \sum_{i=1}^{f^{s,t}} \Omega_{ij}^{ov}$
    1. define continuous function of target and source dofs ($\tilde{J}^s_j(x)$ and $\tilde{J}^t_i(x)$)
        - for FV $\tilde{J}$ = 1 within face and 0 elsewhere
        - for SE - continuous function associated with each dof are computed using a tensor product of polynomials
    2. Enforce that $\sum \tilde{J} = J$
    3. define $R_{ij} = \frac{1}{J^t_i}\int \tilde{J}^s_j(x) \tilde{J}^t_i(x) dx$ (numerical integration can be non-conservative. In that case, we can project coefficients into the space of consrevative maps to recover conservation)
    - see [Ullrich & Taylor 15 Appendix](https://journals.ametsoc.org/view/journals/mwre/143/6/mwr-d-14-00343.1.xml ) for a more concrete example

# ClimaCore
- for cartesian regridding
- R = LinearRemap 
    - `linear_remap_op = overlap_weights / local_weights(target)`
        - overlap_weights (lengths or areas), depending on source and target:
            - 1D: SE{1} > SE{1}
            - 2D: SE{1} > SE{1}

            - 1D: SE{N} > SE{1}
            - 1D: SE{1} > SE{N}
            - where SE{1} ~ FV
        - local_weights
            - wj = space.local_geometry.WJ
- remap!
    - `mul!(vec(parent(target_field)), R.map, vec(parent(source_field)))`

## TempestRemap
- TempestRemap uses a quadrature-based approach to produce a “first guess” operator that is then projected onto the space of conservative and consistent solutions using a novel least-squares formulation. The resulting method **avoids the need for line integrals and can be used to guarantee conservation and consistency** (and, if desired, monotonicity) of the linear map.

# Refs 
- [Ullrich & Taylor 15](https://journals.ametsoc.org/view/journals/mwre/143/6/mwr-d-14-00343.1.xml )
- [Ullrich et al. 16](https://journals.ametsoc.org/view/journals/mwre/144/4/mwr-d-15-0301.1.xml)

Q: in CG cartesian do not need to do the projection to conserative and consistent space because integration is exact, right? What about spherical?
Q: SE{1} in CG caused instability, correct?

# Plan
- Cartesian [Ben,Valeria,Lenka-Jan/Feb] 
    - finish off the current test cases
    - 2D SE cases
- Spherical [Lenka,Valeria,Ben-Jan/Feb]
    - read TempestRemap documentation + familiarise with code
    - get familiar with mesh differences
    - adapter design + implementation
    - prototype julia version [long-term]
