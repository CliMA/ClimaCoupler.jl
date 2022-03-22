# **Regridding & redistribution**
Regridding (or remapping or interpolation) in the coupler ensures that different component models can communicate information (fluxes, states, parameters) with each other with a minimal loss of accuracy at the highest possible speed. 

Depending on the regridded variable, several properties may need to be satisfied:
- conservation = no loss of field mass (e.g. mass, energy)
- consistency = no spurious maxima/minima in data (e.g. mass, energy)
- monotonicity = no spurious sign changes (e.g. mass, energy)
- divergence-free property (e.g. wind vector) - usually need to convert to vorticity/divergence for regridding > solve Poisson eqn to get streamfunction / velocity potential > transform back to wind
- high-order derivatives (e.g. wind stress curl, $\nabla \times (c_D U|U|)$ ) - e.g. capitalize on high-order atmos grid, so calculate there
- discontinuous structure (e.g. precipitation)

Different interpolation methods can be applied to achieve the above properties, for example:
- bilinear and bicubic
    - easy to implement
    - not always conservative
- higher-order patch recovery ([patch](https://www.cesm.ucar.edu/events/workshops/ws.2007/presentations/Neckels_OCN.pdf))
    - higher order accuracy
    - more difficult to implement especially if need conservation
    - a divergence-free preserving extension not fully developed yet
- global line integrals from Gauss theorem 
    - conservative
    - divergence-free preserving
    - BUT global potential function is required 

# Methods in Other Climate Models
##  NOAA-GFDL / ESM4 (FMS)
- [code](https://github.com/NOAA-GFDL/ESM4)
- simple (AMIP): all fluxes are regridded as scalar using first order conservative interpolation 
    - [our notes](esm4.md)
    - also offers second order conservative and bilinear interpolation
- full (CMIP)

## DWD / ICON
- [code](https://code.mpimet.mpg.de/projects/iconpublic/wiki/How%20to%20obtain%20the%20model%20code)
- need licence

## NCAR / CESM
- [code](https://github.com/ESCOMP/CESM)
- [grids](http://esmci.github.io/cime/versions/master/html/index.html)
- [MCT toolkit](https://www.mcs.anl.gov/research/projects/mct/mct_APIs/index.html) -> replaced by [CMEPS](https://github.com/ESCOMP/CMEPS) general mediator
- [our notes](cesm.md)

## mitGCM
- [code](https://github.com/MITgcm/MITgcm)
- conservative (line integrals < Gauss theorem), bicubic

# General Regridding Packages


## OASIS (SPOC)
- [code](https://gitlab.com/cerfacs/oasis3-mct/-/tree/OASIS3-MCT_4.0/examples/spoc/spoc_regridding), [SCRIP doc](https://oasis.cerfacs.fr/wp-content/uploads/sites/114/2021/03/GLOBC_SCRIPusers_1998.pdf)
- interpolation methods
    - second order conservative (weights from line integrals < Gauss theorem)
    - bilinear
    - bicubic
    - nearest-neighbor (e.g. for categorical data)

## TempestRemap
- [code](https://github.com/ClimateGlobalChange/tempestremap)
- [our notes](tempestremap.md)
- interpolation types
    - conservative remap (default) without the need for line integrals - mixture of patching and projection to conservative+consistent solutions
    - inverse distance waving and bilinear - less visible in code (per comm)

# CliMA strategy
- for the initial AMIP milestone, it is sufficient to apply a first order conservative linear remap (e.g. as in TempestRemap or FMS). We thus use TempestRemap for generation of spherical map weights. Since TempestRemap does not directly support Cartesian domains, we develop an additional Cartesian functionality in ClimaCore mimicking methods used in TempestRemap.

## Implementation
- Cartesian
    - [cartesian regridding design](cartesian.md)
    - [ClimaCore regridding implementation ](https://github.com/CliMA/ClimaCore.jl/blob/main/src/Operators/remapping.jl) 
        - remap operator generation
            - `local_weights = space.local_geometry.WJ`
            - `linear_remap_op = overlap_weights / local_weights(target)`
        -  remap! operator application
            - `mul!(vec(parent(target_field)), R.map, vec(parent(source_field)))`
- Spherical
    - [spherical regridding design](spherical.md)
    - [ClimaCoreTempestRemap implementation](https://github.com/CliMA/ClimaCore.jl/tree/main/lib/ClimaCoreTempestRemap) (provides wrappers for [TempestRemap](https://github.com/ClimateGlobalChange/tempestremap) )

 

# CliMA Test Case Plan

## 2D Boundary Regrid

### AMIP + LES
- conservative remapping sufficient for exchage of fields /states:
    - atmos > land: heat/moisture fluxes, wind, momentum fluxes (?)
    - land > atmos: roughness lengths, T_sfc, q_sfc, roughness
    - prescribed ice and ocean > atmos: T_sfc, q_sfc, roughness, albedo
- grid box partition done automatically in Cartesian
- topography

### Dynamic Aquaplanet
- conservative remapping sufficient for exchage of fields /states:
    - atmos > oceananigans 
        - need to map weights to oceananigans grid (need ordering of the nodes)

### Optimization
- more efficient storage of the sparse matrix [in progress]
- grid box splitting
- better treatment of polar vector values (now arbitrarily selected)
- develop divergence-free regridding for vectors
- calculate wind stress curl in (higher order) atmos before sending to ocean (CMIP)
- test best methods for discontinuous output (precip)
- more unit tests (e.g. CC>TR meshes w precision errors)

## 3D Boundary Regrid
### Optimization
- staggering (1d) (climacore)
- staggering (3d) (oceananigans)
- develop 3D wind transform 

# Supporting test cases
- implement in sea breeze LES [in_progress]
    - cartesian regridding (FE <-> FV)
    - ClimaCore + ClimaAtmos + ClimaSim interface (+ Land)
    - topography
- implement in convection 2D LES [not_started]
    - direct comparison with literature
- implement in HS/BCwave GCM w diffusion [in_progress]
    - spherical regridding
    - ClimaCore + ClimaAtmos + ClimaSim interface

# Misc Notes
- other higher than 2nd order regridding may be possible, but for climate modelling unlikely needed, since much larger errors come from elsewhere


# Refs
- [NCAR regridding intro](https://climatedataguide.ucar.edu/climate-data-tools-and-analysis/regridding-overview)
- [ESM notes](https://earthsystemmodeling.org/docs/release/ESMF_6_1_0/esmf_6_1_0_regridding_status.html)
- [Mahadevan et al 2021](https://gmd.copernicus.org/preprints/gmd-2021-323/gmd-2021-323.pdf) - metrics for remapping algorithms

