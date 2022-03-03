# **Regridding & redistribution**
Regridding (or remapping or interpolation) in the coupler ensures that different component models can communicate fields with each other with a minimal loss of information or field properties. These properties include:
- conservation = no loss of field mass (e.g. mass, energy)
- consistency = no spurious maxima/minima in data (e.g. mass, energy)
- monotonicity = no spurious sign changes (e.g. mass, energy)
- divergence-free property (e.g. wind vector) - usually need to convert to vorticity/divergence for regridding > solve Poisson eqn to get streamfunction / velocity potential > transform back to wind
- high-order derivatives (e.g. wind stress curl, $\nabla \times (c_D U|U|)$ ) - e.g. capitalize on high-order atmos grid, so calculate there
- discontinuous structure (e.g. precipitation)

Depending on which variable is being regridded, different properties, and thus interpolation methods, may be desirable. 

The interpolation methods:
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

# Methods in other climate models
## [ESM4 (FMS)](https://github.com/NOAA-GFDL/ESM4)
- simple (AMIP): all fluxes are regridded as scalar using first order conservative interpolation 
    - [more details](fms_coupler.md)
    - also offers second order conservative and bilinear interpolation
- full (CMIP)

## [ICON](https://code.mpimet.mpg.de/projects/iconpublic/wiki/How%20to%20obtain%20the%20model%20code)
- need licence

## [CESM](https://github.com/ESCOMP/CESM)
- [grids](http://esmci.github.io/cime/versions/master/html/index.html)
- [MCT toolkit](https://www.mcs.anl.gov/research/projects/mct/mct_APIs/index.html)
- (need more investigation)

## mitGCM
- conservative (line integrals < Gauss theorem) 
- bicubic

## OASIS (SPOC)
- [code](https://gitlab.com/cerfacs/oasis3-mct/-/tree/OASIS3-MCT_4.0/examples/spoc/spoc_regridding)
- [docs](documentation)
- [SCRIP doc](https://oasis.cerfacs.fr/wp-content/uploads/sites/114/2021/03/GLOBC_SCRIPusers_1998.pdf)
    - second order conservative (weights from line integrals < Gauss theorem)
    - bilinear
    - bicubic
    - nearest-neighbor (e.g. for categorical data)

## [TempestRemap](https://github.com/ClimateGlobalChange/tempestremap)
- conservative remap (default) without the need for line integrals - mixture of patching and projection to conservative+consistent solutions
- inverse distance waving and bilinear - less visible in code (per comm)

# CliMA strategy

- [current implementation docs](climacore_remap_notes.md)

## 2D Boundary Regrid

### AMIP + LES
- fields to exchange (TBC):
    - atmos > land: heat/moisture fluxes, T
    - land > atmos: roughness lengths, T_sfc, q_sfc,...
    - prescribed ice and ocean > atmos: T_sfc, q_sfc, roughness, albedo
- for all conservative remapping sufficient

#### Optimization
- better treatment of polar vector values (now arbitrarily selected)
- develop divergence-free regridding for vectors
- calculate wind stress curl in (higher order) atmos before sending to ocean (CMIP)
- test best methods for discontinuous output (precip)

## 3D Boundary Regrid
#### Optimization
- staggering (1d) (climacore)
- staggering (3d) (oceananigans)
- develop 3D wind transform 

## Supporting test cases
- implement in sea breeze LES [in_progress]
- implement in convection 2D LES [not_started]
- implement in HS GCM [pending]

# Misc Notes
- other higher than 2nd order regridding may be possible, but for climate modelling unlikely needed, since much larger errors come from elsewhere


# Refs
- [NCAR regridding intro](https://climatedataguide.ucar.edu/climate-data-tools-and-analysis/regridding-overview)
- [ESM notes](https://earthsystemmodeling.org/docs/release/ESMF_6_1_0/esmf_6_1_0_regridding_status.html)
- [Mahadevan et al 2021](https://gmd.copernicus.org/preprints/gmd-2021-323/gmd-2021-323.pdf) - metrics for remapping algorithms

