# **Baroclinic Wave + Slab**

# Atmosphere
The momentum equations are in the advective form, and tracers in the consevative form, namely:

- Density:
$$ \frac{\partial \rho}{\partial t} + \nabla \cdot ({\rho \vec{u}})= 0 $$

- Momentum (flux form):
$$ \frac{\partial \vec{u_h}}{\partial t} + \vec{u} \cdot \nabla \vec{u_h} = - \frac{1}{\rho}\nabla_h p
+ \frac{\partial}{\partial z} K_v \frac{\partial}{\partial z} \vec{u_h}
$$
$$ \frac{\partial w}{\partial t} + \vec{u} \cdot \nabla w= 
- \frac{1}{\rho}\frac{\partial p}{\partial z}
- \nabla_z \Phi 
+ \frac{\partial}{\partial z} K_v \frac{\partial}{\partial z} w 
$$

- Total energy:
$$ \frac{\partial \rho e_{tot}}{\partial t} + \nabla \cdot (\rho h_{tot} \vec{u}) = \frac{\partial}{\partial z} K_v \frac{\partial}{\partial z} h_{tot}
$$

where the total specific enthalpy and total specific  energy are
$$ h_{tot} =  e_{tot} + \frac{p}{\rho}  \,\,\,\,\,\,\,\, \,\,\,\,\,\,\,\, e_{tot} = c_v T + \Phi + \frac{1}{2}\vec{u}^2 
$$
(note that $h_{tot} \neq h = c_vT + p/\rho = c_p T$, the specific enthalpy in the thermodynamic sense), $\Phi = gz$ is the geopotential,
$u_h$ is the horizontal velocity vector, $w$ the vertical velocity, $\rho$ the density, $p$ pressure, $K_v$ the vertical diffusivity (assumed constant here). 

## Boundary conditions (BCs)
- We implement BCs similarly to other climate models. 
    - First-order fluxes (i.e., advective fluxes) are always set to zero, corresponding to the *free-slip* and *impenetrable* BC, where:
    $$
    w = 0 \,\,\,\,\,\,\, \partial_t w = 0 \,\,\,\,\,\,\, \nabla \times\vec{u_h}=0 \,\,\,\,\,\,\, \nabla \cdot \vec{\rho u_h}=0  \,\,\,\,\,\,\, \nabla \cdot \rho h_{tot} \vec{u_h}=0
    $$
    - Second-order fluxes (i.e., diffusive fluxes)
        - `No Flux`: By default we have *impenetrable* or *insulating* BCs (no second-order fluxes) at all boundaries. 
        - `Bulk Formula`: Applied to tracers (e.g., temperature and moisture), this imposes a boundary fluxes (e.g., sensible and latent heat) calculated using the bulk aerodynamic formulae  using prescribed surface values of ($T_{sfc}$ and $q_{sfc}^{sat}$).  At the surface, the bulk sensible heat flux formula for total enthalpy essentially replaces the above: 
            $$ (K_v \rho \partial_z h_{tot})_{sfc}$$
            For **total energy**, we have two choices:
            - 1. enthalpy flux: 
                $$ (K_v \rho \partial_z h_{tot})_{sfc} \rightarrow 
                \hat{n} \cdot  \rho C_H ||u||^{1} (h^1- h_{sfc})   
                = F_S
                $$
            - 2. sensible (and latent) heat flux. The sensible heat flux is: 
                $$ (K_v \rho \partial_z h_{tot})_{sfc} \rightarrow 
                \hat{n} \cdot C_H c_{pd} ρ^{1} ||u||^{1} (T^{1} - T_{sfc}) 
                +  \hat{n} \cdot C_H ρ^{1} ||u||^{1} (\Phi^{1} - \Phi_{sfc})      
                = F_S
                $$
                where $^{1}$ corresponds to the lowest model level, $C_H$ is the dimensionless thermal transfer coefficient, $c_{pd}$ is the specific heat capacity for dry air $||u||$ the wind speed. This is the *bulk turbulent sensible heat flux* parameterization, and $F_S$ is positive when atmosphere receives energy from the surface. 
                The contribution of the kinetic energy is usually O(1e4) smaller and is neglected, but it can be added to F_S as:
                $$
                F_{S_{tot}} = F_S + \hat{n} \cdot C_D ρ^{1} ||u||^{1} (\vec{u_h}^{1})^2 
                $$        
         - `Drag Law`: essentially the bulk formula for momentum       
            $$ \frac{\partial}{\partial z} K_v \frac{\partial}{\partial z} \vec{u_h} \rightarrow 
            \hat{n} \cdot C_D ρ^{1} ||u||^{1} \vec{u_h}^{1}          
            = F_M
            $$
        - `Coupled Bulk Formula`: same as `Bulk Formula`, but surface quantities (e.g. $T_{sfc}$) are passed from the state of the neighboring model.

    - The diffusive fluxes are applied via the `vertical_diffusion` ClimaAtmos model sub-component. To apply boundary fluxes without diffusion in the atmospheric interior, the viscosity coefficient needs to be set to zero: $ν = FT(0)$.
    - We use `SurfaceFluxes.jl` to deal with the Monin Obukhov calculations

## Current setup
- total energy and momentum
    - at z=0:
        - ρe: `Coupled Bulk Formula` latent heat flux + sensible heat flux (will be combined into enthalpy flux formulation in SurfaceFluxes.jl)
        - uh: `Drag Law` 
    - the the top:
        - $F_S = F_M = 0$
    - values
        - simple setup: $C_D = C_H = 0.001$, interior diffusivity is set to $\nu = 5$ m^2/s
        - ClimaAtmos setup ∀ p > p_pbl :    $C_D = C_E = 0.0044  exp(-\frac{(p_{pbl} - p)}{p_{strato}}^2)$
            - where $p_{pbl} = 8e4$ Pa, $p_{strato} = 1e4$ Pa
 
- All other boundary fluxes are set to 0.

## Initial conditions
- we initialize with a perturbation in a balanced background state, as in:
https://climate.ucdavis.edu/pubs/UMJS2013QJRMS.pdf

# Heat Slab
The slab solves for temperature in a single layer, whose tendency is the accumulated fluxes divided by the coupling timestep plus a parameterisation of the internal processes, $G$.
$$
\rho c h_s  \, \partial_t T_{sfc} =  - F_{integ}  / \Delta t_{coupler}  
$$

# Distribution
## Julia multithreading
- we can run using multiple threads if the command line argument `enable_threading` is true (e.g., run with `julia --project --threads 8`). 

## MPI via ClimaComms.jl 
- for AMIP we want the surface columns to be on the same processors as the atmos columns. Since all all surface domains will be on masked spheres, all of them can inherit the same distributed horizontal space from the atmos model.
- for this we need to run the `sbatch_job.sh` script, which sets up the `CLIMACORE_DISTRIBUTED` environment variable and job specifications, and runs the coupler_driver with `mpiexec` 

## Regridding
- not needed for AMIP. ClimaCoreTempestRemap can be easily re-introduced for a single processor, but will require more work for MPI runs.  

# Tests
## Conservation
- this uses the `sum` of `ClimaCore/Fields/mapreduce.jl`, which produces a sum weighted by the area Jacobian. 
    - one can easily check this by summing a field of ones using the domain's space:
        ```
        field_of_ones = ones(center_space)
        sum(field_of_ones) ≈ (4*pi*domain_radius^2) * domain_height 
        ```

## Performance
- using `@elapsed` to measure the walltime of the coupling loop
- 1. strong scaling
    - increasing the number of precessing elements (MPI processes or threads)

- 2. weak scaling
    - increasing the number of precessing elements (MPI processes or threads) with job size (vertical resolution)

- 3. comparison to stand-alone atmos 
    - using the original ClimaAtmos driver (using `solve!`)

## Physical correctness
- run the default for 20 days

# Prescribed SST and Sea Ice
- We simply prescribe SSTs from a file as `T_sfc`. As for sea ice, we will follow GFDL's [AMIP setup](https://pcmdi.llnl.gov/mips/amip/home/Documentation/20gfdl.html#RTFToC31) and use prescribed sea ice concentrations and a constant ice thickness, $h_{i} = 2m$ ice thickness, while solving for $T_{sfc}$:
$$
\frac{dT_{sfc}}{dt} = - frac{h_i(F_{atm} - F_{conductive}) / k_i}
$$ 
where
$$
F_{conductive} = \frac{k_i (T_{base} - {T_sfc})}{h_{i}}
$$
with the thermal conductivity of ice, $k_i = 2$ W m$^{-2}$ K$^{-1}$, and $T_{base} = 273.16$ K. For now we use an Euler timestepper (and use $T_{sfc}$ of the previous timestep), though this may be solved implicitly in the future. 

## Data source
- https://gdex.ucar.edu/dataset/158_asphilli.html
    - MODEL.SST.HAD187001-198110.OI198111-202203.nc
    - MODEL.ICE.HAD187001-198110.OI198111-202203.nc
- N.B.: the [pcmdi link](https://pcmdi.llnl.gov/mips/amip/details/amipbc_dwnld.php), used in most AMIP papers, is broken


# NB:
- first coupled iteration does not call rhs!
- slab `T_sfc` gets huge numbers when using `SSPRK33`. ok with `Euler`
- do not init global fields with mpi context

# TODO
- ClimaAtmos: sub in newest CA interface
- interface: 
    - add coupler specific abstractions
- fluxes: re-enable different ways to calculate / accumulate fluxes (at overy coupler timestep; at every atmos timestep via callback; via specification of an additional variable) 
- formalize conservation/physocal/performance tests: add error threshold and exception, interval, show option, and make a general interface for it
- SurfaceFluxes: combine LHF and SHF into enthalpy flux formulation to avoid division by zero
- Temporally varying SSTs/sea ice: https://pcmdi.llnl.gov/mips/amip/details/


# References
- [Kang et al 2021](https://arxiv.org/abs/2101.09263)
- [kth.se blog for strong and weak scaling](https://www.kth.se/blogs/pdc/2018/11/scalability-strong-and-weak-scaling/)




