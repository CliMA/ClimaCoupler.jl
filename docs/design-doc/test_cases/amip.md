# **Atmosphere-only Model Intercomparison Project (AMIP)** 
- simulations using only the atmospheric model component of comprehensive climate models, used to compare atmospheric dynamics in isolation (note that traditionally, land would be part of the atmospheric model, so its processes are included)
- prescribed processes from monthly observations (e.g. 1979-2014 in CMIP6):
    - SSTs and sea ice concentrations 
    - irradiance, orbital forcings
    - CO2 (1 global value), ozone (and other GHG) monthly climatology (spatially varying
    - aerosols (optional)
- topography + land-sea mask (impl details)

# AMIP Cases
- historical
- amip4xCO2 - rapid cloud response (part of CFMIP; Bony et al 11)
- amip4K - dynamical response of atmos on cloud feedback (part of CFMIP; Bony et al 11)
- related individual experiments
    - more targetted prescription of different land components ([Ackerley et al 18](https://gmd.copernicus.org/articles/11/3865/2018/))

# Atmos
- microphysics & clouds
- convection & BL turbulence
- radiation
- aerosols (optional)

# Land 
- implement both bucket and full soil-vegetation-atmosphere transfer scheme (SVAT) - see Bonan, [AMIP II report](https://pcmdi.github.io/mips/amip/DIAGSUBS/sp12.html#Henderson-Sellers%20et%20al.%201995)

# Our Test Case Hierarchy
- SCM
    - boundary fluxes (moisture)
    - soil: (hydro +) energy 
        - bucket (water store)
        - full soil + hydro 
- sea breeze
    - spatial coupling, topography (no canopy)
- gcm
    - slab
    - full land + hydro




# Refs
- [PCMDI requirements](https://pcmdi.github.io/mips/amip/requirements.html) (incl. output diagnostics) 
- [NACCAP](https://www.narccap.ucar.edu/about/index.html) - assessment of regional regional climate models driven by AOGCMs over North America
- [AMIP BCs generator](https://github.com/PCMDI/amipbcs)
- [PCMDI GH](https://pcmdi.github.io/mips/amip/home/overview.html)
- [PCMDI pointer to SST and SI dsta](https://pcmdi.github.io/mips/amip/SST.html)
- [PCMDI recommendded O3 climatology](https://pcmdi.github.io/mips/amip/AMIP2EXPDSN/OZONE/OZONE1/ozone1.html)

# SCRAP
        - Prescribe SSTs and sea ice concentrations using observations (e.g. 1979-2014 in CMIP6), irradiance, orbital, CO2 (1 value), ozone (and other GHG) monthly climatology, topography + land-sea mask, aerosols (optional) 
        - freely evolving: full atmos, full land
        - cases
            - historical
            - amip4xCO2 - rapid cloud response (part of CFMIP; Bony et al 11)
            - amip4K - dynamical response of atmos on cloud feedback (part of CFMIP; Bony et al 11)
            - related individual experiments
                - more targetted prescription of different land components ([Ackerley et al 18](https://gmd.copernicus.org/articles/11/3865/2018/))
        - detailed [PCMDI requirements](https://pcmdi.llnl.gov/mips/amip/requirements.html)
        - implement both bucket and full soil-vegetation-atmosphere transfer scheme (SVAT) - see Bonan, [AMIP II report](https://pcmdi.github.io/mips/amip/DIAGSUBS/sp12.html#Henderson-Sellers%20et%20al.%201995)


- Prescribe SSTs and sea ice concentrations using observations (e.g. 1979-2014 in CMIP6), irradiance, orbital, CO2 (1 value), ozone (and other GHG) monthly climatology, topography + land-sea mask, aerosols (optional) 
- freely evolving: full atmos, full land
- cases
    - historical
    - amip4xCO2 - rapid cloud response (part of CFMIP; Bony et al 11)
    - amip4K - dynamical response of atmos on cloud feedback (part of CFMIP; Bony et al 11)
    - related individual experiments
        - more targetted prescription of different land components ([Ackerley et al 18](https://gmd.copernicus.org/articles/11/3865/2018/))
- detailed [PCMDI requirements](https://pcmdi.llnl.gov/mips/amip/requirements.html)
- chemistry - different levels of complexty exist - CO2 (2.5 W/m2; well mixed), CH4 (1W/m2; not as well mixed), other aerosols, 0_3 (but only significant in strato)

## [Ackerley et al 18](https://gmd.copernicus.org/articles/11/3865/2018/) study
- ACCESS model
- 3.75 lon x 2.5 lat, 38 vert
- Params: ppt, cloud, convection, rad, BL processes, aerosols
- land+soil: Met Office Surface Exchange Scheme (MOSES; Cox et al., 1999; Essery et al., 2001)
- when land active: T_land adjusted implicitly after figuring out changes in moisture availability, and then again for snow melt
- when land prescribed - just use a value like with fixes SSTs
- Subgrid-scale surface heterogeneity-split the grid box into smaller “tiles” of which there are nine different types (for each calculating T, OLR, SI, LH, SH).

# Other notes
- [aogcms](https://www.narccap.ucar.edu/about/aogcms.html) - provides BCs to regional CMs

# Soil processes
- land process modelling intro: https://pcmdi.github.io/mips/amip/DIAGSUBS/sp12.html

- [Henderson-sellers 95](https://journals.ametsoc.org/view/journals/clim/8/5/1520-0442_1995_008_1043_asemtl_2_0_co_2.xml?tab_body=pdf)
- 

## Soil hydrology
- Bucket formulations
    - [Manabe 69](https://journals.ametsoc.org/view/journals/mwre/97/11/1520-0493_1969_097_0775_catoc_2_3_co_2.xml)
        - The moisture field capacity is a spatially uniform 0.15 m, with surface runoff occurring if the predicted soil moisture exceeds this value. Snowmelt contributes to soil moisture, but if snow covers a grid box completely, the permeability of the soil to falling liquid precipitation becomes zero. For partial snow cover, the permeability decreases proportional to increasing snow fraction.
    - more complex buckets (e.g. separate for soil, roots and groundwater - Milly and Shmakin (2002) - bonan)   exist
    - other refs - potential for initial AMIP run: 
        - https://pcmdi.github.io/projects/modeldoc/cmip1/ccsr_tbls.html
        - https://www.osti.gov/servlets/purl/928178/
    - 


