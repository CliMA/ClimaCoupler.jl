# **Test cases**
- See [here](/Users/lenka/ClimaCoupler.jl/experiments/ClimaCore/README.md) for current test case implementations. We broadly follow a two-way advancement in complexity
:
![](figures/testcase_hierarchy.png)

# Milestones
- SCM
    - RCE atmos + land
- cartesian
    - LES + LES
    - LES + Land
- spherical
    - aquaplanet
    - bucket 3-way coupling
    - AMIP simulation
        - Prescribe SSTs and sea ice concentrations using observations (e.g. 1979-2014 in CMIP6), irradiance, orbital, CO2 (1 value), ozone (and other GHG) monthly climatology, topography + land-sea mask, aerosols (optional) 
        - freely evolving: full atmos, full land
        - cases
            - historical
            - amip4xCO2 - rapid cloud response (part of CFMIP; Bony et al 11)
            - amip4K - dynamical response of atmos on cloud feedback (part of CFMIP; Bony et al 11)
            - related individual experiments
                - more targetted prescription of different land components ([Ackerley et al 18](https://gmd.copernicus.org/articles/11/3865/2018/))
        - detailed [PCMDI requirements](https://pcmdi.llnl.gov/mips/amip/requirements.html)
- 

# Refs:
- [2021 single column test case results on Overleaf](https://www.overleaf.com/read/bgfmhgtncpws)