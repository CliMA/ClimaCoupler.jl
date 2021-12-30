# **Test cases**
- See [here](/Users/lenka/ClimaCoupler.jl/experiments/ClimaCore/README.md) for current test case implementations. We broadly follow a two-way advancement in complexity
:
![](figures/testcase_hierarchy.png)

# Milestones
- SCM
    - RCE atmos + land + bucket
- cartesian box
    - Atmos LES + Land
    - Atmos LES + Ocean LES
    - Atmos LES + Ocean LES + Land (full sea breeze) + bucket / SVAT
- sphere
    - aquaplanet with slab
    - small aquaplanet (convection resolving) with oceananigans
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
        - implement both bucket and full soil-vegetation-atmosphere transfer scheme (SVAT) - see Bonan, [AMIP II report](https://pcmdi.github.io/mips/amip/DIAGSUBS/sp12.html#Henderson-Sellers%20et%20al.%201995)
    - CMIP simulation
        - All Atmos, Land and Ocean components interactive
        - sea ice model
        - chemistry model

# Plan / dependencies
- SCM RCE + cartesian box implementations
    - `ClimaAtmos` and `CliamLand` interfaces + compatibility with `Oceananigans` (`ClimaSimulations`) [Akshay/Greg/Toby?]
    - moisture / TD implementation [Toby/Akshay/Charlie?]
    - regridding (almost done) [Dec/Jan]
    - SVAT [Kat]
- AMIP spherical implementations 
    - `ClimaAtmos` and `CliamLand` interfaces + compatibility with `Oceananigans` (`ClimaSimulations`)
    - moisture / TD implementation [Toby/Charlie?]
    - regridding for sphere [Jan/Feb]
    - implemetation of topography [Akshay?]
    - implementation of land/sea mask [Akshay?]
    - SVAT [Kat]


# Refs
- [2021 single column test case results on Overleaf](https://www.overleaf.com/read/bgfmhgtncpws) 
- [Henderson-sellers 95](https://journals.ametsoc.org/view/journals/clim/8/5/1520-0442_1995_008_1043_asemtl_2_0_co_2.xml?tab_body=pdf)
- [Manabe 69](https://journals.ametsoc.org/view/journals/mwre/97/11/1520-0493_1969_097_0775_catoc_2_3_co_2.xml)
