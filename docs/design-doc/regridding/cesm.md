# **CMEPS (CESM coupler)**

- [ESMF](https://earthsystemmodeling.org/) (Earth System Modelling Framework) - layer on top of component models, supplying also links DA and other tools
- provides: interfaces, regridding, parallel comm

- NUOPC (National Unified Operational Prediciton Capability) - layer on top of ESMF for conventions etc

# Couplers 
- MCT coupler (traditionally) - also used in OASIS
- CMEPS (Community Mediator for Earth prediction Systems) 
    - e.g. has NUOPC CESM driver
    - will be default for CESM 2.3

# CMEPS
- [CMEPS tutorial](https://anl.app.box.com/s/bss08sjmndvx3zk53btv58dmgyxjcgsh?page=1
)
- [CMEPS GH](https://github.com/ESCOMP/CMEPS)
- structure
    - ESMF gridding comp
        - NUOPC_driver
        - NUOPC_model
        - NUOPC_mediator 
            - gen remap weights
            - exchange grid option for flux calculation
    - ESMF coupling comp
        - NUOPC_connector
            - transfer data betwen comp and mediator
- [regridding](https://earthsystemmodeling.org/regrid/)
    - various methods supported, modular interface
- no hooks required in models requires
- model handles are in the mediator

# CIME (Common Infrastructure for Modeling the Earth)
-  a Case Control System for configuring, compiling and executing Earth system models  data and stub model components, a driver and associated tools and libraries. 
- [CIME docs](http://esmci.github.io/cime/versions/master/html/index.html)
- [CIME GH](https://github.com/ESMCI/cime)
- a "fork" of the cime repository that has the development version of the nuopc CMEPS driver and mediator 
