"""
    ClimaCouplerModelsExt

Module containing simple component models implemented within ClimaCoupler.jl
and intended to be run as component models in a ClimaCoupler simulation.
This includes slab ocean, prescribed ocean, and prescribed sea ice models.
"""
module ClimaCouplerModelsExt

import ClimaCoupler

# Include the model files
include("ClimaCouplerModelsExt/slab_ocean.jl")
include("ClimaCouplerModelsExt/prescr_ocean.jl")
include("ClimaCouplerModelsExt/prescr_seaice.jl")

end # module ClimaCouplerModelsExt
