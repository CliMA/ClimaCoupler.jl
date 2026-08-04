"""
    Models

Module containing simple component models implemented within ClimaCoupler.jl
and intended to be run as component models in a ClimaCoupler simulation.
This includes slab ocean, prescribed ocean, prescribed sea ice, and
Eisenman-Zhang sea ice models.
"""
module Models

include("Models/slab_ocean.jl")
include("Models/prescr_ocean.jl")
include("Models/prescr_seaice.jl")
include("Models/eisenman_seaice.jl")

end # module Models
