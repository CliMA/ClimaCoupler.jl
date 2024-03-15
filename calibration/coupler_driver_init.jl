include("../experiments/AMIP/components/atmosphere/climaatmos.jl")
include("../experiments/AMIP/components/land/climaland_bucket.jl")
include("../experiments/AMIP/components/ocean/slab_ocean.jl")
include("../experiments/AMIP/components/ocean/prescr_seaice.jl")
include("../experiments/AMIP/components/ocean/eisenman_seaice.jl")

## helper../experiments/AMIP/s for user-specified IO
include("../experiments/AMIP/user_io/user_diagnostics.jl")
include("../experiments/AMIP/user_io/user_logging.jl")