#=
## Logging
When Julia 1.10+ is used interactively, stacktraces contain reduced type information to make them shorter.
Given that ClimaCore objects are heavily parametrized, non-abbreviated stacktraces are hard to read,
so we force abbreviated stacktraces even in non-interactive runs.
(See also `Base.type_limited_string_from_context()`)
=#

redirect_stderr(IOContext(stderr, :stacktrace_types_limited => Ref(true)))

#=
## Package Loading
Import all packages needed to run coupled AMIP simulations.
This file can be included from the REPL to load everything needed
to set up and run a simulation interactively.
=#

using ClimaCoupler

# The Makie plotting stack (CairoMakie, GeoMakie, ...) is intentionally NOT loaded
# here: it is only needed for postprocessing and its load/compile time is large.
# `run_simulation.jl` includes `../load_plotting.jl` right before `postprocess` to
# trigger `ClimaCouplerMakieExt` once the simulation is done. Include that file
# manually if you want plotting in an interactive session.

# Trigger ClimaCouplerClimaLandExt
import ClimaLand

# Trigger ClimaCouplerClimaAtmosExt
import ClimaAtmos

import ClimaComms
ClimaComms.@import_required_backends
