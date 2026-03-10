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
Import all packages needed to run ClimaEarth coupled simulations.
This file can be included from the REPL to load everything needed
to set up and run a simulation interactively.
=#

using ClimaCoupler

# Trigger ClimaCouplerMakieExt
using Makie, GeoMakie, CairoMakie, ClimaCoreMakie, NCDatasets, Poppler_jll

# Trigger ClimaCouplerCMIPExt
import Oceananigans, ClimaOcean, ClimaSeaIce, KernelAbstractions

# Trigger ClimaCouplerClimaLandExt
import ClimaLand

# Trigger ClimaCouplerClimaAtmosExt
import ClimaAtmos
