module Leaderboard

import ClimaAnalysis
import Dates
import NCDatasets
import Interpolations
import CairoMakie
import GeoMakie
import ClimaCoupler
import Statistics: mean
import Artifacts

include(joinpath(pkgdir(ClimaCoupler), "artifacts", "artifact_funcs.jl"))
include("data_sources.jl")
include("utils.jl")
include("compare_with_obs.jl")


end
