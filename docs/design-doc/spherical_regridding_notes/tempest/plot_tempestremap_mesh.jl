# ClimaCore Meshes

using TempestRemap_jll
using NCDatasets

file_name = "/Users/lenkanovak/ClimaCoupler.jl/cubedsphere_6.netcdf"
ds = NCDataset(file_name, "r");

coord = ds["coord"][:]


x = coord[:,1]
y = coord[:,2]
z = coord[:,3]

i= 218 
Makie.linesegments(x[1:i],y[1:i],z[1:i], overdraw = true)
Makie.scatter(y[1:i],z[1:i], overdraw = true) # 3d scatter doesn't work - known bug - https://discourse.julialang.org/t/ploints-not-appearing-in-3d-scatter/64622

# Options
# 1. tempest generates: source, target, overlap meshes + weights file 

