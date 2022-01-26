
# demo of FV > FV regridding using CC meshes and TempestRemap test data

import ClimaCore
using ClimaCore: Geometry, Meshes, Domains, Topologies, Spaces
using NCDatasets
using TempestRemap_jll
using Test
using ClimaCoreTempestRemap


OUTPUT_DIR = mkdir("output_fv")

# input mesh
ne_i = 20
R = 1.0 # unit sphere 
domain = ClimaCore.Domains.SphereDomain(R)
mesh = ClimaCore.Meshes.EquiangularCubedSphere(domain, ne_i) # ne×ne×6 
grid_topology = ClimaCore.Topologies.Topology2D(mesh)
nc_name_in = joinpath(OUTPUT_DIR, "test_in.nc")
#write_exodus(nc_name_in, grid_topology)
run(`$(TempestRemap_jll.GenerateCSMesh_exe()) --res $ne_i --alt --file $nc_name_in`,)

# output mesh
ne_o = 10
R = 1.0 
domain = ClimaCore.Domains.SphereDomain(R)
mesh = ClimaCore.Meshes.EquiangularCubedSphere(domain, ne_o) # ne×ne×6 
grid_topology = ClimaCore.Topologies.Topology2D(mesh)
nc_name_out = joinpath(OUTPUT_DIR, "test_out.nc")
#write_exodus(nc_name_out, grid_topology)
run(`$(TempestRemap_jll.GenerateCSMesh_exe()) --res $ne_o --alt --file $nc_name_out`,)


# overlap mesh
nc_name_ol = joinpath(OUTPUT_DIR, "test_ol.g")
run(`$(TempestRemap_jll.GenerateOverlapMesh_exe()) --a $nc_name_in --b $nc_name_out --out $nc_name_ol`)

# map weights 
nc_name_wgt = joinpath(OUTPUT_DIR, "test_wgt.g")
run(`$(TempestRemap_jll.GenerateOfflineMap_exe()) --in_mesh $nc_name_in --out_mesh $nc_name_out --ov_mesh $nc_name_ol --in_np 1 --out_map $nc_name_wgt`) # FV > FV, parallelized sparse matrix multiply (SparseMatrix.h)

# generate fake input data (replace with CC variable; NB this requires write_exodus_identical() above)
nc_name_data_in = joinpath(OUTPUT_DIR, "Psi_in.nc")
run(`$(GenerateTestData_exe()) --mesh $nc_name_in --test 1 --out $nc_name_data_in`) # var: Psi


# apply map (this to be done by CC at each timestep)
nc_name_data_out = joinpath(OUTPUT_DIR, "Psi_out.nc")
run(`$(TempestRemap_jll.ApplyOfflineMap_exe()) --map $nc_name_wgt --var Psi --in_data $nc_name_data_in --out_data $nc_name_data_out`)

# try FEM


using Plots             
function plot_flatmesh(Psi,nelem)
    plots = []
    tiltes = ["Eq1" "Eq2" "Eq3" "Eq4" "Po1" "Po2"]
    for f in collect(1:1:6)
        Psi_reshape = reshape(Psi[(f-1)*nelem^2+1:f*nelem^2],(nelem,nelem))

        push!(plots, contourf(Psi_reshape))
    end
    plot(plots..., layout = (6), title = tiltes )
end

# viz
ds_indata = NCDataset(nc_name_data_in,"r")
Psi_in = ds_indata["Psi"][:]
close(ds_indata)
ds_inmesh = NCDataset(nc_name_in,"r")
connect1_in = ds_inmesh["connect1"][:]
coord_in = ds_inmesh["coord"][:]
plot_flatmesh(Psi_in,ne_i)
png(joinpath(OUTPUT_DIR,"in.png"))

#map(i -> Psi_in[i], connect1_in)

## get regriddded field and its mapping
ds_outdata = NCDataset(nc_name_data_out,"r")
Psi_out = ds_outdata["Psi"][:] # on unique gll nodes
close(ds_outdata)
ds_outmesh = NCDataset(nc_name_out,"r")
connect1_out = ds_outmesh["connect1"][:]
coord_out = ds_outmesh["coord"][:]
plot_flatmesh(Psi_out,ne_o)
png(joinpath(OUTPUT_DIR,"out.png"))
