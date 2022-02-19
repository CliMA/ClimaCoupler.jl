# demo of FE > FE regridding using CC meshes and Fields

#import ClimaCore
import Pkg; Pkg.add(url="https://github.com/CliMA/ClimaCore.jl",rev="main")

using ClimaCore
using ClimaCore: Geometry, Meshes, Domains, Topologies, Spaces
using NCDatasets
using TempestRemap_jll
using Test

import Pkg; Pkg.add(url="https://github.com/CliMA/ClimaCore.jl.git", subdir="lib/ClimaCoreTempestRemap")
using ClimaCoreTempestRemap

write_exodus = ClimaCoreTempestRemap.write_exodus

nq = 3 # polynomial order (here using same nq for source and target data)

# setup output dir
OUTPUT_DIR = "output_fv"
isdir(OUTPUT_DIR) ? nothing : mkdir(OUTPUT_DIR)

# source grid params
ne_i = 20 # #elements
R = 1.0 # unit sphere 
no_unique_mesh_nodes = (ne_i^2 * 6 * 4 - 8*3) / 4 + 8 
no_unique_gll_nodes = (ne_i^2 * 6 * (nq-1)*(nq-1))  - (8*3 ) / 4 + 8
num_elem = ne_i^2 * 6 

# target grid params
ne_o = 5
R = 1.0 
no_unique_mesh_nodes_o = (ne_o^2 * 6 * 4 - 8*3) / 4 + 8 # = gll nodes if np = 2
no_unique_gll_nodes_o = (ne_o^2 * 6 * (nq-1)*(nq-1))  - (8*3 ) / 4 + 8
num_elem_o = Int(ne_o^2 * 6)

# construct source mesh
domain = ClimaCore.Domains.SphereDomain(R)
mesh_in = ClimaCore.Meshes.EquiangularCubedSphere(domain, ne_i) 
grid_topology_in = ClimaCore.Topologies.Topology2D(mesh_in)
nc_name_in = joinpath(OUTPUT_DIR, "test_in.nc")
write_exodus(nc_name_in, grid_topology_in)
#run(`$(TempestRemap_jll.GenerateCSMesh_exe()) --res $ne_i --alt --file $nc_name_in`,) # if want to generate the same mesh (with different ordering) with TempestRemap

# construct target mesh
domain = ClimaCore.Domains.SphereDomain(R)
mesh_out = ClimaCore.Meshes.EquiangularCubedSphere(domain, ne_o) 
grid_topology_out = ClimaCore.Topologies.Topology2D(mesh_out)
nc_name_out = joinpath(OUTPUT_DIR, "test_out.nc")
write_exodus(nc_name_out, grid_topology_out)
#run(`$(TempestRemap_jll.GenerateCSMesh_exe()) --res $ne_o --alt --file $nc_name_out`,)

# overlap mesh
nc_name_ol = joinpath(OUTPUT_DIR, "test_ol.g")
run(`$(TempestRemap_jll.GenerateOverlapMesh_exe()) --a $nc_name_in --b $nc_name_out --out $nc_name_ol`)

# map weights - TempestRemap splits mesh into GLL nodes (so these are not passed from CC)
nc_name_wgt = joinpath(OUTPUT_DIR, "test_wgt.g")
run(`$(TempestRemap_jll.GenerateOfflineMap_exe()) --in_mesh $nc_name_in --out_mesh $nc_name_out --ov_mesh $nc_name_ol --in_type cgll --out_type cgll --in_np $nq --out_np $nq --out_map $nc_name_wgt`) # GLL > GLL - crashing 

# generate fake input data in exodus format (NB: data overwritten below using CC Field data)
nc_name_data_in = joinpath(OUTPUT_DIR, "Psi_in.nc")
run(`$(GenerateTestData_exe()) --mesh $nc_name_in --test 1 --out $nc_name_data_in --gllint --np $nq`) # default var name = "Psi" ; "a": src Dims, "b": dst dim

# Reset with u values from a CC Field
FT = Float64
quad = Spaces.Quadratures.GLL{nq}()
space = ClimaCore.Spaces.SpectralElementSpace2D(grid_topology_in, quad) #float_type(::AbstractPoint{FT}) where {FT} = FT

coords = ClimaCore.Fields.coordinate_field(space)
u = map(coords) do coord
    ϕ = coord.lat
    uu = cosd(ϕ)
    ClimaCore.Geometry.UVVector(uu, uu)
end

u = u.components.data.:1
u_vals=getfield(u, :values) #IJFH : nq,nq,1,1,nelem

# loop through CC GLL points and generate GLL connectivity matrix
function get_cc_gll_connect(ne_i, nq)

    ntot = ne_i*(nq-1)
    face1 = collect(-ntot/2:1:ntot/2) * collect(-ntot/2:1:ntot/2)'
    face2 = collect(ntot/2:-1:-ntot/2)*collect(-ntot/2:1:ntot/2)'
    
    face3 = collect(ntot/2:-1:-ntot/2) * collect(ntot/2:-1:-ntot/2)'
    face4 = collect(ntot/2:-1:-ntot/2)* collect(ntot/2:-1:-ntot/2)'
    
    face5 = collect(ntot/2:-1:-ntot/2) * collect(-ntot/2:1:ntot/2)' 
    face6 = collect(-ntot/2:1:ntot/2) * collect(-ntot/2:1:ntot/2)' 
    
    face1_c = [ones(ntot+1)*ntot/2,  collect(-ntot/2:1:ntot/2), collect(-ntot/2:1:ntot/2)]
    face2_c = [collect(ntot/2:-1:-ntot/2),  ones(ntot+1)*ntot/2 , collect(-ntot/2:1:ntot/2)]
    
    face3_c = [collect(ntot/2:-1:-ntot/2), collect(ntot/2:-1:-ntot/2), ones(ntot+1)*ntot/2]
    face4_c = [-ones(ntot+1)*ntot/2, collect(ntot/2:-1:-ntot/2), collect(ntot/2:-1:-ntot/2)]
    
    face5_c = [collect(-ntot/2:1:ntot/2), -ones(ntot+1)*ntot/2, collect(ntot/2:-1:-ntot/2) ]
    face6_c = [collect(-ntot/2:1:ntot/2), collect(-ntot/2:1:ntot/2),  -ones(ntot+1)*ntot/2 ]
    
    six_faces_dc = [face1_c,face2_c,face3_c,face4_c,face5_c,face6_c]
    six_faces = [face1,face2,face3,face4,face5,face6]

    connect = zeros(nq, nq, 6 * ne_i * ne_i )

    elem_ct = 0
    coords = []
    for f in collect(1:1:6)
        for e_y in collect(1:1:ne_i)
            for e_x in collect(1:1:ne_i)
                for nq_y in collect(1:1:nq)
                    for nq_x in collect(1:1:nq)
                        if f in [1]
                            x = six_faces_dc[f][1][1]
                            y = six_faces_dc[f][2][(e_x-1)*(nq) + nq_x - (e_x-1)]
                            z = six_faces_dc[f][3][(e_y-1)*(nq) + nq_y - (e_y-1)]
                            if (x,y,z) ∉ coords
                                push!(coords, (x,y,z))
                                elem_ct +=1
                                connect[nq_x, nq_y, (f-1)*ne_i^2 + (e_y-1)*ne_i + e_x]  = elem_ct
                            else
                                ui = findfirst(isequal((x,y,z)), coords)
                                connect[nq_x, nq_y, (f-1)*ne_i^2 + (e_y-1)*ne_i + e_x]  = ui
                            end
                        end
                        if f in [4]
                            x = six_faces_dc[f][1][1]
                            y = six_faces_dc[f][2][(e_y-1)*(nq) + nq_y - (e_y-1)]
                            z = six_faces_dc[f][3][(e_x-1)*(nq) + nq_x - (e_x-1)]
                            if (x,y,z) ∉ coords
                                push!(coords, (x,y,z))
                                elem_ct +=1
                                connect[nq_x, nq_y, (f-1)*ne_i^2 + (e_y-1)*ne_i + e_x]  = elem_ct
                            else
                                ui = findfirst(isequal((x,y,z)), coords)
                                connect[nq_x, nq_y, (f-1)*ne_i^2 + (e_y-1)*ne_i + e_x]  = ui

                            end
                        end
                        if f in [2]
                            x = six_faces_dc[f][1][(e_x-1)*(nq) + nq_x - (e_x-1)]
                            y = six_faces_dc[f][2][1]
                            z = six_faces_dc[f][3][(e_y-1)*(nq) + nq_y - (e_y-1)]
                            if (x,y,z) ∉ coords
                                push!(coords, (x,y,z))
                                elem_ct +=1
                                connect[nq_x, nq_y, (f-1)*ne_i^2 + (e_y-1)*ne_i + e_x]  = elem_ct
                            else
                                ui = findfirst(isequal((x,y,z)), coords)
                                connect[nq_x, nq_y, (f-1)*ne_i^2 + (e_y-1)*ne_i + e_x]  = ui
                                println(ui)
                            end     
                        end
                        if f in [5]
                            x = six_faces_dc[f][1][(e_y-1)*(nq) + nq_y - (e_y-1)]
                            y = six_faces_dc[f][2][1]
                            z = six_faces_dc[f][3][(e_x-1)*(nq) + nq_x - (e_x-1)]
                            if (x,y,z) ∉ coords
                                push!(coords, (x,y,z))
                                elem_ct +=1
                                connect[nq_x, nq_y, (f-1)*ne_i^2 + (e_y-1)*ne_i + e_x]  = elem_ct
                            else
                                ui = findfirst(isequal((x,y,z)), coords)
                                connect[nq_x, nq_y, (f-1)*ne_i^2 + (e_y-1)*ne_i + e_x]  = ui

                            end     
                        end
                        if f in [3]
                            x = six_faces_dc[f][1][(e_x-1)*(nq) + nq_x - (e_x-1)]
                            y = six_faces_dc[f][2][(e_y-1)*(nq) + nq_y - (e_y-1)]
                            z = six_faces_dc[f][3][1]
                            if (x,y,z) ∉ coords
                                push!(coords, (x,y,z))
                                elem_ct +=1
                                connect[nq_x, nq_y, (f-1)*ne_i^2 + (e_y-1)*ne_i + e_x]  = elem_ct
                            else
                                ui = findfirst(isequal((x,y,z)), coords)
                                connect[nq_x, nq_y, (f-1)*ne_i^2 + (e_y-1)*ne_i + e_x]  = ui
                            end              
                        end
                        if f in [6]
                            x = six_faces_dc[f][1][(e_y-1)*(nq) + nq_y - (e_y-1)]
                            y = six_faces_dc[f][2][(e_x-1)*(nq) + nq_x - (e_x-1)]
                            z = six_faces_dc[f][3][1]
                            if (x,y,z) ∉ coords
                                push!(coords, (x,y,z))
                                elem_ct +=1
                                connect[nq_x, nq_y, (f-1)*ne_i^2 + (e_y-1)*ne_i + e_x]  = elem_ct
                            else
                                ui = findfirst(isequal((x,y,z)), coords)
                                connect[nq_x, nq_y, (f-1)*ne_i^2 + (e_y-1)*ne_i + e_x]  = ui
                            end              
                        end
                    end
                end
            end
        end
    end
    return (coords, elem_ct, connect)
end

# generate connectivity for source gll points
coords, elem_ct, conn = get_cc_gll_connect(ne_i, nq)

unwrap_cc_coord2(coord) = [coord[1] , coord[2], coord[3]]

# convert to unique indices as expected by TR
u_vals_uq = map(x ->  getindex(parent(u_vals)[:,:,1,:],unwrap_cc_coord2(findfirst(isequal(x), conn))...), collect(1:1:elem_ct)) #IJFH : nq,nq,1,nelem

ds_indata = NCDataset(nc_name_data_in,"a")
ds_indata["Psi"][:] = u_vals_uq[:] 
close(ds_indata)

# load map
ds_wt = NCDataset("output_fv/test_wgt.g","r")
S = ds_wt["S"][:]
row = ds_wt["row"][:]
col = ds_wt["col"][:]
close(ds_wt)

# apply map, S
u_vals_uq_out = zeros(Int(no_unique_gll_nodes_o))
for (i, val) in enumerate(S)
    u_vals_uq_out[row[i]] += S[i]* u_vals_uq[col[i]]
end

# generate connectivity for target gll points
coords_o, elem_ct_o, conn_o = get_cc_gll_connect(ne_o, nq)

u_vals_out = Array{Float64}(undef, nq, nq, 1,num_elem_o) 
for e in collect(1:1:num_elem_o)
    for nq_y in collect(1:1:nq) 
        for nq_x in collect(1:1:nq) 
            uv_conn = conn_o[nq_y, nq_x, e]
            u_vals_out[nq_x,nq_y,1,e] = u_vals_uq_out[Int(uv_conn)]
        end
    end
end

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

plot_flatmesh(parent(u_vals)[1,1,1,:],ne_i)
png(joinpath(OUTPUT_DIR,"in.png"))

plot_flatmesh(parent(u_vals_out)[1,1,1,:],ne_o)
png(joinpath(OUTPUT_DIR,"out.png"))
# NB: plots will be offset if only plotting the 1st GLL node of each elem



