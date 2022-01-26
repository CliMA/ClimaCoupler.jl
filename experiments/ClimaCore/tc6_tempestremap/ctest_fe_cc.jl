
# demo of FE > FE regridding using CC meshes and Fields

import ClimaCore
using ClimaCore: Geometry, Meshes, Domains, Topologies, Spaces
using NCDatasets
using TempestRemap_jll
using Test
using ClimaCoreTempestRemap

nq = 3

#OUTPUT_DIR = mkdir("output_fe_u_ccidx")

# input mesh
ne_i = 20
R = 1.0 # unit sphere 
domain = ClimaCore.Domains.SphereDomain(R)
mesh_in = ClimaCore.Meshes.EquiangularCubedSphere(domain, ne_i) # ne×ne×6 
grid_topology_in = ClimaCore.Topologies.Topology2D(mesh_in)
nc_name_in = joinpath(OUTPUT_DIR, "test_in.nc")
write_exodus(nc_name_in, grid_topology_in)
#run(`$(TempestRemap_jll.GenerateCSMesh_exe()) --res $ne_i --alt --file $nc_name_in`,)

# output mesh
ne_o = 5
R = 1.0 
domain = ClimaCore.Domains.SphereDomain(R)
mesh_out = ClimaCore.Meshes.EquiangularCubedSphere(domain, ne_o) # ne×ne×6 
grid_topology_out = ClimaCore.Topologies.Topology2D(mesh_out)
nc_name_out = joinpath(OUTPUT_DIR, "test_out.nc")
write_exodus(nc_name_out, grid_topology_out)
#run(`$(TempestRemap_jll.GenerateCSMesh_exe()) --res $ne_o --alt --file $nc_name_out`,)


# overlap mesh
nc_name_ol = joinpath(OUTPUT_DIR, "test_ol.g")
run(`$(TempestRemap_jll.GenerateOverlapMesh_exe()) --a $nc_name_in --b $nc_name_out --out $nc_name_ol`)

# map weights 
nc_name_wgt = joinpath(OUTPUT_DIR, "test_wgt.g")
run(`$(TempestRemap_jll.GenerateOfflineMap_exe()) --in_mesh $nc_name_in --out_mesh $nc_name_out --ov_mesh $nc_name_ol --in_type cgll --out_type cgll --in_np $nq --out_np $nq --out_map $nc_name_wgt`) # GLL > GLL - crashing 

# generate fake input data (replace with CC variable; NB this requires write_exodus_identical() above)
nc_name_data_in = joinpath(OUTPUT_DIR, "Psi_in.nc")
run(`$(GenerateTestData_exe()) --mesh $nc_name_in --test 1 --out $nc_name_data_in --gllint --np $nq`) # var: Psi

### Reset with u values
# try FEM
no_unique_mesh_nodes = (ne_i^2 * 6 * 4 - 8*3) / 4 + 8 # = gll nodes if np = 2
no_unique_gll_nodes = (ne_i^2 * 6 * (nq-1)*(nq-1))  - (8*3 ) / 4 + 8
num_elem = ne_i^2 * 6 
#a > src Dims
#b > dst dim

FT = Float64
quad = Spaces.Quadratures.GLL{nq}()
space = ClimaCore.Spaces.SpectralElementSpace2D(grid_topology_in, quad) #float_type(::AbstractPoint{FT}) where {FT} = FT

coords = ClimaCore.Fields.coordinate_field(space)
u = map(coords) do coord
    # u0 = 20.0
    # α0 = 45.0
    # ϕ = coord.lat
    # λ = coord.long

    # uu = u0 * (cosd(α0) * cosd(ϕ) + sind(α0) * cosd(λ) * sind(ϕ))
    # uv = -u0 * sind(α0) * sind(λ)
    ϕ = coord.lat
    uu = cosd(ϕ)
    ClimaCore.Geometry.UVVector(uu, uu)
end

u = u.components.data.:1
u_vals=getfield(u, :values) #IJFH : 4,4,1,1,216

# loop through CC GLL points and generate GLL connectivity matrix
function idx(ne_i, nq)

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
coords, elem_ct, conn = idx(ne_i, nq)

unwrap_cc_coord2(coord) = [coord[1] , coord[2], coord[3]]

# convert to unique indices as expected by TR
u_vals_uq = map(x ->  getindex(parent(u_vals)[:,:,1,:],unwrap_cc_coord2(findfirst(isequal(x), conn))...), collect(1:1:elem_ct))

ds_indata = NCDataset(nc_name_data_in,"a")
ds_indata["Psi"][:] = u_vals_uq[:] 
close(ds_indata)

# apply map (this to be done by CC at each timestep)
nc_name_data_out = joinpath(OUTPUT_DIR, "Psi_out.nc")
run(`$(TempestRemap_jll.ApplyOfflineMap_exe()) --map $nc_name_wgt --var Psi --in_data $nc_name_data_in --out_data $nc_name_data_out`)

# load and convert output
no_unique_mesh_nodes_o = (ne_o^2 * 6 * 4 - 8*3) / 4 + 8 # = gll nodes if np = 2
no_unique_gll_nodes_o = (ne_o^2 * 6 * (nq-1)*(nq-1))  - (8*3 ) / 4 + 8
num_elem_o = Int(ne_o^2 * 6)
ds_outdata = NCDataset(nc_name_data_out,"r")

# generate connectivity for target gll points
coords_o, elem_ct_o, conn_o = idx(ne_o, nq)

u_vals_out = Array{Float64}(undef, nq, nq, 1,num_elem_o) 
for e in collect(1:1:num_elem_o)
    for nq_y in collect(1:1:nq) 
        for nq_x in collect(1:1:nq) 
            uv_conn = conn_o[nq_y, nq_x, e]
            println(uv_conn)
            u_vals_out[nq_x,nq_y,1,e] = ds_outdata["Psi"][:][Int(uv_conn)]
        end
    end
end

close(ds_outdata)

plot_flatmesh(parent(u_vals)[1,1,1,:],ne_i)
png(joinpath(OUTPUT_DIR,"in.png"))

plot_flatmesh(parent(u_vals_out)[1,1,1,:],ne_o)
png(joinpath(OUTPUT_DIR,"out.png"))




