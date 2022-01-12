# **TempestRemap plugin**

# TempestRemap_jll
- normally would have to install TempestRemap and configure, the generate meshes / overlap meshes / weights
	- e.g.: `./GenerateCSMesh --res 6 --alt --file gravitySam.000000.3d.cubedSphere_6.netcdf --out_format Netcdf4`
	- TempestRemap_jll simplifies this, so no need to download / configure the original tempest remap 
- need Julia 1.7 
- [this PR](https://github.com/JuliaPackaging/Yggdrasil/pull/4174)
- [info on JLL pkgs](https://docs.binarybuilder.org/stable/jll/)
```
julia> using TempestRemap_jll

julia> run(`$(TempestRemap_jll.GenerateCSMesh_exe()) --res 6 --file cubedsphere_6.netcdf  --out_format Netcdf4`);
Parameters:
  --res <integer> [6] 
  --file <string> ["cubedsphere_6.netcdf"] 
  --out_format <string> ["Netcdf4"] 
  --alt <bool> [false] 
=========================================================
..Generating mesh with resolution [6]
..Writing mesh to file [cubedsphere_6.netcdf] 
Nodes per element
..Block 1 (4 nodes): 216
..Mesh generator exited successfully
=========================================================

shell> ncdump -h cubedsphere_6.netcdf
netcdf cubedsphere_6 {
dimensions:
	len_string = 33 ;
	len_line = 81 ;
	four = 4 ;
	time_step = UNLIMITED ; // (0 currently)
	num_dim = 3 ;
	num_nodes = 218 ;
	num_elem = 216 ;
	num_qa_rec = 1 ;
	num_el_blk = 1 ;
	num_el_in_blk1 = 216 ;
	num_nod_per_el1 = 4 ;
	num_att_in_blk1 = 1 ;
variables:
	double time_whole(time_step) ;
	char qa_records(num_qa_rec, four, len_string) ;
	char coor_names(num_dim, len_string) ;
	char eb_names(num_el_blk, len_string) ;
	int eb_status(num_el_blk) ;
	int eb_prop1(num_el_blk) ;
		eb_prop1:name = "ID" ;
	double attrib1(num_el_in_blk1, num_att_in_blk1) ;
	int connect1(num_el_in_blk1, num_nod_per_el1) ;
		connect1:elem_type = "SHELL4" ;
	int global_id1(num_el_in_blk1) ;
	int edge_type1(num_el_in_blk1, num_nod_per_el1) ;
	double coord(num_dim, num_nodes) ;

// global attributes:
		:api_version = 5.f ;
		:version = 5.f ;
		:floating_point_word_size = 8 ;
		:file_size = 0 ;
		:title = "tempest(cubedsphere_6.netcdf) 01/06/2022: 08:44:30" ;
}
```

# TempestRemap features
## Mesh generation
- Source/target meshes
	- cubed sphere
	- lat-lon
	- geodesic
- Overlap mesh
	- collects intersections of the two meshes

## Offline Linear Weight Map Generation
- Types 
	- FV
	- continuous FE
	- discontinuous FE

## Offline map application
- use to check with ClimeCore application

## Alternatives
- [Conduit](https://llnl-conduit.readthedocs.io/en/latest/blueprint_mesh.html) - JSON + binary
- OASIS regridding - but quite clunky

# Refs 
- [TempestRemap docs](https://github.com/ClimateGlobalChange/tempestremap)
- [Ullrich & Taylor 15](https://journals.ametsoc.org/view/journals/mwre/143/6/mwr-d-14-00343.1.xml )
- [Ullrich et al. 16](https://journals.ametsoc.org/view/journals/mwre/144/4/mwr-d-15-0301.1.xml)

# Questions
Q: in CG cartesian do not need to do the projection to conserative and consistent space because integration is exact, right? What about spherical?
Q: SE{1} in CG caused instability, correct?