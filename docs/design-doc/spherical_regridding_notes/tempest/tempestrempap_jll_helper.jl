
run(`$(TempestRemap_jll.GenerateCSMesh_exe()) --res 64 --alt --file cubedsphere_64.g`)
run(`$(TempestRemap_jll.GenerateRLLMesh_exe()) --lon 256 --lat 128 --file latlonsphere_256.g`)
run(`$(TempestRemap_jll.GenerateOverlapMesh_exe()) --a cubedsphere_64.g --b latlonsphere_256.g --out overlapsphere.g`)
# FV - FV:
run(`$(TempestRemap_jll.GenerateOfflineMap_exe()) --in_mesh cubedsphere_64.g --out_mesh latlonsphere_256.g --ov_mesh overlapsphere.g --in_np 1 --out_map weights.nc`)
