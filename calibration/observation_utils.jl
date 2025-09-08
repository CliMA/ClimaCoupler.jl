import ClimaCoupler

# include(joinpath(pkgdir(ClimaCoupler), "experiments/ClimaEarth/user_io/arg_parsing.jl"))
# TODO: This is dumb, since I need to create this CoupledSimulation object

"""
    resampled_lonlat(config_file)

Return a function to resample longitude and latitudes according to the model
grid specified by `config_file`.
"""
function resampled_lonlat(config_file)
    config_dict = get_coupler_config_dict(config_file)
    cs = ClimaCoupler.Interfacer.CoupledSimulation(config_dict)
    center_space = cs.model_sims.atmos_sim.domain.center_space
    (lon_nlevels, lat_nlevels, z_nlevels) =
        ClimaDiagnostics.Writers.default_num_points(center_space)
    longitudes = range(-180, 180, lon_nlevels)
    latitudes = range(-90, 90, lat_nlevels)
    stretch = center_space.grid.vertical_grid.topology.mesh.stretch
    # TODO: Account for stretch for 3D variables and interpolate to pressure?
    z_levels = range(dz_bottom, Spaces.z_max(center_space), z_nlevels)
    return var -> resampled_to(var; lon = longitudes, lat = latitudes)
end
