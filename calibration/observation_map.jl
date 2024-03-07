### Place holder for NCEP data from the ClimaCoupler outputs
#

function observation_map(iteration)
    experiment_id = "amip_coupled"
    config = YAML.load_file(joinpath("experiments", experiment_id, "ekp_config.yml"))
    output_dir = config["output_dir"]
    ensemble_size = config["ensemble_size"]
    model_output = "wa_inst.nc"
    dims = 1
    G_ensemble = Array{Float64}(undef, dims..., ensemble_size)
    for m in 1:ensemble_size
        member_path =
            TOMLInterface.path_to_ensemble_member(output_dir, iteration, m)
        ta = ncread(joinpath(member_path, model_output), "wa")
        G_ensemble[:, m] = process_member_data(ta)
    end
    return G_ensemble
end

function process_member_data(wa; output_variance = false)
    # Cut off first 120 days to get equilibrium, take second level slice
    level_slice = 2
    wa_second_height = wa[3:size(wa)[1], :, :, level_slice]
    # Average over long and latitude
    area_avg_wa_second_height =
        longitudinal_avg(latitudinal_avg(wa_second_height))
    observation = Float64[area_avg_wa_second_height[3]]
    if !(output_variance)
        return observation
    else
        variance = Matrix{Float64}(undef, 1, 1)
        variance[1] = var(area_avg_wa_second_height)
        return (; observation, variance)
    end
end