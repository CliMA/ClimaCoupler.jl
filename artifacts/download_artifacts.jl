include(joinpath(@__DIR__, "artifact_funcs.jl"))

# Trigger download if data doesn't exist locally
function trigger_download()
    @info "sst dataset path: `$(sst_dataset_path())`"
    @info "sic dataset path: `$(sic_dataset_path())`"
    @info "co2 dataset path: `$(co2_dataset_path())`"
    @info "mask dataset path: `$(mask_dataset_path())`"
    return nothing
end
trigger_download()
