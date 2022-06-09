using Pkg.Artifacts
using ArtifactWrappers
sst_dataset = ArtifactWrapper(
    @__DIR__,
    isempty(get(ENV, "CI", "")),
    "sst",
    ArtifactFile[ArtifactFile(
        url = "https://caltech.box.com/shared/static/8gd3wjq2dbrnzuv8pd5ww3a54h0kkcz8.nc",
        filename = "sst.nc",
    ),],
)
sst_dataset_path = get_data_folder(sst_dataset)
sst_data = joinpath(sst_dataset_path, "sst.nc")

sic_dataset = ArtifactWrapper(
    @__DIR__,
    isempty(get(ENV, "CI", "")),
    "sic",
    ArtifactFile[ArtifactFile(
        url = "https://caltech.box.com/shared/static/n0omgqkmnwpr9gylhixew8ywb4psgvj4.nc",
        filename = "sic.nc",
    ),],
)
sic_dataset_path = get_data_folder(sic_dataset)
sic_data = joinpath(sic_dataset_path, "sic.nc")

mask_dataset = ArtifactWrapper(
    @__DIR__,
    isempty(get(ENV, "CI", "")),
    "mask",
    ArtifactFile[ArtifactFile(
        url = "https://caltech.box.com/shared/static/vubmq84nhvbgdqayezguf3i1w6nqtwvu.ncc",
        filename = "seamask.nc",
    ),],
)
mask_dataset_path = get_data_folder(mask_dataset)
mask_data = joinpath(mask_dataset_path, "seamask.nc")
