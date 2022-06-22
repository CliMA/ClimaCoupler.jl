using Pkg.Artifacts
import ArtifactWrappers as AW

import ClimaCoupler

function sst_dataset_path()
    sst_dataset = AW.ArtifactWrapper(
        @__DIR__,
        isempty(get(ENV, "CI", "")),
        "sst",
        AW.ArtifactFile[AW.ArtifactFile(
            url = "https://caltech.box.com/shared/static/8gd3wjq2dbrnzuv8pd5ww3a54h0kkcz8.nc",
            filename = "sst.nc",
        ),],
    )
    return AW.get_data_folder(sst_dataset)
end

function sic_dataset_path()
    sic_dataset = AW.ArtifactWrapper(
        @__DIR__,
        isempty(get(ENV, "CI", "")),
        "sic",
        AW.ArtifactFile[AW.ArtifactFile(
            url = "https://caltech.box.com/shared/static/n0omgqkmnwpr9gylhixew8ywb4psgvj4.nc",
            filename = "sic.nc",
        ),],
    )
    return AW.get_data_folder(sic_dataset)
end

function mask_dataset_path()
    mask_dataset = AW.ArtifactWrapper(
        @__DIR__,
        isempty(get(ENV, "CI", "")),
        "mask",
        AW.ArtifactFile[AW.ArtifactFile(
            url = "https://caltech.box.com/shared/static/vubmq84nhvbgdqayezguf3i1w6nqtwvu.ncc",
            filename = "seamask.nc",
        ),],
    )
    return AW.get_data_folder(mask_dataset)
end

sst_data = joinpath(sst_dataset_path(), "sst.nc")
sic_data = joinpath(sic_dataset_path(), "sic.nc")
mask_data = joinpath(mask_dataset_path(), "seamask.nc")
