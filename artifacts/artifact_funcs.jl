import ArtifactWrappers as AW

function sst_dataset_path()
    sst_dataset = AW.ArtifactWrapper(
        @__DIR__,
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
        "sic",
        AW.ArtifactFile[AW.ArtifactFile(
            url = "https://caltech.box.com/shared/static/n0omgqkmnwpr9gylhixew8ywb4psgvj4.nc",
            filename = "sic.nc",
        ),],
    )
    return AW.get_data_folder(sic_dataset)
end

function co2_dataset_path()
    co2_dataset = AW.ArtifactWrapper(
        @__DIR__,
        "co2",
        AW.ArtifactFile[AW.ArtifactFile(
            url = "https://caltech.box.com/shared/static/xg028wnsn57wam6euwrh98fe43ibei8g.nc",
            filename = "mauna_loa_co2.nc",
        ),],
    )
    return AW.get_data_folder(co2_dataset)
end

function mask_dataset_path()
    mask_dataset = AW.ArtifactWrapper(
        @__DIR__,
        "land_mask",
        AW.ArtifactFile[AW.ArtifactFile(
            url = "https://caltech.box.com/shared/static/vubmq84nhvbgdqayezguf3i1w6nqtwvu.ncc",
            filename = "seamask.nc",
        ),],
    )
    return AW.get_data_folder(mask_dataset)
end
