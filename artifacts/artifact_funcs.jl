import ArtifactWrappers as AW

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

"""
    artifact_data(datapath_full, filename)

Returns input dataset at datapath_full
"""
function artifact_data(datapath_full, filename)
    datafile_truncated = joinpath(datapath_full, string(lowercase(filename), ".nc"))
    return datafile_truncated
end

"""
    artifact_data(datapath_full, filename, varname, datapath_trunc, date0, t_start, t_end, comms_ctx)

Truncates given data set, and constructs a new dataset containing only
the dates needed and stores it in datapath_trunc
"""
function artifact_data(datapath_full, filename, varname, datapath_trunc, date0, t_start, t_end, comms_ctx)
    datafile = joinpath(datapath_full, string(lowercase(filename), ".nc"))
    datafile_truncated =
        Regridder.truncate_dataset(datafile, filename, varname, datapath_trunc, date0, t_start, t_end, comms_ctx)
    return datafile_truncated
end
