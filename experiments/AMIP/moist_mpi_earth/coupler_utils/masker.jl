function LandSeaMask(FT, infile, varname, boundary_space; outfile = "land_sea_cgll.nc", threshold = 0.7)
    weightfile, datafile_cgll = ncreader_rll_to_cgll_from_space(infile, varname, boundary_space, outfile = outfile)
    mask = ncreader_cgll_sparse_to_field(datafile_cgll, varname, weightfile, (Int(1),), boundary_space)[1]
    mask = swap_space!(mask, boundary_space) # needed if we are reading from previous run
    return mask
end

combine_surface(FT, mask, sfc_1, sfc_2, sfc_3, value1 = -0.5, value2 = 0.5) =
    (mask < FT(value1) ? sfc_3 : FT(0)) +
    ((mask >= FT(value1) && (mask <= FT(value2))) ? sfc_2 : FT(0)) +
    (mask > FT(value2) ? sfc_1 : FT(0))

"""
apply_mask(T, mask, condition, field; value = 0.5) 

"""
apply_mask(T, mask, condition, field; value = 0.5) = condition(mask, value) ? field : T(0.0)

"""
    time_slice_ncfile(sic_data, time_idx = 1)
- slices a dataset at time index `time_idx` and saves it under `sic_data_slice`. Used for more efficient regridding of mask, SST and SIC files. 
"""
function time_slice_ncfile(sic_data, time_idx = 1)
    sic_data_slice = sic_data[1:(end - 3)] * "_one_time.nc"
    isfile(sic_data_slice) ? run(`rm $sic_data_slice`) : nothing
    NCDataset(sic_data) do ds
        write(sic_data_slice, ds, idimensions = Dict("time" => time_idx:time_idx))
    end
    return sic_data_slice
end
