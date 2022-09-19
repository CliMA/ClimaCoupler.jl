function LandSeaMask(FT, infile, varname, boundary_space; outfile = "land_sea_cgll.nc", mono = false, threshold = 0.7)
    weightfile, datafile_cgll =
        ncreader_rll_to_cgll_from_space(infile, varname, boundary_space, outfile = outfile, mono = mono)
    mask = ncreader_cgll_sparse_to_field(datafile_cgll, varname, weightfile, (Int(1),), boundary_space)[1]
    mask = swap_space!(mask, boundary_space) # needed if we are reading from previous run
    return mono ? mask : binary_mask.(mask, threshold = threshold)
end

"""
    binary_mask(var::FT; threshold = 0.5)

Converts a number to 1 or 0 of the same type, based on a threashold. 
"""

binary_mask(var::FT; threshold = 0.5) where {FT} = (var - FT(threshold)) > FT(0) ? FT(1) : FT(0)

"""
    combine_surfaces!(combined_field::Fields.Field, masks::NamedTuple, fields::NamedTuple)

Sums Field objects in `fields` weighted by the respective masks.
"""
function combine_surfaces!(combined_field::Fields.Field, masks::NamedTuple, fields::NamedTuple)
    combined_field .= eltype(combined_field)(0)
    for surface_name in propertynames(fields) # could use dot here?
        field_no_nans = nans_to_zero.(getproperty(fields, surface_name))  # TODO: performance analysis / alternatives
        combined_field .+= getproperty(masks, surface_name) .* field_no_nans
    end

end
nans_to_zero(v) = isnan(v) ? FT(0) : v

"""
    update_masks(cs)

Updates dynamically changing masks. 
"""
function update_masks(cs)

    # dynamic masks
    ice_d = cs.model_sims.ice_sim.integrator.p.ice_mask
    FT = eltype(ice_d)

    # static mask
    land_s = cs.surface_masks.land

    cs.surface_masks.ice .= min.(ice_d .+ land_s, FT(1)) .- land_s
    cs.surface_masks.ocean .= (FT(1) .- cs.surface_masks.ice .- land_s)

    @assert minimum(cs.surface_masks.ice) >= FT(0)
    @assert minimum(cs.surface_masks.land) >= FT(0)
    @assert minimum(cs.surface_masks.ocean) >= FT(0)

end

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
