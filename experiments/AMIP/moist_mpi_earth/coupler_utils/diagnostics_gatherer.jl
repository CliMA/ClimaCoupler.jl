# output_dumper.jl
# dumps specified output on model grid

"""
    save_output_func(diags::NamedTuple, date::DateTime, output_dir::String)

Saves all entries in `diags` in separate HDF5 files in `output_dir`. 
It can be callled by @calendar_callback to dump snapshots at a particular time frequency
"""
function save_hdf5(comms_ctx, diags::NamedTuple, date::DateTime, output_dir::String; name_tag = "")
    diags_names = propertynames(diags)
    diags_values = map(x -> getproperty(diags, x), diags_names)

    @info "Coupler: dumping data on $date"

    # convert to 1st of month (convention for monthly data)
    date_m1 = date - Dates.Day(1)
    year_month = Dates.DateTime(Dates.yearmonth(date_m1)[1], Dates.yearmonth(date_m1)[2])

    for (name, values) in zip(diags_names, diags_values)
        output_file = joinpath(output_dir, name_tag * "$name.monthly_" * string(year_month) * ".hdf5")
        hdfwriter = InputOutput.HDF5Writer(output_file, comms_ctx)
        InputOutput.HDF5.write_attribute(hdfwriter.file, "unix time", Dates.datetime2unix(year_month)) # TODO: a better way to write metadata, CMIP convention
        InputOutput.write!(hdfwriter, values, string(name))
        Base.close(hdfwriter)
    end
    return nothing
end

"""
    accumulate_diags(diags::Union{NamedTuple, Fields.Field}, diags_cache::NamedTuple)

Accumulates all entries in `diags` and saves the result in `diags_cache.fields`, while increasing the `diags_cache.ct` counter. 
"""
function accumulate_diags(diags::Union{NamedTuple, Fields.Field}, diags_cache::NamedTuple)
    diags_names = propertynames(diags_cache.fields)

    for (v, name) in zip(diags_cache.fields, diags_names)
        v .+= getproperty(diags, name) # @. breaks for NamedTuples in test
    end
    diags_cache.ct[1] += FT(1)

    return nothing
end

"""
    collect_diags(cs, diags_names)

collect diagnostics in diags names
"""
function collect_diags(cs, diags_names)

    diags = (;)

    for name in diags_names
        diags = (; diags..., zip((name,), (get_var(cs, Val(name)),))...)
    end

    return diags
end
