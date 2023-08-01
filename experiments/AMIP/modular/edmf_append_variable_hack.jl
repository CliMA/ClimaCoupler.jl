
function add_atmos_cache(sim, input_dir = "data/", t = 2678400)

    input_file = joinpath(input_dir, "checkpoint", "checkpoint_" * Interfacer.name(sim) * "_$t.hdf5")

    hdfreader = InputOutput.HDF5Reader(input_file)
    Y = InputOutput.read_field(hdfreader, "model_state")
    t = InputOutput.HDF5.read_attribute(hdfreader.file, "time")
    FT = eltype(Y)
    sgs_new = (;ρatke = FT(0))
    Yc_new = map(nt -> (; nt..., sgs⁰= sgs_new), Y.c)
    Y = Fields.FieldVector(; c=Yc_new, f=Y.f)
    # This has to be a new file that doesn't exist yet
    output_dir = mkpath(joinpath(input_dir, "edmf_cache_hack"))
    output_file = joinpath(output_dir, "checkpoint", "checkpoint_" * Interfacer.name(sim) * "_$t.hdf5")
    if @isdefined(hdfwriter)
        Base.close(hdfwriter)
    end
    hdfwriter = InputOutput.HDF5Writer(output_file, ClimaComms.SingletonCommsContext())
    InputOutput.HDF5.write_attribute(hdfwriter.file, "time", t)
    InputOutput.write!(hdfwriter, Y, "model_state")
    @info "Added sgs⁰ to Atmos restart. New file is $output_file"
    Base.close(hdfwriter)
end

add_atmos_cache(atmos_sim, input_dir, t)


# orinput_dir = "data/"
# t = 2678400

# input_file = joinpath(input_dir, "checkpoint", "checkpoint_" * Interfacer.name(atmos_sim) * "_$t.hdf5")
# hdfreader = InputOutput.HDF5Reader(input_file)
# Y = InputOutput.read_field(hdfreader, "model_state")
# t = InputOutput.HDF5.read_attribute(hdfreader.file, "time")
# FT = eltype(Y)
# Yc_new = map(nt -> (; nt..., sgs⁰=(;ρatke = FT(0))), Y.c)
# Y = Fields.FieldVector(; c=Yc_new, f=Y.f)
# # This has to be a new file that doesn't exist yet
# output_dir = mkpath(joinpath(input_dir, "edmf_cache_hack"))
# output_file = joinpath(output_dir, "checkpoint", "checkpoint_" * Interfacer.name(atmos_sim) * "_$t.hdf5")
# if @isdefined(hdfwriter)
#     Base.close(hdfwriter)
# end
# hdfwriter = InputOutput.HDF5Writer(output_file, ClimaComms.SingletonCommsContext())
# InputOutput.HDF5.write_attribute(hdfwriter.file, "time", t)
# InputOutput.write!(hdfwriter, Y, "model_state")
# @info "Added sgs⁰ to Atmos restart. New file is $output_file"
