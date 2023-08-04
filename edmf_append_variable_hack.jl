
using ClimaCore: InputOutput, Fields
using ClimaComms

function append_atmos_variable(;sim_name = "ClimaAtmosSimulation", input_dir = "data/", output_dir = "data_new/", t = 2678400)

    input_file = joinpath(input_dir, "checkpoint", "checkpoint_" * sim_name * "_$t.hdf5")
    output_file = joinpath(output_dir, "checkpoint", "checkpoint_" * sim_name * "_$t.hdf5")

    hdfreader = InputOutput.HDF5Reader(input_file)
    Y = InputOutput.read_field(hdfreader, "model_state")
    t = InputOutput.HDF5.read_attribute(hdfreader.file, "time")
    FT = eltype(Y)
    sgs_new = (;ρatke = FT(0))
    Yc_new = map(nt -> (; nt..., sgs⁰= sgs_new), Y.c)
    Y = Fields.FieldVector(; c=Yc_new, f=Y.f)
    if @isdefined(hdfwriter)
        Base.close(hdfwriter)
    end
    hdfwriter = InputOutput.HDF5Writer(output_file, ClimaComms.SingletonCommsContext())
    InputOutput.HDF5.write_attribute(hdfwriter.file, "time", t)
    InputOutput.write!(hdfwriter, Y, "model_state")
    @info "Added sgs⁰ to Atmos restart. New file is $output_file"
    Base.close(hdfwriter)

end

append_atmos_variable(input_dir = ENV["RESTART_DIR"], output_dir = ENV["NEW_RESTART_DIR"], t = ENV["RESTART_T"])

