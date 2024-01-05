using ClimaComms
using ClimaCoupler: TestHelper, Regridder, BCReader
using ClimaCoupler
using Dates

# get artifacts
include(joinpath(pkgdir(ClimaCoupler), "artifacts", "artifact_funcs.jl"))
land_mask_data = joinpath(mask_dataset_path(), "seamask.nc")
sst_data = joinpath(sst_dataset_path(), "sst.nc")

# set up comms context for MPI
comms_ctx = ClimaComms.MPICommsContext()
pid, nprocs = ClimaComms.init(comms_ctx)

# more setup
FT = Float64
REGRID_DIR = joinpath(pkgdir(ClimaCoupler), "debug", "regrid")
boundary_space = TestHelper.create_space(FT; comms_ctx)
mono_surface = false
date0 = DateTime("19790101", dateformat"yyyymmdd")
varname = "SST"

# get land fraction
land_fraction =
    FT.(
        Regridder.land_fraction(
            FT,
            REGRID_DIR,
            comms_ctx,
            land_mask_data,
            "LSMASK",
            boundary_space,
            mono = mono_surface,
        )
    )

# regrid all times, save data to HDF5 files, and get BCFileInfo object
bcf_info = BCReader.bcfile_info_init(
    FT,
    REGRID_DIR,
    sst_data,
    varname,
    boundary_space,
    comms_ctx,
    interpolate_daily = true,
    scaling_function = nothing, ## convert to fraction
    land_fraction = land_fraction,
    date0 = date0,
    mono = mono_surface,
)

# extract info from BCFileInfo object
(; bcfile_dir, comms_ctx, hd_outfile_root, varname, all_dates, scaling_function) = bcf_info
midmonth_idx0 = bcf_info.segment_idx0[1]

Regridder.read_from_hdf5(bcfile_dir, hd_outfile_root, all_dates[Int(midmonth_idx0)], varname, comms_ctx)
