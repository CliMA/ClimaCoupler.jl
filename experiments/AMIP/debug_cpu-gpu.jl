using ClimaComms
using ClimaCore: InputOutput, Fields
using Statistics

# atmos rho e
varname = "atmos_ρe_tot"
comms_ctx = ClimaComms.SingletonCommsContext()

cpu_filename = "experiments/AMIP/cpu-gpu_temporalmap/cpu_temporalmap_atmos_ρe_tot_1979-01-01T00-03-20.hdf5"
cpu_hdf5reader = InputOutput.HDF5Reader(cpu_filename, comms_ctx)
cpu_atmos_ρe_tot = InputOutput.read_field(cpu_hdf5reader, varname)
close(cpu_hdf5reader)

gpu_filename = "experiments/AMIP/cpu-gpu_temporalmap/gpu_temporalmap_atmos_ρe_tot_1979-01-01T00-03-20.hdf5"
gpu_hdf5reader = InputOutput.HDF5Reader(gpu_filename, comms_ctx)
gpu_atmos_ρe_tot = InputOutput.read_field(gpu_hdf5reader, varname)
close(gpu_hdf5reader)

cpu_atmos_ρe_tot == gpu_atmos_ρe_tot # false
abs(mean(parent(cpu_atmos_ρe_tot)) - mean(parent(gpu_atmos_ρe_tot))) < eps(Float64) # false temporal map (3.64e-12), true function

# atmos rho q
varname = "atmos_ρq_tot"
comms_ctx = ClimaComms.SingletonCommsContext()

cpu_filename = "experiments/AMIP/cpu_function_atmos_ρq_tot_1979-01-01T00-03-20.hdf5"
cpu_hdf5reader = InputOutput.HDF5Reader(cpu_filename, comms_ctx)
cpu_atmos_ρq_tot = InputOutput.read_field(cpu_hdf5reader, varname)
close(cpu_hdf5reader)

gpu_filename = "experiments/AMIP/gpu_function_atmos_ρq_tot_1979-01-01T00-03-20.hdf5"
gpu_hdf5reader = InputOutput.HDF5Reader(gpu_filename, comms_ctx)
gpu_atmos_ρq_tot = InputOutput.read_field(gpu_hdf5reader, varname)
close(gpu_hdf5reader)

cpu_atmos_ρq_tot == gpu_atmos_ρq_tot # false temporal map, true function
abs(mean(parent(cpu_atmos_ρq_tot)) - mean(parent(gpu_atmos_ρq_tot))) < eps(Float64) # true

# land temp
varname = "land_T"
comms_ctx = ClimaComms.SingletonCommsContext()

cpu_filename = "experiments/AMIP/cpu_function_land_T_1979-01-01T00-03-20.hdf5"
cpu_hdf5reader = InputOutput.HDF5Reader(cpu_filename, comms_ctx)
cpu_land_T = InputOutput.read_field(cpu_hdf5reader, varname)
close(cpu_hdf5reader)

gpu_filename = "experiments/AMIP/gpu_function_land_T_1979-01-01T00-03-20.hdf5"
gpu_hdf5reader = InputOutput.HDF5Reader(gpu_filename, comms_ctx)
gpu_land_T = InputOutput.read_field(gpu_hdf5reader, varname)
close(gpu_hdf5reader)

cpu_land_T == gpu_land_T # false
abs(mean(parent(cpu_land_T)) - mean(parent(gpu_land_T))) < eps(Float64) # true

# land water
varname = "land_W"
comms_ctx = ClimaComms.SingletonCommsContext()

cpu_filename = "experiments/AMIP/cpu_function_land_W_1979-01-01T00-03-20.hdf5"
cpu_hdf5reader = InputOutput.HDF5Reader(cpu_filename, comms_ctx)
cpu_land_W = InputOutput.read_field(cpu_hdf5reader, varname)
close(cpu_hdf5reader)

gpu_filename = "experiments/AMIP/gpu_function_land_W_1979-01-01T00-03-20.hdf5"
gpu_hdf5reader = InputOutput.HDF5Reader(gpu_filename, comms_ctx)
gpu_land_W = InputOutput.read_field(gpu_hdf5reader, varname)
close(gpu_hdf5reader)

cpu_land_W == gpu_land_W # false
abs(mean(parent(cpu_land_W)) - mean(parent(gpu_land_W))) < eps(Float64) # true

# ocean temp
varname = "ocean_T_sfc"
comms_ctx = ClimaComms.SingletonCommsContext()

cpu_filename = "experiments/AMIP/cpu_function_ocean_T_sfc_1979-01-01T00-03-20.hdf5"
cpu_hdf5reader = InputOutput.HDF5Reader(cpu_filename, comms_ctx)
cpu_ocean_T_sfc = InputOutput.read_field(cpu_hdf5reader, varname)
close(cpu_hdf5reader)

gpu_filename = "experiments/AMIP/gpu_function_ocean_T_sfc_1979-01-01T00-03-20.hdf5"
gpu_hdf5reader = InputOutput.HDF5Reader(gpu_filename, comms_ctx)
gpu_ocean_T_sfc = InputOutput.read_field(gpu_hdf5reader, varname)
close(gpu_hdf5reader)

cpu_ocean_T_sfc == gpu_ocean_T_sfc # false
abs(mean(parent(cpu_ocean_T_sfc)) - mean(parent(gpu_ocean_T_sfc))) < eps(Float64) # true temporal map, false function (1.13e-13)
