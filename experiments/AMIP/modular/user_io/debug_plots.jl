using Plots
using ClimaCorePlots
using Printf
using ClimaCoupler.Interfacer: ComponentModelSimulation, SurfaceModelSimulation

# plotting functions for the coupled simulation
"""
    debug(cs::CoupledSimulation, dir = "debug")

Plot the fields of a coupled simulation and save plots to a directory.
"""
function debug(cs::CoupledSimulation, dir = "debug")
    mkpath(dir)
    @info "plotting debug in " * dir
    for sim in cs.model_sims
        debug(sim, dir)
    end
    debug(cs.fields, dir)
end

"""
    debug(cs_fields::NamedTuple, dir)

Plot useful coupler fields (in `field_names`) and save plots to a directory.
"""
function debug(cs_fields::NamedTuple, dir)
    field_names = (:F_turb_energy, :F_turb_moisture, :P_liq, :T_S, :ρ_sfc, :q_sfc)
    all_plots = []
    for field_name in field_names
        field = getproperty(cs_fields, field_name)
        push!(all_plots, Plots.plot(field, title = string(field_name) * print_extrema(field)))
    end
    fig = Plots.plot(all_plots..., size = (1500, 800))
    Plots.png(joinpath(dir, "debug_coupler"))
end

"""
    debug(sim::ComponentModelSimulation, dir)

Plot the fields of a component model simulation and save plots to a directory.
"""
function debug(sim::ComponentModelSimulation, dir)

    field_names = plot_field_names(sim)

    all_plots = []
    for field_name in field_names
        field = get_field(sim, Val(field_name))
        push!(all_plots, Plots.plot(field, title = string(field_name) * print_extrema(field)))
    end
    fig = Plots.plot(all_plots..., size = (1500, 800))
    Plots.png(joinpath(dir, "debug_$(name(sim))"))

end

"""
    print_extrema(field::ClimaCore.Fields.Field)

Return the minimum and maximum values of a field as a string.
"""
function print_extrema(field::ClimaCore.Fields.Field)
    ext_vals = extrema(field)
    min = @sprintf("%.2E", ext_vals[1])
    max = @sprintf("%.2E", ext_vals[2])
    return " [$min, $max]"
end

# below are additional fields specific to this experiment (ourside of the required coupler fields) that we are interested in plotting for debugging purposes

# additional ClimaAtmos model debug fields
function get_field(sim::ClimaAtmosSimulation, ::Val{:w})
    w_c = ones(sim.domain.face_space.horizontal_space)
    parent(w_c) .= parent(Fields.level(Geometry.WVector.(sim.integrator.u.f.u₃), 5 .+ half))
    return w_c
end
get_field(sim::ClimaAtmosSimulation, ::Val{:ρq_tot}) = sim.integrator.u.c.ρq_tot
get_field(sim::ClimaAtmosSimulation, ::Val{:ρe_tot}) = sim.integrator.u.c.ρe_tot

# additional BucketSimulation debug fields
get_field(sim::BucketSimulation, ::Val{:σS}) = sim.integrator.u.bucket.σS
get_field(sim::BucketSimulation, ::Val{:Ws}) = sim.integrator.u.bucket.Ws
get_field(sim::BucketSimulation, ::Val{:W}) = sim.integrator.u.bucket.W

# currently selected plot fields
plot_field_names(sim::SurfaceModelSimulation) = (:area_fraction, :surface_temperature, :surface_humidity)
plot_field_names(sim::BucketSimulation) =
    (:area_fraction, :surface_temperature, :surface_humidity, :air_density, :σS, :Ws, :W)
plot_field_names(sim::ClimaAtmosSimulation) = (:w, :ρq_tot, :ρe_tot, :liquid_precipitation, :snow_precipitation)


#=
# offline debuging
using ClimaCore
using ClimaCore: InputOutput, Fields, Geometry
using Plots
using ClimaCorePlots
using ClimaCore.Utilities: half, PlusHalf
import ClimaAtmos.Parameters as CAP


include("../components/atmosphere/climaatmos_init.jl")

# monthly mean data
output_dir = "/central/scratch/esm/slurm-buildkite/climacoupler-longruns/239/climacoupler-longruns/experiments/AMIP/modular/output/amip/new_target" # hdf5
output_dir_atmos = "/central/scratch/esm/slurm-buildkite/climacoupler-longruns/239/climacoupler-longruns/new_target/" # hdf5

# get the monthly mean data from coupler diagnostics
vname = "u"
filename = joinpath(output_dir,"$vname.monthly_mean_3d_.1979-01-01T00:02:30.hdf5")
reader = InputOutput.HDF5Reader(filename)
monthly_mean = InputOutput.read_field(reader, vname)
close(reader)

# get atmos data from coupler restarts


# get atmos data from atmos restarts
filename_atmos = joinpath(output_dir_atmos,"day10.0.hdf5")
reader_atmos = InputOutput.HDF5Reader(filename_atmos)
Y = InputOutput.read_field(reader_atmos, "Y")
close(reader_atmos)

face_space = axes(Y.f)

all_plots = []
push!(all_plots, Plots.plot(Y.c.ρq_tot, title = "ρq_tot", size = (1500, 800)))
push!(all_plots, Plots.plot(Y.c.ρe_tot, title = "ρe_tot", size = (1500, 800)))
push!(all_plots, Plots.plot(ClimaCore.Geometry.UVVector.(Y.c.uₕ).components.data.:1, title = "uₕ", size = (1500, 800)))
push!(all_plots, Plots.plot(ClimaCore.Geometry.WVector.(Y.f.u₃).components.data.:1, title = "w", size = (1500, 800), level = half))
fig = Plots.plot(all_plots..., size = (1500, 800))
Plots.png("atmos_states")

# coordinates of extrema
using ClimaCore: Spaces
function get_max_coords(field)
    max_val, max_idx = findmax(parent(field))
    @info "Printing coordinates of max: $max_val"
    lat = Fields.coordinate_field(field).lat
    println("lat: "*"$(parent(lat)[max_idx])")
    lon = Fields.coordinate_field(field).long
    println("lon: "*"$(parent(lon)[max_idx])")
    if axes(field) isa Spaces.ExtrudedFiniteDifferenceSpace
        z =  Fields.coordinate_field(field).z
        println("lz: " * "$(parent(z)[max_idx])")
    end
end
=#