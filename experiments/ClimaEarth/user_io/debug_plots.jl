import Plots
import Printf
import ClimaComms
import ClimaCore as CC
import ClimaCorePlots
import ClimaCoupler: Interfacer
import ClimaAtmos as CA
import StaticArrays

include("plot_helper.jl")

# plotting functions for the coupled simulation
"""
    debug(cs::Interfacer.CoupledSimulation, dir = "debug", cs_fields_ref = nothing)

Plot the fields of a coupled simulation and save plots to a directory.
"""
function debug(cs::Interfacer.CoupledSimulation, dir = "debug", cs_fields_ref = nothing)
    isdir(dir) || mkpath(dir)
    @info "plotting debug in " * dir
    for sim in cs.model_sims
        debug(sim, dir)
    end
    debug(cs.fields, dir, cs_fields_ref)
end

"""
    debug(cs_fields::NamedTuple, dir, cs_fields_ref = nothing)

Plot useful coupler fields (in `field_names`) and save plots to a directory.

If `cs_fields_ref` is provided (e.g., using a copy of cs.fields from the initialization),
plot the anomalies of the fields with respect to `cs_fields_ref`.
"""
function debug(cs_fields::NamedTuple, dir, cs_fields_ref = nothing)
    field_names = (
        :surface_direct_albedo,
        :surface_diffuse_albedo,
        :F_radiative,
        :F_turb_energy,
        :F_turb_moisture,
        :F_turb_ρτxz,
        :F_turb_ρτyz,
        :P_liq,
        :P_snow,
        :T_S,
        :ρ_sfc,
        :q_sfc,
        :beta,
        :z0b_S,
        :z0m_S,
        :radiative_energy_flux_toa,
    )
    all_plots = []

    for field_name in field_names
        field = getproperty(cs_fields, field_name)

        # Copy field onto cpu space if necessary
        cpu_field = to_cpu(field)

        push!(all_plots, Plots.plot(cpu_field, title = string(field_name) * print_extrema(field)))
    end
    Plots.plot(all_plots..., size = (1500, 800))
    Plots.png(joinpath(dir, "debug_coupler"))

    # plot anomalies if a reference cs.fields, `cs_fields_ref`, are provided
    if !isnothing(cs_fields_ref)
        all_plots = []
        for field_name in field_names
            field = getproperty(cs_fields, field_name)
            # Copy field onto cpu space if necessary
            cpu_field = to_cpu(field)

            push!(
                all_plots,
                Plots.plot(
                    cpu_field .- getproperty(cs_fields_ref, field_name),
                    title = string(field_name) * print_extrema(field),
                ),
            )
        end
        Plots.plot(all_plots..., size = (1500, 800))
        Plots.png(joinpath(dir, "debug_coupler_amomalies"))
    end
end

"""
    debug(sim::Interfacer.ComponentModelSimulation, dir)

Plot the fields of a component model simulation and save plots to a directory.
"""
function debug(sim::Interfacer.ComponentModelSimulation, dir)
    field_names = plot_field_names(sim)

    all_plots = []
    for field_name in field_names
        field = Interfacer.get_field(sim, Val(field_name))

        # Copy field onto cpu space if necessary
        cpu_field = to_cpu(field)

        push!(all_plots, Plots.plot(cpu_field, title = string(field_name) * print_extrema(field)))
    end
    fig = Plots.plot(all_plots..., size = (1500, 800))
    Plots.png(joinpath(dir, "debug_$(Interfacer.name(sim))"))
end

"""
    print_extrema(field::CC.Fields.Field)

Return the minimum and maximum values of a field as a string.
"""
function print_extrema(field::Union{CC.Fields.Field, Vector, StaticArrays.SVector})
    ext_vals = extrema(field)
    min = Printf.@sprintf("%.2E", ext_vals[1])
    max = Printf.@sprintf("%.2E", ext_vals[2])
    return " [$min, $max]"
end

# below are additional fields specific to this experiment (ourside of the required coupler fields) that we are interested in plotting for debugging purposes

# additional ClimaAtmos model debug fields
function Interfacer.get_field(sim::ClimaAtmosSimulation, ::Val{:w})
    w_c = ones(CC.Spaces.horizontal_space(sim.domain.face_space))
    parent(w_c) .= parent(CC.Fields.level(CC.Geometry.WVector.(sim.integrator.u.f.u₃), 5 .+ CC.Utilities.half))
    return w_c
end
specific_humidity(::CA.DryModel, integrator) = [eltype(integrator.u)(0)]
specific_humidity(::Union{CA.EquilMoistModel, CA.NonEquilMoistModel}, integrator) = integrator.u.c.ρq_tot
Interfacer.get_field(sim::ClimaAtmosSimulation, ::Val{:ρq_tot}) =
    specific_humidity(sim.integrator.p.atmos.moisture_model, sim.integrator)
Interfacer.get_field(sim::ClimaAtmosSimulation, ::Val{:ρe_tot}) = sim.integrator.u.c.ρe_tot

# additional BucketSimulation debug fields
Interfacer.get_field(sim::BucketSimulation, ::Val{:σS}) = sim.integrator.u.bucket.σS
Interfacer.get_field(sim::BucketSimulation, ::Val{:Ws}) = sim.integrator.u.bucket.Ws
Interfacer.get_field(sim::BucketSimulation, ::Val{:W}) = sim.integrator.u.bucket.W

# currently selected plot fields
plot_field_names(sim::Interfacer.SurfaceModelSimulation) = (:area_fraction, :surface_temperature, :surface_humidity)
plot_field_names(sim::BucketSimulation) =
    (:area_fraction, :surface_temperature, :surface_humidity, :air_density, :σS, :Ws, :W)
plot_field_names(sim::ClimaAtmosSimulation) = (:w, :ρq_tot, :ρe_tot, :liquid_precipitation, :snow_precipitation)
