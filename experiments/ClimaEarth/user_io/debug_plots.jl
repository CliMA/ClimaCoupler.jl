import Plots
import Printf
import ClimaComms
@static pkgversion(ClimaComms) >= v"0.6" && ClimaComms.@import_required_backends
import ClimaCore as CC
import ClimaCorePlots
import ClimaCoupler: Interfacer, ConservationChecker
import ClimaAtmos as CA
import StaticArrays

include("plot_helper.jl")

"""
    plot_global_conservation(
        cc::AbstractConservationCheck,
        coupler_sim::Interfacer.CoupledSimulation,
        softfail = false;
        figname1 = "total.png",
        figname2 = "total_log.png",
    )

Creates two plots of the globally integrated quantity (energy, ``\\rho e``):
1. global quantity of each model component as a function of time,
relative to the initial value;
2. fractional change in the sum of all components over time on a log scale.
"""
function plot_global_conservation(
    cc::ConservationChecker.AbstractConservationCheck,
    coupler_sim::Interfacer.CoupledSimulation,
    softfail = false;
    figname1 = "total.png",
    figname2 = "total_log.png",
)

    model_sims = coupler_sim.model_sims
    ccs = cc.sums

    days = collect(1:length(ccs[1])) * coupler_sim.Δt_cpl / 86400

    # evolution of energy of each component relative to initial value
    total = ccs.total  # total

    var_name = Interfacer.name(cc)
    cum_total = [0.0]

    Plots.plot(
        days,
        total .- total[1],
        label = "total",
        xlabel = "time [days]",
        ylabel = "$var_name: (t) - (t=0)",
        linewidth = 3,
    )
    for sim in model_sims
        sim_name = Interfacer.name(sim)
        global_field = getproperty(ccs, Symbol(sim_name))
        diff_global_field = (global_field .- global_field[1])
        Plots.plot!(days, diff_global_field[1:length(days)], label = sim_name)
        cum_total .+= abs.(global_field[end])
    end
    if cc isa ConservationChecker.EnergyConservationCheck
        global_field = ccs.toa_net_source
        diff_global_field = (global_field .- global_field[1])
        Plots.plot!(days, diff_global_field[1:length(days)], label = "toa_net")
        cum_total .+= abs.(global_field[end])
    end
    Plots.savefig(figname1)

    # use the cumulative global sum at the final time step as a reference for the error calculation
    rse = abs.((total .- total[1]) ./ cum_total)

    # evolution of log error of total
    Plots.plot(days, log.(rse), label = "rs error", xlabel = "time [days]", ylabel = "log( |x(t) - x(t=0)| / Σx(t=T) )")
    Plots.savefig(figname2)

    # check that the relative error is small (TODO: reduce this to sqrt(eps(FT)))
    if !softfail
        @info typeof(cc)
        @info rse[end]
        @assert rse[end] < 5e-3
    end
end


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
