import Printf
import ClimaCore as CC
import Makie
import ClimaCoreMakie
import CairoMakie
import ClimaCoupler: Interfacer, ConservationChecker
import ClimaAtmos as CA
import StaticArrays

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

    days = collect(1:length(ccs[1])) * float(coupler_sim.Δt_cpl) / 86400

    # evolution of energy of each component relative to initial value
    total = ccs.total  # total

    var_name = Interfacer.name(cc)
    cum_total = [0.0]
    f = Makie.Figure()
    ax = Makie.Axis(f[1, 1], xlabel = "time [days]", ylabel = "$var_name: (t) - (t=0)")
    Makie.lines!(ax, days, total .- total[1], label = "total"; linewidth = 3)
    for sim in model_sims
        sim_name = Interfacer.name(sim)
        global_field = getproperty(ccs, Symbol(sim_name))
        diff_global_field = (global_field .- global_field[1])
        Makie.lines!(ax, days, diff_global_field[1:length(days)], label = sim_name)
        cum_total .+= abs.(global_field[end])
    end
    if cc isa ConservationChecker.EnergyConservationCheck
        global_field = ccs.toa_net_source
        diff_global_field = (global_field .- global_field[1])
        Makie.lines!(ax, days, diff_global_field[1:length(days)], label = "toa_net")
        cum_total .+= abs.(global_field[end])
    end
    Makie.axislegend(ax, position = :lb)
    Makie.save(figname1, f)

    # use the cumulative global sum at the final time step as a reference for the error calculation
    rse = abs.((total .- total[1]) ./ cum_total)
    l_rse = log.(rse)
    # evolution of log error of total
    lp = Makie.lines(days, l_rse, label = "rs error")
    lp.axis.xlabel = "time [days]"
    lp.axis.ylabel = "log( |x(t) - x(t=0)| / Σx(t=T) )"
    filter_l_rse = filter(x -> !isinf(x) && !isnan(x), l_rse)
    !isempty(filter_l_rse) && Makie.ylims!(minimum(filter_l_rse), maximum(filter_l_rse))
    Makie.axislegend(position = :lt)
    Makie.save(figname2, lp)

    # check that the relative error is small (TODO: reduce this to sqrt(eps(FT)))
    if !softfail
        @info typeof(cc)
        @info rse[end]
        @assert rse[end] < 0.0055
    end
end


# plotting functions for the coupled simulation
"""
    debug(cs::Interfacer.CoupledSimulation, dir = "debug", cs_fields_ref = nothing)

Plot the fields of a coupled simulation and save plots to a directory.
"""
function debug(cs::Interfacer.CoupledSimulation, dir = "debug", cs_fields_ref = nothing)
    isdir(dir) || mkpath(dir)
    @info "plotting debug in $dir"
    for sim in cs.model_sims
        debug(sim, dir)
    end
    debug(cs.fields, dir, cs_fields_ref)
end

"""
    debug(cs_fields::CC.Fields.Field, dir, cs_fields_ref = nothing)

Plot useful coupler fields (in `field_names`) and save plots to a directory.

If `cs_fields_ref` is provided (e.g., using a copy of cs.fields from the initialization),
plot the anomalies of the fields with respect to `cs_fields_ref`.
"""
function debug(cs_fields::CC.Fields.Field, dir, cs_fields_ref = nothing)
    field_names = propertynames(cs_fields)
    fig = Makie.Figure(size = (1500, 800))
    min_square_len = ceil(Int, sqrt(length(field_names)))
    for i in 1:min_square_len, j in 1:min_square_len
        field_index = (i - 1) * min_square_len + j
        if field_index <= length(field_names)
            field_name = field_names[field_index]
            field = getproperty(cs_fields, field_name)

            # Copy field onto cpu space if necessary
            cpu_field = CC.to_cpu(field)

            # ClimaCoreMakie doesn't support NaNs/Infs, so we substitute them with 100max
            FT = CC.Spaces.undertype(axes(cpu_field))
            isinvalid = (x) -> isnan(x) || isinf(x)

            field_valid_min, field_valid_max = extrema(map(x -> isinvalid(x) ? FT(0) : x, parent(cpu_field)))
            map!(x -> isinvalid(x) ? 100field_valid_max : x, parent(cpu_field), parent(cpu_field))
            colorrange =
                isapprox(field_valid_min, field_valid_max) ? (field_valid_min - 1, field_valid_max + 1) :
                (field_valid_min, field_valid_max)
            ax = Makie.Axis(fig[i, j * 2 - 1], title = string(field_name) * print_extrema(field))
            hm = ClimaCoreMakie.fieldheatmap!(ax, cpu_field, colorrange = colorrange)
            Makie.Colorbar(fig[i, j * 2], hm)
        end
    end
    Makie.save(joinpath(dir, "debug_coupler.png"), fig)

    # plot anomalies if a reference cs.fields, `cs_fields_ref`, are provided
    if !isnothing(cs_fields_ref)
        for i in 1:min_square_len, j in 1:min_square_len
            field_index = (i - 1) * min_square_len + j
            if field_index <= length(field_names)
                field_name = field_names[field_index]
                field = getproperty(cs_fields, field_name)
                # Copy field onto cpu space if necessary

                cpu_field = CC.to_cpu(field)

                # ClimaCoreMakie doesn't support NaNs/Infs, so we substitute them with 100max
                FT = CC.Spaces.undertype(axes(cpu_field))
                isinvalid = (x) -> isnan(x) || isinf(x)

                field_valid_min, field_valid_max = extrema(map(x -> isinvalid(x) ? FT(0) : x, parent(cpu_field)))
                map!(x -> isinvalid(x) ? 100field_valid_max : x, parent(cpu_field), parent(cpu_field))
                colorrange =
                    isapprox(field_valid_min, field_valid_max) ? (field_valid_min - 1, field_valid_max + 1) :
                    (field_valid_min, field_valid_max)

                ax = Makie.Axis(fig[i, j * 2 - 1], title = string(field_name) * print_extrema(field))
                hm = ClimaCoreMakie.fieldheatmap!(
                    ax,
                    cpu_field .- getproperty(cs_fields_ref, field_name),
                    colorrange = colorrange,
                )
                Makie.Colorbar(fig[i, j * 2], hm)
            end
        end
        Makie.save(joinpath(dir, "debug_coupler_anomalies.png"), fig)
    end
end

"""
    debug(sim::Interfacer.ComponentModelSimulation, dir)

Plot the fields of a component model simulation and save plots to a directory.
"""
function debug(sim::Interfacer.ComponentModelSimulation, dir)
    field_names = plot_field_names(sim)
    fig = Makie.Figure(size = (1500, 800))
    min_square_len = ceil(Int, sqrt(length(field_names)))
    for i in 1:min_square_len, j in 1:min_square_len
        field_index = (i - 1) * min_square_len + j
        if field_index <= length(field_names)
            field_name = field_names[field_index]
            field = Interfacer.get_field(sim, Val(field_name))
            ax = Makie.Axis(fig[i, j * 2 - 1], title = string(field_name) * print_extrema(field))
            if field isa AbstractArray
                lin = Makie.lines!(ax, Array(field))
            else
                # Copy field onto cpu space if necessary
                cpu_field = CC.to_cpu(field)
                if cpu_field isa CC.Fields.ExtrudedCubedSphereSpectralElementField3D
                    cpu_field = CC.Fields.level(cpu_field, 1)
                end
                FT = CC.Spaces.undertype(axes(cpu_field))
                vals = CC.Fields.field_values(cpu_field)
                isinvalid = (x) -> isnan(x) || isinf(x)
                # ClimaCoreMakie doesn't support NaNs or Infs, so we make them
                # stand out by replacing them with 100max

                field_valid_min, field_valid_max = extrema(map(x -> isinvalid(x) ? FT(0) : x, parent(cpu_field)))
                map!(x -> isinvalid(x) ? 100field_max : x, parent(cpu_field), parent(cpu_field))
                colorrange =
                    isapprox(field_valid_min, field_valid_max) ? (field_valid_min - 1, field_valid_max + 1) :
                    (field_valid_min, field_valid_max)
                hm = ClimaCoreMakie.fieldheatmap!(ax, cpu_field, colorrange = colorrange)
                Makie.Colorbar(fig[i, j * 2], hm)
            end
        end
    end
    Makie.save(joinpath(dir, "debug_$(Interfacer.name(sim)).png"), fig)
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

# additional ClimaLand model debug fields
Interfacer.get_field(sim::ClimaLandSimulation, ::Val{:soil_water}) = sim.integrator.u.soil.ϑ_l
Interfacer.get_field(sim::ClimaLandSimulation, ::Val{:soil_ice}) = sim.integrator.u.soil.θ_i
Interfacer.get_field(sim::ClimaLandSimulation, ::Val{:soil_energy}) = sim.integrator.u.soil.ρe_int
Interfacer.get_field(sim::ClimaLandSimulation, ::Val{:canopy_temp}) = sim.integrator.u.canopy.energy.T
Interfacer.get_field(sim::ClimaLandSimulation, ::Val{:canopy_water}) = sim.integrator.u.canopy.hydraulics.ϑ_l.:1
Interfacer.get_field(sim::ClimaLandSimulation, ::Val{:snow_energy}) = sim.integrator.u.snow.U
Interfacer.get_field(sim::ClimaLandSimulation, ::Val{:snow_water_equiv}) = sim.integrator.u.snow.S
Interfacer.get_field(sim::ClimaLandSimulation, ::Val{:snow_liquid_water}) = sim.integrator.u.snow.S_l

# currently selected plot fields
plot_field_names(sim::Interfacer.SurfaceModelSimulation) = (:area_fraction, :surface_temperature)
plot_field_names(sim::ClimaLandSimulation) = (
    :area_fraction,
    :surface_direct_albedo,
    :surface_diffuse_albedo,
    :surface_temperature,
    :soil_water,
    :soil_ice,
    :soil_energy,
    :canopy_temp,
    :canopy_water,
    :snow_energy,
    :snow_water_equiv,
    :snow_liquid_water,
)
plot_field_names(sim::BucketSimulation) = (:area_fraction, :surface_temperature, :σS, :Ws, :W)
plot_field_names(sim::ClimaAtmosSimulation) = (:w, :ρq_tot, :ρe_tot, :liquid_precipitation, :snow_precipitation)
