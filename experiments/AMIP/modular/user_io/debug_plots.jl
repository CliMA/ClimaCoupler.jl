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
function get_field(atmos_sim::ClimaAtmosSimulation, ::Val{:F_radiative_sfc})
    radiation = atmos_sim.integrator.p.radiation_model
    FT = eltype(atmos_sim.integrator.u)
    # save radiation source
    if radiation != nothing
        face_space = axes(atmos_sim.integrator.u.f)
        z = parent(Fields.coordinate_field(face_space).z)
        Δz_top = round(FT(0.5) * (z[end, 1, 1, 1, 1] - z[end - 1, 1, 1, 1, 1]))
        n_faces = length(z[:, 1, 1, 1, 1])

        LWd_TOA = Fields.level(
            CA.RRTMGPI.array2field(FT.(atmos_sim.integrator.p.radiation_model.face_lw_flux_dn), face_space),
            half,
        )
        LWu_TOA = Fields.level(
            CA.RRTMGPI.array2field(FT.(atmos_sim.integrator.p.radiation_model.face_lw_flux_up), face_space),
            half,
        )
        SWd_TOA = Fields.level(
            CA.RRTMGPI.array2field(FT.(atmos_sim.integrator.p.radiation_model.face_sw_flux_dn), face_space),
            half,
        )
        SWu_TOA = Fields.level(
            CA.RRTMGPI.array2field(FT.(atmos_sim.integrator.p.radiation_model.face_sw_flux_up), face_space),
            half,
        )

        return @. -(LWd_TOA + SWd_TOA - LWu_TOA - SWu_TOA) # [W/m^2]
    else
        return FT(0)
    end
end

get_field(sim::ClimaAtmosSimulation, ::Val{:surface_temperature}) =  TD.air_temperature.(thermo_params, sim.integrator.p.sfc_conditions.ts)

# additional BucketSimulation debug fields
get_field(sim::BucketSimulation, ::Val{:σS}) = sim.integrator.u.bucket.σS
get_field(sim::BucketSimulation, ::Val{:Ws}) = sim.integrator.u.bucket.Ws
get_field(sim::BucketSimulation, ::Val{:W}) = sim.integrator.u.bucket.W

# currently selected plot fields
plot_field_names(sim::SurfaceModelSimulation) = (:area_fraction, :surface_temperature, :surface_humidity)
plot_field_names(sim::BucketSimulation) =
    (:area_fraction, :surface_temperature, :surface_humidity, :air_density, :σS, :Ws, :W)
plot_field_names(sim::ClimaAtmosSimulation) = (:w, :ρq_tot, :ρe_tot, :liquid_precipitation, :snow_precipitation, :F_radiative_sfc, :surface_temperature)
