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

# reset land
function reset_land!(land_sim; W = 0, Ws = 0, σS = 0, T = 310, P_liq = 0, P_snow = 0, R_n = 0, evaporation = 0, T_sfc = 310, q_sfc = 0, turbulent_energy_flux = 0, α_sfc = 0.5, ρ_sfc = 1)
    land_sim.integrator.u.bucket.W .= W
    land_sim.integrator.u.bucket.Ws .= Ws
    land_sim.integrator.u.bucket.σS .=  σS
    land_sim.integrator.u.bucket.T .=  T

    land_sim.integrator.p.bucket.P_liq .= P_liq # 0.001 :( - increases W
    land_sim.integrator.p.bucket.P_snow .= P_snow # 0.001 :( - increases σS
    land_sim.integrator.p.bucket.R_n .= R_n # :) 100 - decreases T
    land_sim.integrator.p.bucket.evaporation .= evaporation # :| 0.001 decreases W (how can this be negative? - I guess this should be tapered by beta? Negativity checks)
    land_sim.integrator.p.bucket.T_sfc .= T_sfc
    land_sim.integrator.p.bucket.q_sfc .= q_sfc
    land_sim.integrator.p.bucket.turbulent_energy_flux .= turbulent_energy_flux # :) 100 - decreases T
    land_sim.integrator.p.bucket.α_sfc .= α_sfc
    land_sim.integrator.p.bucket.ρ_sfc .= ρ_sfc
end


# init land
reset_land!(land_sim; W=0.0, σS = 0.01)
debug(cs)

# flux calculation
land_sim.integrator.p.bucket.evaporation .= 0.001 .* ClimaLSM.surface_evaporative_scaling(land_sim.model, land_sim.integrator.u, land_sim.integrator.p)

# step land
step!(land_sim, land_sim.integrator.t + Δt_cpl)

debug(cs)
=#