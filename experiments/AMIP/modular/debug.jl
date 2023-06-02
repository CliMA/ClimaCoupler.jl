using Plots
using ClimaCorePlots


function debug(sim::Union{SlabOceanSimulation,PrescribedIceSimulation})
    @info "debug" * name(sim)
    map(x -> println(zip(x, getproperty(sim.integrator.u, x))), propertynames(sim.integrator.u))
    @info "    -> state extrema"
    map(x -> println(zip(x, extrema(getproperty(sim.integrator.u, x)))), propertynames(sim.integrator.u))
    @info "    -> cache extrema"
    # map(x -> println(zip(x, extrema(getproperty(sim.integrator.p, x)))), propertynames(sim.integrator.p))
    println(zip(:F_aero, extrema(sim.integrator.p.F_aero)))
    println(zip(:F_rad, extrema(sim.integrator.p.F_rad)))

    all_plots = []
    push!(all_plots, Plots.plot(sim.integrator.u.T_sfc, title = "T_sfc"))
    push!(all_plots, Plots.plot(sim.integrator.p.F_aero, title = "F_aero"))
    push!(all_plots, Plots.plot(sim.integrator.p.F_rad, title = "F_rad"))
    fig = Plots.plot(all_plots..., size = (1500, 800))
    Plots.png("debug_" * name(sim))
end


using Plots
using ClimaCorePlots
function debug(sim::ClimaAtmosSimulation)
    @info "debug" * name(sim)
    map(x -> println(zip(x, getproperty(sim.integrator.u, x))), propertynames(sim.integrator.u))
    @info "    -> state extrema"
    #map(x -> println(zip(x, extrema(getproperty(sim.integrator.u, x)))), propertynames(sim.integrator.u))
    println(zip(:ρe_tot, extrema(sim.integrator.u.c.ρe_tot)))
    println(zip(:ρq_tot, extrema(sim.integrator.u.c.ρq_tot)))
    @info "    -> cache extrema"
    println(zip(:dif_flux_energy_bc, extrema(sim.integrator.p.dif_flux_energy_bc.components.data.:1)))
    println(zip(:dif_flux_uₕ_bc, extrema(sim.integrator.p.dif_flux_uₕ_bc.components.data.:1)))
    println(zip(:dif_flux_ρq_tot_bc, extrema(sim.integrator.p.dif_flux_ρq_tot_bc.components.data.:1)))
    println(zip(:col_integrated_rain, extrema(sim.integrator.p.col_integrated_rain)))

    all_plots = []
    push!(all_plots, Plots.plot(sim.integrator.u.c.ρq_tot, title = "ρq_tot"))
    push!(all_plots, Plots.plot(sim.integrator.p.col_integrated_rain, title = "col_integrated_rain"))
    push!(all_plots, Plots.plot(sim.integrator.p.dif_flux_energy_bc.components.data.:1, title = "dif_flux_energy_bc"))
    push!(all_plots, Plots.plot(sim.integrator.p.dif_flux_ρq_tot_bc.components.data.:1, title = "dif_flux_ρq_tot_bc"))
    fig = Plots.plot(all_plots..., size = (1500, 800))
    Plots.png("debug_" * name(sim))

    all_plots = []
    push!(all_plots, Plots.plot(sim.Y_init.c.ρq_tot, title = "ρq_tot"))
    push!(all_plots, Plots.plot(sim.Y_init.c.ρe_tot, title = "ρe_tot"))
    fig = Plots.plot(all_plots..., size = (1500, 400))
    Plots.png("debug_init_" * name(sim))
end


function debug(cs)
    map(x -> println(zip(x, extrema(getproperty(cs.fields, x)))), propertynames(cs.fields))
    all_plots = []
    push!(all_plots, Plots.plot(cs.fields.F_lhf, title = "F_lhf"))
    push!(all_plots, Plots.plot(cs.fields.F_shf, title = "F_shf"))
    push!(all_plots, Plots.plot(cs.fields.P_liq, title = "P_liq"))
    push!(all_plots, Plots.plot(cs.fields.T_S, title = "T_S"))
    push!(all_plots, Plots.plot(cs.surface_fractions.land, title = "fractions.land"))
    push!(all_plots, Plots.plot(cs.surface_fractions.ocean, title = "cs.surface_fractions.ocean"))
    push!(all_plots, Plots.plot(cs.surface_fractions.ice, title = "cs.surface_fractions.ice"))
    fig = Plots.plot(all_plots..., size = (1500, 800))
    Plots.png("debug_coupler")
end


function debug(sim::BucketSimulation)
    @info "debug" * name(sim)
    map(x -> println(zip(x, getproperty(sim.integrator.u.bucket, x))), propertynames(sim.integrator.u.bucket))
    @info "    -> state extrema"
    map(x -> println(zip(x, extrema(getproperty(sim.integrator.u.bucket, x)))), propertynames(sim.integrator.u.bucket))
    @info "    -> cache extrema"
    map(x -> println(zip(x, extrema(getproperty(sim.integrator.p.bucket, x)))), propertynames(sim.integrator.p.bucket))

    all_plots = []
    push!(all_plots, Plots.plot(sim.integrator.u.bucket.T, title = "T"))
    push!(all_plots, Plots.plot(sim.integrator.u.bucket.W, title = "W"))
    push!(all_plots, Plots.plot(sim.integrator.u.bucket.σS, title = "σS"))
    push!(all_plots, Plots.plot(sim.integrator.u.bucket.Ws, title = "Ws"))
    push!(all_plots, Plots.plot(sim.integrator.p.bucket.q_sfc, title = "q_sfc"))
    push!(all_plots, Plots.plot(sim.integrator.p.bucket.T_sfc, title = "T_sfc"))
    push!(all_plots, Plots.plot(sim.integrator.p.bucket.α_sfc, title = "α_sfc"))
    push!(all_plots, Plots.plot(sim.integrator.p.bucket.ρ_sfc, title = "ρ_sfc"))
    push!(all_plots, Plots.plot(sim.integrator.p.bucket.P_snow, title = "P_snow"))
    push!(all_plots, Plots.plot(sim.integrator.p.bucket.P_liq, title = "P_liq"))
    fig = Plots.plot(all_plots..., size = (1500, 800))
    Plots.png("debug_" * name(sim))
end

debug(atmos_sim)
debug(land_sim)
cs.mode.name == "amip" ? debug(ice_sim) : nothing
cs.mode.name == "slabplanet" ? debug(ocean_sim) : nothing
debug(cs)