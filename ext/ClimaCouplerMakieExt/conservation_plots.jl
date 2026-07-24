import Makie
import ClimaCoupler: Interfacer, ConservationChecker, Plotting

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

Conservation checks are available for energy and water, and can be enabled by
running a Slabplanet simulation with `energy_check` set to true.

If `softfail` is false, asserts that the relative error in conservation of the
provided quantity is less than a pre-determined threshold. This argument
is controlled by the `conservation_softfail` simulation flag.
"""
function Plotting.plot_global_conservation(
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

    var_name = nameof(cc)
    cum_total = [0.0]
    f = Makie.Figure()
    ax = Makie.Axis(f[1, 1], xlabel = "time [days]", ylabel = "$var_name: (t) - (t=0)")
    Makie.lines!(ax, days, total .- total[1], label = "total"; linewidth = 3)
    for sim in model_sims
        sim_name = nameof(sim)
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
    l_rse_valid = filter(x -> !isinf(x) && !isnan(x), l_rse)
    if !isempty(l_rse_valid)
        y_min = minimum(l_rse_valid)
        y_max = maximum(l_rse_valid)
        if y_min != y_max
            Makie.ylims!(y_min, y_max)
        else
            # If all values are the same, add a small padding to avoid Makie error
            padding = max(abs(y_min) * 0.01, 0.1)
            Makie.ylims!(y_min - padding, y_max + padding)
        end
    end
    Makie.axislegend(position = :lt)
    Makie.save(figname2, lp)

    # check that the relative error is small (TODO: reduce this to sqrt(eps(FT)))
    if !softfail
        @info typeof(cc)
        @info rse[end]
        @assert rse[end] < 0.055
    end
end
