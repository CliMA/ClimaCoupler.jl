include("coupler_driver_init.jl")
include("coupler_parse_args.jl")

include("coupler_component_init.jl")
cs = ClimaCoupler.Interfacer.CoupledSimulation{FT}(
    comms_ctx,
    dates,
    boundary_space,
    coupler_fields,
    config_dict,
    conservation_checks,
    [tspan[1], tspan[2]],
    atmos_sim.integrator.t,
    Δt_cpl,
    (; land = land_fraction, ocean = zeros(boundary_space), ice = zeros(boundary_space)),
    model_sims,
    mode_specifics,
    diagnostics,
    callbacks,
    dir_paths,
);

turbulent_fluxes = nothing
if config_dict["turb_flux_partition"] == "PartitionedStateFluxes"
    turbulent_fluxes = PartitionedStateFluxes()
elseif config_dict["turb_flux_partition"] == "CombinedStateFluxes"
    turbulent_fluxes = CombinedStateFluxes()
else
    error("turb_flux_partition must be either PartitionedStateFluxes or CombinedStateFluxes")
end

# 1) coupler combines surface states and calculates rho_sfc using surface and atmos variables
update_surface_fractions!(cs)
import_combined_surface_fields!(cs.fields, cs.model_sims, cs.boundary_space, turbulent_fluxes)
import_atmos_fields!(cs.fields, cs.model_sims, cs.boundary_space, turbulent_fluxes)
update_model_sims!(cs.model_sims, cs.fields, turbulent_fluxes)

# 2) each surface component model calculates its own vapor specific humidity (q_sfc)
# TODO: the q_sfc calculation follows the design of the bucket q_sfc, but it would be neater to abstract this from step!
step!(land_sim, Δt_cpl)
step!(ocean_sim, Δt_cpl)
step!(ice_sim, Δt_cpl)

# 3) coupler re-imports updated surface fields and calculates turbulent fluxes, while updating atmos sfc_conditions
if turbulent_fluxes isa CombinedStateFluxes
    # calculate fluxes using combined surface states on the atmos grid
    import_combined_surface_fields!(cs.fields, cs.model_sims, cs.boundary_space, turbulent_fluxes) # i.e. T_sfc, albedo, z0, beta, q_sfc
    combined_turbulent_fluxes!(cs.model_sims, cs.fields, turbulent_fluxes) # this updates the atmos thermo state, sfc_ts
elseif turbulent_fluxes isa PartitionedStateFluxes
    # calculate turbulent fluxes in surface models and save the weighted average in coupler fields
    partitioned_turbulent_fluxes!(cs.model_sims, cs.fields, cs.boundary_space, MoninObukhovScheme(), thermo_params)

    # update atmos sfc_conditions for surface temperature
    # TODO: this is hard coded and needs to be simplified (need CA modification)
    new_p = get_new_cache(atmos_sim, cs.fields)
    CA.SurfaceConditions.update_surface_conditions!(atmos_sim.integrator.u, new_p, atmos_sim.integrator.t) # sets T_sfc (but SF calculation not necessary - CA)
    atmos_sim.integrator.p.precomputed.sfc_conditions .= new_p.precomputed.sfc_conditions
end

# 4) given the new sfc_conditions, atmos calls the radiative flux callback
reinit_model_sims!(cs.model_sims) # NB: for atmos this sets a nonzero radiation flux

# 5) coupler re-imports updated atmos fluxes (radiative fluxes for both `turbulent_fluxes` types
# and also turbulent fluxes if `turbulent_fluxes isa CombinedStateFluxes`,
# and sends them to the surface component models. If `turbulent_fluxes isa PartitionedStateFluxes`
# atmos receives the turbulent fluxes from the coupler.
import_atmos_fields!(cs.fields, cs.model_sims, cs.boundary_space, turbulent_fluxes)
update_model_sims!(cs.model_sims, cs.fields, turbulent_fluxes)

#=
## Coupling Loop
=#
function solve_coupler!(cs::ClimaCoupler.Interfacer.CoupledSimulation)
    @info "Starting coupling loop"

    (; model_sims, Δt_cpl, tspan) = cs
    (; atmos_sim, land_sim, ocean_sim, ice_sim) = model_sims

    ## step in time
    walltime = @elapsed for t in ((tspan[1] + Δt_cpl):Δt_cpl:tspan[end])

        cs.dates.date[1] = current_date(cs, t) # if not global, `date` is not updated.

        ## print date on the first of month
        if cs.dates.date[1] >= cs.dates.date1[1]
            ClimaComms.iamroot(comms_ctx) ? @show(cs.dates.date[1]) : nothing
        end

        if cs.mode.name == "amip"

            ## monthly read of boundary condition data for SST and SIC and CO2
            if cs.dates.date[1] >= next_date_in_file(cs.mode.SST_info)
                update_midmonth_data!(cs.dates.date[1], cs.mode.SST_info)
            end
            SST_current = interpolate_midmonth_to_daily(cs.dates.date[1], cs.mode.SST_info)
            update_field!(ocean_sim, Val(:surface_temperature), SST_current)

            if cs.dates.date[1] >= next_date_in_file(cs.mode.SIC_info)
                update_midmonth_data!(cs.dates.date[1], cs.mode.SIC_info)
            end
            SIC_current =
                get_ice_fraction.(interpolate_midmonth_to_daily(cs.dates.date[1], cs.mode.SIC_info), mono_surface)
            update_field!(ice_sim, Val(:area_fraction), SIC_current)

            if cs.dates.date[1] >= next_date_in_file(cs.mode.CO2_info)
                update_midmonth_data!(cs.dates.date[1], cs.mode.CO2_info)
            end
            CO2_current = interpolate_midmonth_to_daily(cs.dates.date[1], cs.mode.CO2_info)
            update_field!(atmos_sim, Val(:co2_gm), CO2_current)

            ## calculate and accumulate diagnostics at each timestep
            ClimaComms.barrier(comms_ctx)
            accumulate_diagnostics!(cs)

            ## save and reset monthly averages
            save_diagnostics(cs)

        end

        ## compute global energy
        !isnothing(cs.conservation_checks) ? check_conservation!(cs) : nothing

        ## run component models sequentially for one coupling timestep (Δt_cpl)
        ClimaComms.barrier(comms_ctx)
        update_surface_fractions!(cs)
        update_model_sims!(cs.model_sims, cs.fields, turbulent_fluxes)

        ## step sims
        step_model_sims!(cs.model_sims, t)

        ## exchange combined fields and (if specified) calculate fluxes using combined states
        import_combined_surface_fields!(cs.fields, cs.model_sims, cs.boundary_space, turbulent_fluxes) # i.e. T_sfc, albedo, z0, beta
        if turbulent_fluxes isa CombinedStateFluxes
            combined_turbulent_fluxes!(cs.model_sims, cs.fields, turbulent_fluxes) # this updates the surface thermo state, sfc_ts, in ClimaAtmos (but also unnecessarily calculates fluxes)
        elseif turbulent_fluxes isa PartitionedStateFluxes
            # calculate turbulent fluxes in surfaces and save the weighted average in coupler fields
            partitioned_turbulent_fluxes!(cs.model_sims, cs.fields, cs.boundary_space, MoninObukhovScheme(), thermo_params)

            # update atmos sfc_conditions for surface temperature - TODO: this needs to be simplified (need CA modification)
            new_p = get_new_cache(atmos_sim, cs.fields)
            CA.SurfaceConditions.update_surface_conditions!(atmos_sim.integrator.u, new_p, atmos_sim.integrator.t) # to set T_sfc (but SF calculation not necessary - CA modification)
            atmos_sim.integrator.p.precomputed.sfc_conditions .= new_p.precomputed.sfc_conditions
        end

        import_atmos_fields!(cs.fields, cs.model_sims, cs.boundary_space, turbulent_fluxes) # radiative and/or turbulent

        ## callback to update the fist day of month if needed (for BCReader)
        trigger_callback!(cs, cs.callbacks.update_firstdayofmonth!)

        ## callback to checkpoint model state
        trigger_callback!(cs, cs.callbacks.checkpoint)

    end
    @show walltime

    return cs
end
