"""
    calculate_and_send_turbulent_fluxes!(cs::)

Calculates surface fluxes using adapter function `get_surface_fluxes_point!`
from ClimaAtmos that calls `SurfaceFluxes.jl`. The coupler updates in
atmos model cache fluxes at each coupling timestep.

- TODO: generalize interface for regridding and take land state out of atmos's integrator.p

The current setup calculates the aerodynamic fluxes in the coupler and assumes no regridding is needed.
(NB: Radiation surface fluxes are calculated by the atmosphere.)

"""

function calculate_and_send_turbulent_fluxes!(model_sims, fields, boundary_space, surface_scheme, thermo_params)

    atmos_sim = model_sims.atmos_sim;
    csf = fields
    FT = eltype(csf[1])

    # reset coupler fields (TODO: add flux accumulation)
    csf.F_ρτxz .*= FT(0)
    csf.F_ρτyz .*= FT(0)
    csf.F_shf .*= FT(0)
    csf.F_lhf .*= FT(0)
    csf.F_evap .*= FT(0)

    # for sim in model_sims
    #     if sim isa SurfaceModelSimulation
    #         # primarily to allow rho_sfc calculation
    #         # TODO: absorb in to colidx once aux_update in ClimaLSM allows this
    #         extra_aux_update(sim, thermo_params, get_field(atmos_sim, Val(:thermo_state_int)) )
    #     end
    # end

    # iterate over all columns (when regridding, this will need to change)
    Fields.bycolumn(boundary_space) do colidx
        # atmos state of center level 1
        z_int = get_field(atmos_sim, Val(:height_int), colidx)
        uₕ_int = get_field(atmos_sim, Val(:uv_int), colidx)
        thermo_state_int = get_field(atmos_sim, Val(:thermo_state_int), colidx)

        z_sfc = get_field(atmos_sim, Val(:height_sfc), colidx)

        for sim in model_sims
            # iterate over all surface models with non-zero area fractions
            if sim isa SurfaceModelSimulation
                area_fraction = get_field(sim, Val(:area_fraction), colidx)
                area_mask  = Regridder.binary_mask.(area_fraction, threshold = eps())

                if !iszero(parent(area_mask))

                    thermo_state_sfc = surface_thermo_state(sim, thermo_params, thermo_state_int, colidx)

                    # set inputs based on whether the surface_scheme is MOST or bulk
                    inputs = surface_inputs(
                        surface_scheme,
                        thermo_state_sfc,
                        thermo_state_int,
                        uₕ_int,
                        z_int,
                        z_sfc,
                        get_scheme_specific_properties(surface_scheme, sim, colidx)...,
                        )

                    # update fluxes in the coupler
                    surface_params = get_surface_params(atmos_sim)
                    F_ρτxz, F_ρτyz, F_shf, F_lhf, F_evap = get_surface_fluxes_point!(inputs, surface_params)

                    fields = (; F_ρτxz = F_ρτxz, F_ρτyz = F_ρτyz, F_shf = F_shf, F_lhf = F_lhf, F_evap = F_evap)

                    # update the fluxes of this surface model
                    update_turbulent_fluxes_point!(sim, fields, colidx)

                    # add the flux contributing from this surface
                    @. csf.F_ρτxz[colidx] += F_ρτxz * area_fraction * area_mask
                    @. csf.F_ρτyz[colidx] += F_ρτyz * area_fraction * area_mask
                    @. csf.F_shf[colidx] += F_shf * area_fraction * area_mask
                    @. csf.F_lhf[colidx] += F_lhf * area_fraction * area_mask
                    @. csf.F_evap[colidx] += F_evap * area_fraction * area_mask
                    # @. csf.ρ_sfc[colidx] += ρ_sfc
                    # @. csf.q_sfc[colidx] += q_sfc

                end
            end
        end

    end

    # update atmos fluxes (TODO: include to the above loop, with atmos_flux_reset)
    for sim in model_sims
        if sim isa AtmosModelSimulation
            Fields.bycolumn(boundary_space) do colidx
                coupler_fields = (; F_ρτxz = csf.F_ρτxz[colidx], F_ρτyz = csf.F_ρτyz[colidx], F_shf = csf.F_shf[colidx], F_lhf = csf.F_lhf[colidx], F_evap = csf.F_evap[colidx])
                update_turbulent_fluxes_point!(sim, coupler_fields, colidx)
            end
        end
    end
    # TODO: add allowable bounds here, check explicitly that all fluxes are equal

    # check_field = zeros(boundary_space)
    # for sim in model_sims
    #     if sim isa SurfaceModelSimulation
    #         check_field .+= get_sensible_heat_flux(sim)
    #     end
    # end
    # @assert(extrema(check_field .- get_sensible_heat_flux(atmos_sim)) ≈ (0.0, 0.0))

end

function get_scheme_specific_properties(::BulkScheme, sim, colidx)
    Ch = get_field(sim, Val(:heat_transfer_coefficient), colidx)
    Cd = get_field(sim, Val(:beta), colidx)
    beta = get_field(sim, Val(:drag_coefficient), colidx)
    FT = eltype(Ch)
    return (; z0b = FT(0), z0m = FT(0), Ch = Ch, Cd = Cd, beta = beta, gustiness = FT(1))
end
function surface_inputs(
    ::BulkScheme,
    thermo_state_sfc,
    thermo_state_int,
    uₕ_int,
    z_int,
    z_sfc,
    z0b,
    z0m,
    Ch,
    Cd,
    beta,
    gustiness,
)
    FT = Spaces.undertype(axes(z_sfc))

    # wrap state values
    return @. SF.Coefficients(
        SF.InteriorValues(z_int, uₕ_int, thermo_state_int), # state_in                              # state_sfc
        SF.SurfaceValues(                                  # state_sfc
            z_sfc,
            StaticArrays.SVector(FT(0), FT(0)),
            thermo_state_sfc,
        ),
        Cd,                                     # Cd
        Ch,                                     # Ch
        z0m,                                             # z0m
        z0b,                                             # z0b
        gustiness,                                             # gustiness
        beta,                                   # beta
    )
end

function get_scheme_specific_properties(::MoninObukhovScheme, sim, colidx)
    z0m = get_field(sim, Val(:z0m), colidx)
    z0b = get_field(sim, Val(:z0b), colidx)
    beta = get_field(sim, Val(:drag_coefficient), colidx)
    return (; z0b = z0b, z0m = z0m, Ch = FT(0), Cd = FT(0), beta= beta, gustiness = FT(1))
end

function surface_inputs(
    ::MoninObukhovScheme,
    thermo_state_sfc,
    thermo_state_int,
    uₕ_int,
    z_int,
    z_sfc,
    z0b,
    z0m,
    Ch,
    Cd,
    beta,
    gustiness,
)
    FT = Spaces.undertype(axes(z_sfc))

    # wrap state values

    return @. SF.ValuesOnly(
        SF.InteriorValues(z_int, uₕ_int, thermo_state_int), # state_in
        SF.SurfaceValues(                                  # state_sfc
            z_sfc,
            StaticArrays.SVector(FT(0), FT(0)),
            thermo_state_sfc,
        ),
        z0m,                                    # z0m
        z0b,                                    # z0b
        FT(-1),                                            # L_MO_init
        gustiness,                                             # gustiness
        beta                                   # beta
    )

end


function surface_thermo_state(sim::SurfaceModelSimulation, thermo_params, thermo_state_int, colidx)
    @warn("Simulation " * name(sim) * " uses the default thermo (saturated) surface state", maxlog = 10)
    T_sfc = get_field(sim, Val(:air_temperature), colidx) #
    ρ_sfc = extrapolate_ρ_to_sfc.(thermo_params, thermo_state_int, T_sfc) # ideally the # calculate elsewhere, here just getter...
    q_sfc = TD.q_vap_saturation_generic.(thermo_params, T_sfc, ρ_sfc, TD.Liquid()) # default = saturated
    @. TD.PhaseEquil_ρTq.(thermo_params, ρ_sfc, T_sfc, q_sfc)
end

