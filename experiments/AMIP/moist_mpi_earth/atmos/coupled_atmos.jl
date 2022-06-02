using ClimaCore.Geometry: ⊗

# These functions are copied from the ClimaAtmos driver and modified for the coupler (TODO: update when interface is stable)
function vertical_diffusion_boundary_layer_coupled_tendency!(Yₜ, Y, p, t)
    ᶜρ = Y.c.ρ
    (; ᶜp, ᶠv_a, ᶠz_a, ᶠK_E) = p # assume ᶜts and ᶜp have been updated
    (; dif_flux_energy, dif_flux_ρq_tot, dif_flux_uₕ) = p
    ᶠgradᵥ = Operators.GradientC2F() # apply BCs to ᶜdivᵥ, which wraps ᶠgradᵥ

    Fields.field_values(ᶠv_a) .= Fields.field_values(Spaces.level(Y.c.uₕ, 1)) .* one.(Fields.field_values(ᶠz_a)) # TODO: fix VIJFH copyto! to remove this
    @. ᶠK_E = eddy_diffusivity_coefficient(norm(ᶠv_a), ᶠz_a, ᶠinterp(ᶜp)) #* FT(0)

    # Total Energy
    if :ρe in propertynames(Y.c)
        ᶜdivᵥ = Operators.DivergenceF2C(
            top = Operators.SetValue(Geometry.WVector(FT(0))),
            bottom = Operators.SetValue(.-dif_flux_energy),
        )
        @. Yₜ.c.ρe += ᶜdivᵥ(ᶠK_E * ᶠinterp(ᶜρ) * ᶠgradᵥ((Y.c.ρe + ᶜp) / ᶜρ))
    end

    # Liquid Mass and Total Mass
    if :ρq_tot in propertynames(Y.c)
        ᶜdivᵥ = Operators.DivergenceF2C(
            top = Operators.SetValue(Geometry.WVector(FT(0))),
            bottom = Operators.SetValue(.-dif_flux_ρq_tot),
        )
        @. Yₜ.c.ρq_tot += ᶜdivᵥ(ᶠK_E * ᶠinterp(ᶜρ) * ᶠgradᵥ(Y.c.ρq_tot / ᶜρ))
        @. Yₜ.c.ρ += ᶜdivᵥ(ᶠK_E * ᶠinterp(ᶜρ) * ᶠgradᵥ(Y.c.ρq_tot / ᶜρ))
    end

    # Momentum
    if :uₕ in propertynames(Y.c)
        ᶜdivᵥ = Operators.DivergenceF2C(
            top = Operators.SetValue(Geometry.Contravariant3Vector(FT(0)) ⊗ Geometry.Covariant12Vector(FT(0), FT(0))),
            bottom = Operators.SetValue(.-dif_flux_uₕ),
        )

        @. Yₜ.c.uₕ += ᶜdivᵥ(ᶠK_E * ᶠgradᵥ(Y.c.uₕ))
    end

end

function vertical_diffusion_boundary_layer_coupled_cache(Y; Cd = FT(0.0014), Ch = FT(0.0014))
    ᶠz_a = similar(Y.f, FT)
    z_bottom = Spaces.level(Fields.coordinate_field(Y.c).z, 1)
    Fields.field_values(ᶠz_a) .= Fields.field_values(z_bottom) .* one.(Fields.field_values(ᶠz_a))
    # TODO: fix VIJFH copyto! to remove the one.(...)

    if :ρq_tot in propertynames(Y.c)
        dif_flux_ρq_tot = similar(z_bottom, Geometry.WVector{FT}) #ones(axes(z_bottom))
    else
        dif_flux_ρq_tot = Ref(Geometry.WVector(FT(0)))
    end

    dif_flux_uₕ =
        Geometry.Contravariant3Vector.(zeros(axes(z_bottom))) .⊗
        Geometry.Covariant12Vector.(zeros(axes(z_bottom)), zeros(axes(z_bottom)))

    if (:ρq_liq in propertynames(Y.c) && :ρq_ice in propertynames(Y.c) && :ρq_tot in propertynames(Y.c))
        ts_type = TD.PhaseNonEquil{FT}
    elseif :ρq_tot in propertynames(Y.c)
        ts_type = TD.PhaseEquil{FT}
    else
        ts_type = TD.PhaseDry{FT}
    end
    coef_type = SF.Coefficients{
        FT,
        SF.InteriorValues{FT, Tuple{FT, FT}, ts_type},
        SF.SurfaceValues{FT, Tuple{FT, FT}, TD.PhaseEquil{FT}},
    }

    return (;
        ᶠv_a = similar(Y.f, eltype(Y.c.uₕ)),
        ᶠz_a,
        ᶠK_E = similar(Y.f, FT),
        flux_coefficients = similar(z_bottom, coef_type),
        dif_flux_energy = similar(z_bottom, Geometry.WVector{FT}),
        dif_flux_ρq_tot,
        dif_flux_uₕ,
        ∂F_aero∂T_sfc = zeros(axes(z_bottom)),
        Cd,
        Ch,
    )
end

function rrtmgp_model_coupled_cache(
    Y,
    params;
    radiation_mode = ClearSkyRadiation(),
    interpolation = BestFit(),
    bottom_extrapolation = SameAsInterpolation(),
    idealized_insolation = true,
    idealized_h2o = false,
)
    bottom_coords = Fields.coordinate_field(Spaces.level(Y.c, 1))
    if eltype(bottom_coords) <: Geometry.LatLongZPoint
        latitude = field2array(bottom_coords.lat)
    else
        latitude = field2array(zero(bottom_coords.z)) # flat space is on Equator
    end
    input_data = rrtmgp_artifact("atmos_state", "clearsky_as.nc")
    if radiation_mode isa GrayRadiation
        kwargs = (;
            lapse_rate = 3.5,
            optical_thickness_parameter = (@. (
                (300 + 60 * (FT(1 / 3) - sind(latitude)^2)) / 200
            )^4 - 1),
        )
    else
        # the pressure and ozone concentrations are provided for each of 100
        # sites, which we average across
        n = input_data.dim["layer"]
        input_center_pressure =
            vec(mean(reshape(input_data["pres_layer"][:, :], n, :); dims = 2))
        # the first values along the third dimension of the ozone concentration
        # data are the present-day values
        input_center_volume_mixing_ratio_o3 =
            vec(mean(reshape(input_data["ozone"][:, :, 1], n, :); dims = 2))

        # interpolate the ozone concentrations to our initial pressures (set the
        # kinetic energy to 0 when computing the pressure using total energy)
        pressure2ozone =
            Spline1D(input_center_pressure, input_center_volume_mixing_ratio_o3)
        if :ρθ in propertynames(Y.c)
            ᶜts = @. thermo_state_ρθ(Y.c.ρθ, Y.c, params)
        elseif :ρe in propertynames(Y.c)
            ᶜΦ = FT(Planet.grav(params)) .* Fields.coordinate_field(Y.c).z
            ᶜts = @. thermo_state_ρe(Y.c.ρe, Y.c, 0, ᶜΦ, params)
        elseif :ρe_int in propertynames(Y.c)
            ᶜts = @. thermo_state_ρe_int(Y.c.ρe_int, Y.c, params)
        end
        ᶜp = @. TD.air_pressure(params, ᶜts)
        center_volume_mixing_ratio_o3 = field2array(@. FT(pressure2ozone(ᶜp)))

        # the first value for each global mean volume mixing ratio is the
        # present-day value
        input_vmr(name) =
            input_data[name][1] * parse(FT, input_data[name].attrib["units"])
        kwargs = (;
            use_global_means_for_well_mixed_gases = true,
            center_volume_mixing_ratio_h2o = NaN, # initialize in tendency
            center_volume_mixing_ratio_o3,
            volume_mixing_ratio_co2 = input_vmr("carbon_dioxide_GM"),
            volume_mixing_ratio_n2o = input_vmr("nitrous_oxide_GM"),
            volume_mixing_ratio_co = input_vmr("carbon_monoxide_GM"),
            volume_mixing_ratio_ch4 = input_vmr("methane_GM"),
            volume_mixing_ratio_o2 = input_vmr("oxygen_GM"),
            volume_mixing_ratio_n2 = input_vmr("nitrogen_GM"),
            volume_mixing_ratio_ccl4 = input_vmr("carbon_tetrachloride_GM"),
            volume_mixing_ratio_cfc11 = input_vmr("cfc11_GM"),
            volume_mixing_ratio_cfc12 = input_vmr("cfc12_GM"),
            volume_mixing_ratio_cfc22 = input_vmr("hcfc22_GM"),
            volume_mixing_ratio_hfc143a = input_vmr("hfc143a_GM"),
            volume_mixing_ratio_hfc125 = input_vmr("hfc125_GM"),
            volume_mixing_ratio_hfc23 = input_vmr("hfc23_GM"),
            volume_mixing_ratio_hfc32 = input_vmr("hfc32_GM"),
            volume_mixing_ratio_hfc134a = input_vmr("hfc134a_GM"),
            volume_mixing_ratio_cf4 = input_vmr("cf4_GM"),
            volume_mixing_ratio_no2 = 1e-8, # not available in input_data
            latitude,
        )
        if !(radiation_mode isa ClearSkyRadiation)
            error("rrtmgp_model_cache not yet implemented for $radiation_mode")
        end
    end

    if requires_z(interpolation) || requires_z(bottom_extrapolation)
        kwargs = (;
            kwargs...,
            center_z = field2array(Fields.coordinate_field(Y.c).z),
            face_z = field2array(Fields.coordinate_field(Y.f).z),
        )
    end

    if idealized_insolation
        # perpetual equinox with no diurnal cycle
        solar_zenith_angle = FT(π) / 3
        weighted_irradiance =
            @. 1360 * (1 + FT(1.2) / 4 * (1 - 3 * sind(latitude)^2)) /
               (4 * cos(solar_zenith_angle))
    else
        solar_zenith_angle = weighted_irradiance = NaN # initialized in tendency
    end

    if idealized_h2o && radiation_mode isa GrayRadiation
        error("idealized_h2o cannot be used with GrayRadiation")
    end

    # surface_emissivity and surface_albedo are provided for each of 100 sites,
    # which we average across
    rrtmgp_model = RRTMGPModel(
        params;
        FT = Float64,
        ncol = length(Spaces.all_nodes(axes(Spaces.level(Y.c, 1)))),
        domain_nlay = Spaces.nlevels(axes(Y.c)),
        radiation_mode,
        interpolation,
        bottom_extrapolation,
        add_isothermal_boundary_layer = true,
        center_pressure = NaN, # initialized in tendency
        center_temperature = NaN, # initialized in tendency
        surface_temperature = (@. 29 * exp(-(latitude / 26)^2 / 2) + 271),
        surface_emissivity = mean(input_data["surface_emissivity"]),
        direct_sw_surface_albedo = mean(input_data["surface_albedo"]),
        diffuse_sw_surface_albedo = mean(input_data["surface_albedo"]),
        solar_zenith_angle,
        weighted_irradiance,
        kwargs...,
    )
    close(input_data)
    return (;
        ᶜT = similar(Y.c, FT),
        ᶜvmr_h2o = similar(Y.c, FT),
        insolation_tuple = similar(Spaces.level(Y.c, 1), Tuple{FT, FT, FT}),
        zenith_angle = similar(Spaces.level(Y.c, 1), FT),
        weighted_irradiance = similar(Spaces.level(Y.c, 1), FT),
        ᶠradiation_flux = similar(Y.f, Geometry.WVector{FT}),
        idealized_insolation,
        idealized_h2o,
        rrtmgp_model,
    )
end

# TODO: flip order so that NamedTuple() is fallback.
additional_cache(Y, params, dt; use_tempest_mode = false) = merge(
    hyperdiffusion_cache(Y; κ₄ = FT(2e17), use_tempest_mode),
    sponge ? rayleigh_sponge_cache(Y, dt) : NamedTuple(),
    isnothing(microphy) ? NamedTuple() : zero_moment_microphysics_cache(Y),
    isnothing(forcing) ? NamedTuple() : held_suarez_cache(Y),
    isnothing(rad) ? NamedTuple() : rrtmgp_model_coupled_cache(Y, params; radiation_mode, idealized_h2o),
    vert_diff ? vertical_diffusion_boundary_layer_coupled_cache(Y) : NamedTuple(),
    (;
        tendency_knobs = (;
            hs_forcing = forcing == "held_suarez",
            microphy_0M = microphy == "0M",
            rad_flux = !isnothing(rad),
            vert_diff,
            hyperdiff,
        )
    ),
)

additional_tendency!(Yₜ, Y, p, t) = begin
    (; rad_flux, vert_diff, hs_forcing) = p.tendency_knobs
    (; microphy_0M, hyperdiff) = p.tendency_knobs
    hyperdiff && hyperdiffusion_tendency!(Yₜ, Y, p, t)
    sponge && rayleigh_sponge_tendency!(Yₜ, Y, p, t)
    hs_forcing && held_suarez_tendency!(Yₜ, Y, p, t)
    vert_diff && vertical_diffusion_boundary_layer_coupled_tendency!(Yₜ, Y, p, t)
    microphy_0M && zero_moment_microphysics_tendency!(Yₜ, Y, p, t)
    rad_flux && rrtmgp_model_tendency!(Yₜ, Y, p, t)
end

parsed_args["microphy"] = "0M"
parsed_args["forcing"] = nothing
parsed_args["idealized_h2o"] = false
parsed_args["vert_diff"] = true
parsed_args["rad"] = "gray"
parsed_args["hyperdiff"] = true
parsed_args["config"] = "sphere"
parsed_args["moist"] = "equil"
