
# get the paths to the necessary data files: land-sea mask, sst map, sea ice concentration
include(joinpath(pkgdir(ClimaCoupler), "artifacts", "artifact_funcs.jl"))
sst_data = joinpath(sst_dataset_path(), "sst.nc")
sic_data = joinpath(sic_dataset_path(), "sic.nc")
co2_data = joinpath(co2_dataset_path(), "mauna_loa_co2.nc")
land_mask_data = joinpath(mask_dataset_path(), "seamask.nc")

atmos_sim = atmos_init(FT, config_dict_atmos);
thermo_params = get_thermo_params(atmos_sim) # TODO: this should be shared by all models

#=
We use a common `Space` for all global surfaces. This enables the MPI processes to operate on the same columns in both
the atmospheric and surface components, so exchanges are parallelized. Note this is only possible when the
atmosphere and surface are of the same horizontal resolution.
=#
## init a 2D boundary space at the surface
boundary_space = Spaces.horizontal_space(atmos_sim.domain.face_space)

# init land-sea fraction
land_fraction =
    FT.(
        Regridder.land_fraction(
            FT,
            REGRID_DIR,
            comms_ctx,
            land_mask_data,
            "LSMASK",
            boundary_space,
            mono = mono_surface,
        )
    )

@info mode_name
if mode_name == "amip"
    @info "AMIP boundary conditions - do not expect energy conservation"

    ## land
    land_sim = bucket_init(
        FT,
        tspan,
        config_dict["land_domain_type"],
        config_dict["land_albedo_type"],
        config_dict["land_temperature_anomaly"],
        comms_ctx,
        REGRID_DIR;
        dt = Δt_cpl,
        space = boundary_space,
        saveat = saveat,
        area_fraction = land_fraction,
        date_ref = date0,
        t_start = t_start,
    )

    ## ocean
    SST_info = bcfile_info_init(
        FT,
        REGRID_DIR,
        sst_data,
        "SST",
        boundary_space,
        comms_ctx,
        interpolate_daily = true,
        scaling_function = clean_sst, ## convert to Kelvin
        land_fraction = land_fraction,
        date0 = date0,
        mono = mono_surface,
    )

    update_midmonth_data!(date0, SST_info)
    SST_init = interpolate_midmonth_to_daily(date0, SST_info)
    ocean_sim = SurfaceStub((;
        T_sfc = SST_init,
        ρ_sfc = ClimaCore.Fields.zeros(boundary_space),
        z0m = FT(1e-3),
        z0b = FT(1e-3),
        beta = FT(1),
        α = FT(0.06),
        area_fraction = (FT(1) .- land_fraction),
        phase = TD.Liquid(),
        thermo_params = thermo_params,
    ))

    ## sea ice
    SIC_info = bcfile_info_init(
        FT,
        REGRID_DIR,
        sic_data,
        "SEAICE",
        boundary_space,
        comms_ctx,
        interpolate_daily = true,
        scaling_function = clean_sic, ## convert to fraction
        land_fraction = land_fraction,
        date0 = date0,
        mono = mono_surface,
    )
    update_midmonth_data!(date0, SIC_info)
    SIC_init = interpolate_midmonth_to_daily(date0, SIC_info)
    ice_fraction = get_ice_fraction.(SIC_init, mono_surface)
    ice_sim = ice_init(
        FT;
        tspan = tspan,
        dt = Δt_cpl,
        space = boundary_space,
        saveat = saveat,
        area_fraction = ice_fraction,
        thermo_params = thermo_params,
    )

    ## CO2 concentration
    CO2_info = bcfile_info_init(
        FT,
        REGRID_DIR,
        co2_data,
        "co2",
        boundary_space,
        comms_ctx,
        interpolate_daily = true,
        land_fraction = ones(boundary_space),
        date0 = date0,
        mono = mono_surface,
    )

    update_midmonth_data!(date0, CO2_info)
    CO2_init = interpolate_midmonth_to_daily(date0, CO2_info)
    update_field!(atmos_sim, Val(:co2_gm), CO2_init)

    mode_specifics = (; name = mode_name, SST_info = SST_info, SIC_info = SIC_info, CO2_info = CO2_info)

elseif mode_name in ("slabplanet", "slabplanet_aqua", "slabplanet_terra")

    land_fraction = mode_name == "slabplanet_aqua" ? land_fraction .* 0 : land_fraction
    land_fraction = mode_name == "slabplanet_terra" ? land_fraction .* 0 .+ 1 : land_fraction

    ## land
    land_sim = bucket_init(
        FT,
        tspan,
        config_dict["land_domain_type"],
        config_dict["land_albedo_type"],
        config_dict["land_temperature_anomaly"],
        comms_ctx,
        REGRID_DIR;
        dt = Δt_cpl,
        space = boundary_space,
        saveat = saveat,
        area_fraction = land_fraction,
        date_ref = date0,
        t_start = t_start,
    )

    ## ocean
    ocean_sim = ocean_init(
        FT;
        tspan = tspan,
        dt = Δt_cpl,
        space = boundary_space,
        saveat = saveat,
        area_fraction = (FT(1) .- land_fraction), ## NB: this ocean fraction includes areas covered by sea ice (unlike the one contained in the cs)
        thermo_params = thermo_params,
        evolving = evolving_ocean,
    )

    ## sea ice (here set to zero area coverage)
    ice_sim = SurfaceStub((;
        T_sfc = ClimaCore.Fields.ones(boundary_space),
        ρ_sfc = ClimaCore.Fields.zeros(boundary_space),
        z0m = FT(0),
        z0b = FT(0),
        beta = FT(1),
        α = FT(1),
        area_fraction = ClimaCore.Fields.zeros(boundary_space),
        phase = TD.Ice(),
        thermo_params = thermo_params,
    ))

    mode_specifics = (; name = mode_name, SST_info = nothing, SIC_info = nothing)
end

#=
## Coupler Initialization
The coupler needs to contain exchange information, manage the calendar and be able to access all component models. It can also optionally
save online diagnostics. These are all initialized here and saved in a global `CouplerSimulation` struct, `cs`.
=#

## coupler exchange fields
coupler_field_names = (
    :T_S,
    :z0m_S,
    :z0b_S,
    :ρ_sfc,
    :q_sfc,
    :albedo,
    :beta,
    :F_turb_energy,
    :F_turb_moisture,
    :F_turb_ρτxz,
    :F_turb_ρτyz,
    :F_radiative,
    :P_liq,
    :P_snow,
    :F_radiative_TOA,
    :P_net,
)
coupler_fields =
    NamedTuple{coupler_field_names}(ntuple(i -> ClimaCore.Fields.zeros(boundary_space), length(coupler_field_names)))

## model simulations
model_sims = (atmos_sim = atmos_sim, ice_sim = ice_sim, land_sim = land_sim, ocean_sim = ocean_sim);

## dates
dates = (; date = [date], date0 = [date0], date1 = [Dates.firstdayofmonth(date0)], new_month = [false])

#=
### Online Diagnostics
User can write custom diagnostics in the `user_diagnostics.jl`.
=#
monthly_3d_diags = init_diagnostics(
    (:T, :u, :q_tot, :q_liq_ice),
    atmos_sim.domain.center_space;
    save = Monthly(),
    operations = (; accumulate = TimeMean([Int(0)])),
    output_dir = COUPLER_OUTPUT_DIR,
    name_tag = "monthly_mean_3d_",
)

monthly_2d_diags = init_diagnostics(
    (:precipitation_rate, :toa_fluxes, :T_sfc, :tubulent_energy_fluxes),
    boundary_space;
    save = Monthly(),
    operations = (; accumulate = TimeMean([Int(0)])),
    output_dir = COUPLER_OUTPUT_DIR,
    name_tag = "monthly_mean_2d_",
)

diagnostics = (monthly_3d_diags, monthly_2d_diags)

#=
## Initialize Conservation Checks
=#
## init conservation info collector
conservation_checks = nothing
if energy_check
    @assert(
        mode_name[1:10] == "slabplanet" && !CA.is_distributed(ClimaComms.context(boundary_space)),
        "Only non-distributed slabplanet allowable for energy_check"
    )
    conservation_checks = (; energy = EnergyConservationCheck(model_sims), water = WaterConservationCheck(model_sims))
end

dir_paths = (; output = COUPLER_OUTPUT_DIR, artifacts = COUPLER_ARTIFACTS_DIR)
checkpoint_cb =
    HourlyCallback(dt = FT(480), func = checkpoint_sims, ref_date = [dates.date[1]], active = hourly_checkpoint) # 20 days
update_firstdayofmonth!_cb =
    MonthlyCallback(dt = FT(1), func = update_firstdayofmonth!, ref_date = [dates.date1[1]], active = true) # for BCReader
callbacks = (; checkpoint = checkpoint_cb, update_firstdayofmonth! = update_firstdayofmonth!_cb)
