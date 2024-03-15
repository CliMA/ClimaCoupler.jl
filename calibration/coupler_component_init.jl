#=
## Data File Paths
The data files are downloaded from the `ClimaCoupler` artifacts directory. If the data files are not present, they are downloaded from the
original sources.
=#

include(joinpath(pkgdir(ClimaCoupler), "artifacts", "artifact_funcs.jl"))
sst_data = joinpath(sst_dataset_path(), "sst.nc")
sic_data = joinpath(sic_dataset_path(), "sic.nc")
co2_data = joinpath(co2_dataset_path(), "mauna_loa_co2.nc")
land_mask_data = joinpath(mask_dataset_path(), "seamask.nc")

#=
## Component Model Initialization
Here we set initial and boundary conditions for each component model. Each component model is required to have an `init` function that
returns a `ComponentModelSimulation` object (see `Interfacer` docs for more details).
=#

#=
### Atmosphere
This uses the `ClimaAtmos.jl` model, with parameterization options specified in the `config_dict_atmos` dictionary.
=#

## init atmos model component
atmos_sim = atmos_init(FT, config_dict_atmos);
thermo_params = get_thermo_params(atmos_sim) # TODO: this should be shared by all models #610

#=
### Boundary Space
We use a common `Space` for all global surfaces. This enables the MPI processes to operate on the same columns in both
the atmospheric and surface components, so exchanges are parallelized. Note this is only possible when the
atmosphere and surface are of the same horizontal resolution.
=#

## init a 2D boundary space at the surface
boundary_space = ClimaCore.Spaces.horizontal_space(atmos_sim.domain.face_space) # TODO: specify this in the coupler and pass it to all component models #665

#=
### Land-sea Fraction
This is a static field that contains the area fraction of land and sea, ranging from 0 to 1. If applicable, sea ice is included in the sea fraction. at this stage.
Note that land-sea area fraction is different to the land-sea mask, which is a binary field (masks are used internally by the coupler to indicate passive cells that are not populated by a given component model).
=#

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

#=
### Surface Models: AMIP and SlabPlanet Modes
Both modes evolve `ClimaLand.jl`'s bucket model.

In the `AMIP` mode, all ocean properties are prescribed from a file, while sea-ice temperatures are calculated using observed
SIC and assuming a 2m thickness of the ice.

In the `SlabPlanet` mode, all ocean and sea ice are dynamical models, namely thermal slabs, with different parameters. We have several `SlabPlanet` versions
- `slabplanet` = land + slab ocean
- `slabplanet_aqua` = slab ocean everywhere
- `slabplanet_terra` = land everywhere
- `slabplanet_eisenman` = land + slab ocean + slab sea ice with an evolving thickness
=#

ClimaComms.iamroot(comms_ctx) ? @info(mode_name) : nothing
if mode_name == "amip"
    ClimaComms.iamroot(comms_ctx) ? @info("AMIP boundary conditions - do not expect energy conservation") : nothing

    ## land model
    land_sim = bucket_init(
        FT,
        tspan,
        config_dict["land_domain_type"],
        config_dict["land_albedo_type"],
        config_dict["land_temperature_anomaly"],
        REGRID_DIR;
        dt = Δt_cpl,
        space = boundary_space,
        saveat = saveat,
        area_fraction = land_fraction,
        date_ref = date0,
        t_start = t_start,
    )

    ## ocean stub
    SST_info = bcfile_info_init(
        FT,
        REGRID_DIR,
        sst_data,
        "SST",
        boundary_space,
        comms_ctx,
        interpolate_daily = true,
        scaling_function = scale_sst, ## convert to Kelvin
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

    ## sea ice model
    SIC_info = bcfile_info_init(
        FT,
        REGRID_DIR,
        sic_data,
        "SEAICE",
        boundary_space,
        comms_ctx,
        interpolate_daily = true,
        scaling_function = scale_sic, ## convert to fraction
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

    ## CO2 concentration from temporally varying file
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
    update_field!(atmos_sim, Val(:co2), CO2_init)

    mode_specifics = (; name = mode_name, SST_info = SST_info, SIC_info = SIC_info, CO2_info = CO2_info)

elseif mode_name in ("slabplanet", "slabplanet_aqua", "slabplanet_terra")


    land_fraction = mode_name == "slabplanet_aqua" ? land_fraction .* 0 : land_fraction
    land_fraction = mode_name == "slabplanet_terra" ? land_fraction .* 0 .+ 1 : land_fraction

    ## land model
    land_sim = bucket_init(
        FT,
        tspan,
        config_dict["land_domain_type"],
        config_dict["land_albedo_type"],
        config_dict["land_temperature_anomaly"],
        REGRID_DIR;
        dt = Δt_cpl,
        space = boundary_space,
        saveat = saveat,
        area_fraction = land_fraction,
        date_ref = date0,
        t_start = t_start,
    )

    ## ocean model
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

    ## sea ice stub (here set to zero area coverage)
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

elseif mode_name == "slabplanet_eisenman"

    ## land model
    land_sim = bucket_init(
        FT,
        tspan,
        config_dict["land_domain_type"],
        config_dict["land_albedo_type"],
        config_dict["land_temperature_anomaly"],
        REGRID_DIR;
        dt = Δt_cpl,
        space = boundary_space,
        saveat = saveat,
        area_fraction = land_fraction,
        date_ref = date0,
        t_start = t_start,
    )

    ## ocean stub (here set to zero area coverage)
    ocean_sim = ocean_init(
        FT;
        tspan = tspan,
        dt = Δt_cpl,
        space = boundary_space,
        saveat = saveat,
        area_fraction = ClimaCore.Fields.zeros(boundary_space), # zero, since ML is calculated below
        thermo_params = thermo_params,
    )

    ## sea ice + ocean model
    ice_sim = eisenman_seaice_init(
        FT,
        tspan,
        space = boundary_space,
        area_fraction = (FT(1) .- land_fraction),
        dt = Δt_cpl,
        saveat = saveat,
        thermo_params = thermo_params,
    )

    mode_specifics = (; name = mode_name, SST_info = nothing, SIC_info = nothing)
end