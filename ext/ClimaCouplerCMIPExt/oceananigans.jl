import ClimaComms
import SurfaceFluxes as SF
import Thermodynamics as TD
import ClimaOcean.EN4: download_dataset

include("climaocean_helpers.jl")
include("ocean_skin.jl")
import Dates

"""
    OceananigansSimulation{SIM, A, OPROP, REMAP, SIC}

The ClimaCoupler simulation object used to run with Oceananigans.
This type is used by the coupler to indicate that this simulation
is a surface/ocean simulation for dispatch.

It contains the following objects:
- `ocean::SIM`: The Oceananigans simulation object.
- `area_fraction::A`: A ClimaCore Field representing the surface area fraction of this component model on the exchange grid.
- `ocean_properties::OPROP`: A NamedTuple of ocean properties and parameters (including COARE3 roughness params).
- `remapping::REMAP`: Objects needed to remap from the exchange (spectral) grid to Oceananigans spaces.
- `ice_concentration::SIC`: An Oceananigans Field representing the sea ice concentration on the ocean/sea ice grid.
"""
struct OceananigansSimulation{SIM, A, OPROP, REMAP, SIC} <:
       Interfacer.AbstractOceanSimulation
    ocean::SIM
    area_fraction::A
    ocean_properties::OPROP
    remapping::REMAP
    ice_concentration::SIC
end

"""
    Interfacer.OceanSimulation(::Type{FT}, ::Val{:oceananigans}; kwargs...)

Extension of the generic OceanSimulation constructor for the Oceananigans ocean model.
FT is accepted for consistency with other ocean models but is not used.
"""
function Interfacer.OceanSimulation(::Type{FT}, ::Val{:oceananigans}; kwargs...) where {FT}
    return OceananigansSimulation(FT; kwargs...)
end

"""
    OceananigansSimulation(; kwargs...)

Creates an OceananigansSimulation object containing a model, an integrator, and
a surface area fraction field.
This type is used to indicate that this simulation is an ocean simulation for
dispatch in coupling.

# Required keyword arguments
- `boundary_space`: The boundary space of the coupled simulation
- `start_date`: Start date for the simulation
- `stop_date`: Stop date for the simulation
- `output_dir`: Directory for output files
- `ice_model`: Ice model type (Val type)

# Optional keyword arguments
- `dt`: Time step (default: `nothing`)
- `comms_ctx`: Communication context (default: `ClimaComms.context()`)
- `coupled_param_dict`: Coupled parameter dictionary (default: created from `area_fraction`)

Specific details about the default model configuration
can be found in the documentation for `ClimaOcean.ocean_simulation`.
"""
function OceananigansSimulation(
    ::Type{FT};
    boundary_space,
    start_date,
    tspan,
    output_dir,
    ice_model,
    dt = nothing,
    comms_ctx = ClimaComms.context(),
    coupled_param_dict = CP.create_toml_dict(FT),
    use_iterative_ocean_skin = false,
    skin_layer_thickness = FT(1e-3),
    thermal_diffusivity = FT(1.4e-7),
    extra_kwargs...,
) where {FT}
    arch = comms_ctx.device isa ClimaComms.CUDADevice ? OC.GPU() : OC.CPU()
    OC.Oceananigans.defaults.FloatType = FT

    # Compute stop_date for oceananigans (needed for EN4 data retrieval)
    stop_date = start_date + Dates.Second(float(tspan[2] - tspan[1]))

    # Use Float64 for the ocean to avoid precision issues
    FT_ocean = Float64
    OC.Oceananigans.defaults.FloatType = FT_ocean

    # Retrieve EN4 data (monthly)
    # (It requires username and password)
    dates = range(start_date, step = Dates.Month(1), stop = stop_date)
    en4_temperature = CO.Metadata(:temperature; dates, dataset = CO.EN4.EN4Monthly())
    en4_salinity = CO.Metadata(:salinity; dates, dataset = CO.EN4.EN4Monthly())
    download_dataset(en4_temperature)
    download_dataset(en4_salinity)

    # Set up ocean grid (1 degree)
    resolution_points = (360, 160, 32)
    Nz = last(resolution_points)
    depth = 4000 # meters
    z = OC.ExponentialDiscretization(Nz, -depth, 0; scale = 0.85 * depth)

    # Regular LatLong because we know how to do interpolation there
    underlying_grid = OC.LatitudeLongitudeGrid(
        arch;
        size = resolution_points,
        longitude = (-180, 180),
        latitude = (-80, 80),   # NOTE: Don't goo to high up when using LatLongGrid, or the cells will be too small
        z,
        halo = (7, 7, 7),
    )

    bottom_height = CO.regrid_bathymetry(
        underlying_grid;
        minimum_depth = 30,
        interpolation_passes = 20,
        major_basins = 1,
    )

    grid = OC.ImmersedBoundaryGrid(
        underlying_grid,
        OC.GridFittedBottom(bottom_height);
        active_cells_map = true,
    )

    # Restore the ocean to the EN4 state periodically if running for more than one month and not using ClimaSeaIce
    use_restoring =
        start_date + Dates.Month(1) < stop_date && ice_model != Val(:clima_seaice)

    if use_restoring
        # When we use EN4 data, the forcing takes care of everything, including
        # the initial conditions
        restoring_rate = 1 / (3 * 86400)
        mask = CO.LinearlyTaperedPolarMask(
            southern = (-80, -70),
            northern = (70, 90),
            z = (z(1), 0),
        )

        forcing_T = CO.DatasetRestoring(en4_temperature, grid; mask, rate = restoring_rate)
        forcing_S = CO.DatasetRestoring(en4_salinity, grid; mask, rate = restoring_rate)
        forcing = (T = forcing_T, S = forcing_S)
    else
        forcing = (;)
    end

    # Create ocean simulation
    free_surface = OC.SplitExplicitFreeSurface(grid; substeps = 70)
    momentum_advection = OC.VectorInvariant()
    tracer_advection = OC.WENO(order = 5)
    horizontal_viscosity = OC.HorizontalScalarBiharmonicDiffusivity(ν = 1e11)

    # Use Float32 for the vertical mixing parameters to avoid parameter memory limits
    vertical_mixing = OC.CATKEVerticalDiffusivity(Float32)

    Δt = isnothing(dt) ? CO.OceanSimulations.estimate_maximum_Δt(grid) : dt
    ocean = CO.ocean_simulation(
        grid;
        Δt,
        forcing,
        momentum_advection,
        tracer_advection,
        free_surface,
        closure = (horizontal_viscosity, vertical_mixing),
    )

    # Set initial condition to EN4 state estimate at start_date
    OC.set!(ocean.model, T = en4_temperature[1], S = en4_salinity[1])

    long_cc = OC.λnodes(grid, OC.Center(), OC.Center(), OC.Center())
    lat_cc = OC.φnodes(grid, OC.Center(), OC.Center(), OC.Center())

    # TODO: Go from 0 to Nx+1, Ny+1 (for halos) (for LatLongGrid)

    # Construct a remapper from the exchange grid to `Center, Center` fields
    long_cc = reshape(long_cc, length(long_cc), 1)
    lat_cc = reshape(lat_cc, 1, length(lat_cc))
    target_points_cc = @. CC.Geometry.LatLongPoint(lat_cc, long_cc)
    # TODO: We can remove the `nothing` after CC > 0.14.33
    remapper_cc = CC.Remapping.Remapper(boundary_space, target_points_cc, nothing)

    # Construct two 2D Center/Center fields to use as scratch space while remapping
    scratch_cc1 = OC.Field{OC.Center, OC.Center, Nothing}(grid)
    scratch_cc2 = OC.Field{OC.Center, OC.Center, Nothing}(grid)

    # Construct two scratch arrays to use while remapping
    # We get the array type and dimensions from the remapper object to maintain consistency
    ArrayType = ClimaComms.array_type(remapper_cc.space)
    interpolated_values_dim..., _buffer_length = size(remapper_cc._interpolated_values)

    scratch_arr1 = ArrayType(zeros(FT, interpolated_values_dim...))
    scratch_arr2 = ArrayType(zeros(FT, interpolated_values_dim...))
    scratch_arr3 = ArrayType(zeros(FT, interpolated_values_dim...))

    # Allocate space for a Field of UVVectors, which we need for remapping momentum fluxes
    temp_uv_vec = CC.Fields.Field(CC.Geometry.UVVector{FT}, boundary_space)

    remapping = (;
        remapper_cc,
        scratch_cc1,
        scratch_cc2,
        scratch_arr1,
        scratch_arr2,
        scratch_arr3,
        temp_uv_vec,
    )

    # COARE3 roughness params (allocated once, reused each timestep)
    coare3_roughness_params = CC.Fields.Field(SF.COARE3RoughnessParams{FT}, boundary_space)
    coare3_roughness_params .= SF.COARE3RoughnessParams{FT}()

    # Get some ocean properties and parameters (including COARE3 roughness params)
    ocean_properties = (;
        reference_density = 1020,
        heat_capacity = 3991,
        σ = coupled_param_dict["stefan_boltzmann_constant"],
        C_to_K = coupled_param_dict["temperature_water_freeze"],
        use_iterative_ocean_skin,
        skin_layer_thickness,
        thermal_diffusivity,
        coare3_roughness_params,
    )

    # Before version 0.96.22, the NetCDFWriter was broken on GPU
    if arch isa OC.CPU
        # Save all tracers and velocities to a NetCDF file at daily frequency
        outputs = merge(ocean.model.tracers, ocean.model.velocities)
        surface_writer = OC.NetCDFWriter(
            ocean.model,
            outputs;
            schedule = OC.TimeInterval(86400), # Daily output
            filename = joinpath(output_dir, "ocean_diagnostics.nc"),
            indices = (:, :, grid.Nz),
            overwrite_existing = true,
            array_type = Array{FT},
        )
        free_surface_writer = OC.NetCDFWriter(
            ocean.model,
            (; η = ocean.model.free_surface.displacement);
            schedule = OC.TimeInterval(3600), # hourly snapshots
            filename = joinpath(output_dir, "ocean_free_surface.nc"),
            overwrite_existing = true,
            array_type = Array{FT},
        )
        Tflux = ocean.model.tracers.T.boundary_conditions.top.condition
        Sflux = ocean.model.tracers.S.boundary_conditions.top.condition
        uflux = ocean.model.velocities.u.boundary_conditions.top.condition
        vflux = ocean.model.velocities.v.boundary_conditions.top.condition
        fluxes_writer = OC.NetCDFWriter(
            ocean.model,
            (; Tflux, Sflux, uflux, vflux);
            schedule = OC.TimeInterval(3600), # hourly snapshots
            filename = joinpath(output_dir, "ocean_fluxes.nc"),
            overwrite_existing = true,
            array_type = Array{FT},
        )

        ocean.output_writers[:surface] = surface_writer
        ocean.output_writers[:free_surface] = free_surface_writer
        ocean.output_writers[:fluxes] = fluxes_writer
    end

    # Initialize with 0 ice concentration; this will be updated in `resolve_area_fractions!`
    # if the ocean is coupled to a non-prescribed sea ice model.
    ice_concentration = OC.Field{OC.Center, OC.Center, Nothing}(grid)

    # Create a dummy area fraction that will get overwritten in `update_surface_fractions!`
    area_fraction = ones(boundary_space)

    return OceananigansSimulation(
        ocean,
        area_fraction,
        ocean_properties,
        remapping,
        ice_concentration,
    )
end

"""
    FieldExchanger.resolve_area_fractions!(ocean_sim, ice_sim, land_fraction)

Ensure the ocean and ice area fractions are consistent with each other.
This matters in the case of a LatitudeLongitudeGrid, which is only
defined between -80 and 80 degrees latitude. In this case, we set the ice
and ocean area fractions to 0 and the land fraction to 1 on [78°S, 90°S]
and on [78°N, 90°N]. Note the overlap between 78° and 80° to avoid any gaps
introduced by the cubed sphere not aligning with latitude lines.

This function also updates the ice concentration field in the ocean simulation
so that it can be used for weighting flux updates.
"""
function FieldExchanger.resolve_area_fractions!(
    ocean_sim::OceananigansSimulation,
    ice_sim,
    land_fraction,
)
    if ocean_sim.ocean.model.grid.underlying_grid isa OC.LatitudeLongitudeGrid
        ocean_fraction = Interfacer.get_field(ocean_sim, Val(:area_fraction))
        ice_fraction = Interfacer.get_field(ice_sim, Val(:area_fraction))

        # Create a "polar" mask that's 1 at latitudes in [-90, -80] and [80, 90] degrees
        boundary_space = axes(ocean_fraction)
        FT = CC.Spaces.undertype(boundary_space)
        lat = CC.Fields.coordinate_field(boundary_space).lat
        polar_mask = CC.Fields.zeros(boundary_space)
        polar_mask .= abs.(lat) .>= FT(78)

        # Set land fraction to 1 and ice/ocean fraction to 0 where polar_mask is 1
        @. land_fraction = ifelse.(polar_mask == FT(1), FT(1), land_fraction)
        @. ice_fraction = ifelse.(polar_mask == FT(1), FT(0), ice_fraction)
        @. ocean_fraction = ifelse.(polar_mask == FT(1), FT(0), ocean_fraction)
    end

    # Update the ice concentration field in the ocean simulation
    ice_sim isa ClimaSeaIceSimulation && (
        ocean_sim.ice_concentration .=
            Interfacer.get_field(ice_sim, Val(:ice_concentration))
    )
    return nothing
end

###############################################################################
### Functions required by ClimaCoupler.jl for a AbstractSurfaceSimulation
###############################################################################

# Timestep the simulation forward to time `t`
Interfacer.step!(sim::OceananigansSimulation, t) =
    OC.time_step!(sim.ocean, float(t) - sim.ocean.model.clock.time)

Interfacer.get_field(sim::OceananigansSimulation, ::Val{:area_fraction}) = sim.area_fraction

# TODO: Better values for this
Interfacer.get_field(sim::OceananigansSimulation, ::Val{:roughness_model}) = :coare3
Interfacer.get_field(sim::OceananigansSimulation, ::Val{:coare3_roughness_params}) =
    sim.ocean_properties.coare3_roughness_params
Interfacer.get_field(sim::OceananigansSimulation, ::Val{:roughness_buoyancy}) =
    eltype(sim.ocean.model)(5.8e-5)
Interfacer.get_field(sim::OceananigansSimulation, ::Val{:roughness_momentum}) =
    eltype(sim.ocean.model)(5.8e-5)
Interfacer.get_field(sim::OceananigansSimulation, ::Val{:emissivity}) =
    eltype(sim.ocean.model)(0.97)
Interfacer.get_field(sim::OceananigansSimulation, ::Val{:surface_direct_albedo}) =
    eltype(sim.ocean.model)(0.011)
Interfacer.get_field(sim::OceananigansSimulation, ::Val{:surface_diffuse_albedo}) =
    eltype(sim.ocean.model)(0.069)

# NOTE: This is 3D, but it will be remapped to 2D
Interfacer.get_field(sim::OceananigansSimulation, ::Val{:surface_temperature}) =
    sim.ocean.model.tracers.T + sim.ocean_properties.C_to_K # convert from Celsius to Kelvin

Interfacer.get_field(sim::OceananigansSimulation, ::Val{:use_iterative_ocean_skin}) =
    sim.ocean_properties.use_iterative_ocean_skin
Interfacer.get_field(sim::OceananigansSimulation, ::Val{:skin_layer_thickness}) =
    sim.ocean_properties.skin_layer_thickness
Interfacer.get_field(sim::OceananigansSimulation, ::Val{:thermal_diffusivity}) =
    sim.ocean_properties.thermal_diffusivity

"""
    FluxCalculator.compute_surface_fluxes!(csf, sim::OceananigansSimulation, atmos_sim, thermo_params)

When `use_iterative_ocean_skin` is true, solve for the ocean skin temperature (DiffusiveFlux,
following ClimaOcean's SkinTemperature{<:DiffusiveFlux}) using SurfaceFluxes.jl's
`T_sfc_guess` and `update_T_sfc` callback: a single call to SurfaceFluxes.surface_fluxes
with an `update_T_sfc` callback that returns T_sfc satisfying the ocean flux balance
at each Monin-Obukhov iteration step. Otherwise invokes the base implementation.
"""
function FluxCalculator.compute_surface_fluxes!(
    csf,
    sim::OceananigansSimulation,
    atmos_sim::Interfacer.AbstractAtmosSimulation,
    thermo_params,
)
    if !Interfacer.get_field(sim, Val(:use_iterative_ocean_skin))
        return invoke(
            FluxCalculator.compute_surface_fluxes!,
            (Any, Interfacer.AbstractSurfaceSimulation, Interfacer.AbstractAtmosSimulation, Any),
            csf,
            sim,
            atmos_sim,
            thermo_params,
        )
    end

    boundary_space = axes(csf)
    surface_fluxes_params = FluxCalculator.get_surface_params(atmos_sim)
    uv_int = StaticArrays.SVector.(csf.u_int, csf.v_int)
    thermo_state_atmos =
        TD.PhaseNonEquil_ρTq.(
            thermo_params,
            csf.ρ_atmos,
            csf.T_atmos,
            TD.PhasePartition.(csf.q_atmos),
        )

    # Bulk ocean T (top cell) and skin params on boundary space
    Interfacer.get_field!(csf.scalar_temp1, sim, Val(:surface_temperature))
    T_ocean = csf.scalar_temp1
    δ = Interfacer.get_field(sim, Val(:skin_layer_thickness))
    κ = Interfacer.get_field(sim, Val(:thermal_diffusivity))
    α = Interfacer.get_field(sim, Val(:surface_direct_albedo))
    ε = Interfacer.get_field(sim, Val(:emissivity))
    σ = sim.ocean_properties.σ
    ρ = sim.ocean_properties.reference_density
    c = sim.ocean_properties.heat_capacity

    # Initial surface state for q_vap_sfc_guess (T_sfc_guess = T_ocean passed separately)
    ρ_sfc =
        SF.surface_density.(
            surface_fluxes_params,
            csf.T_atmos,
            csf.ρ_atmos,
            T_ocean,
            csf.height_int .- csf.height_sfc,
            csf.q_atmos,
            0,
            0,
        )
    FluxCalculator.compute_surface_humidity!(csf.scalar_temp3, T_ocean, ρ_sfc, thermo_params)
    q_sfc = csf.scalar_temp3
    thermo_state_sfc =
        TD.PhaseNonEquil_ρTq.(thermo_params, ρ_sfc, T_ocean, TD.PhasePartition.(q_sfc))
    Interfacer.get_field!(csf.scalar_temp2, sim, Val(:roughness_momentum))
    z0m = csf.scalar_temp2
    Interfacer.get_field!(csf.scalar_temp3, sim, Val(:roughness_buoyancy))
    z0b = csf.scalar_temp3
    config =
        SF.SurfaceFluxConfig.(
            SF.ConstantRoughnessParams.(z0m, z0b),
            SF.ConstantGustinessSpec.(ones(boundary_space)),
        )

    # Single SurfaceFluxes call with T_sfc_guess and update_T_sfc: skin T_sfc solved inside MOST iteration
    ocean_update_T_sfc_callbacks =
        ocean_update_T_sfc_callback.(
            csf.SW_d,
            csf.LW_d,
            α,
            ε,
            σ,
            T_ocean,
            δ,
            κ,
            Ref(ρ),
            Ref(c),
        )
    fluxes =
        FluxCalculator.get_surface_fluxes.(
            surface_fluxes_params,
            uv_int,
            thermo_state_atmos,
            csf.height_int,
            StaticArrays.SVector.(0, 0),
            thermo_state_sfc,
            csf.height_sfc,
            0,
            config;
            T_sfc_guess = T_ocean,
            update_T_sfc = ocean_update_T_sfc_callbacks,
        )

    (; F_turb_ρτxz, F_turb_ρτyz, F_sh, F_lh, F_turb_moisture, L_MO, ustar, buoyancy_flux) =
        fluxes
    Interfacer.get_field!(csf.scalar_temp1, sim, Val(:area_fraction))
    area_fraction = csf.scalar_temp1

    @. F_turb_ρτxz = ifelse(area_fraction ≈ 0, zero(F_turb_ρτxz), F_turb_ρτxz)
    @. F_turb_ρτyz = ifelse(area_fraction ≈ 0, zero(F_turb_ρτyz), F_turb_ρτyz)
    @. F_sh = ifelse(area_fraction ≈ 0, zero(F_sh), F_sh)
    @. F_lh = ifelse(area_fraction ≈ 0, zero(F_lh), F_lh)
    @. F_turb_moisture = ifelse(area_fraction ≈ 0, zero(F_turb_moisture), F_turb_moisture)
    @. L_MO = ifelse(area_fraction ≈ 0, zero(L_MO), L_MO)
    @. ustar = ifelse(area_fraction ≈ 0, zero(ustar), ustar)
    @. buoyancy_flux = ifelse(area_fraction ≈ 0, zero(buoyancy_flux), buoyancy_flux)

    fields = (; F_turb_ρτxz, F_turb_ρτyz, F_lh, F_sh, F_turb_moisture)
    FluxCalculator.update_turbulent_fluxes!(sim, fields)

    @. csf.F_turb_ρτxz += F_turb_ρτxz * area_fraction
    @. csf.F_turb_ρτyz += F_turb_ρτyz * area_fraction
    @. csf.F_lh += F_lh * area_fraction
    @. csf.F_sh += F_sh * area_fraction
    @. csf.F_turb_moisture += F_turb_moisture * area_fraction
    @. csf.L_MO += ifelse(isinf(L_MO), L_MO, L_MO * area_fraction)
    @. csf.ustar += ustar * area_fraction
    @. csf.buoyancy_flux += buoyancy_flux * area_fraction
    return nothing
end

"""
    FluxCalculator.update_turbulent_fluxes!(sim::OceananigansSimulation, fields)

Update the turbulent fluxes in the simulation using the values computed at this time step.
These include latent heat flux, sensible heat flux, momentum fluxes, and moisture flux.

Rather than setting the surface fluxes and overwriting previous values, this function adds only
the contributions from the turbulent fluxes. `update_sim!` sets the surface fluxes due to
radiation and precipitation. Additional contributions may be made in `ocean_seaice_fluxes!`.
An exception is the momentum fluxes, which are set directly here since they are not updated
in `update_sim!`.

A note on sign conventions:
SurfaceFluxes and Oceananigans both use the convention that a positive flux is an upward flux.
No sign change is needed during the exchange, except for moisture/salinity fluxes:
SurfaceFluxes provides moisture moving from atmosphere to ocean as a negative flux at the surface,
and Oceananigans represents moisture moving from atmosphere to ocean as a positive salinity flux,
so a sign change is needed when we convert from moisture to salinity flux.
"""
function FluxCalculator.update_turbulent_fluxes!(sim::OceananigansSimulation, fields)
    (; F_lh, F_sh, F_turb_ρτxz, F_turb_ρτyz, F_turb_moisture) = fields
    (; reference_density, heat_capacity) = sim.ocean_properties
    grid = sim.ocean.model.grid
    ice_concentration = sim.ice_concentration

    # Convert the momentum fluxes from contravariant to Cartesian basis
    contravariant_to_cartesian!(sim.remapping.temp_uv_vec, F_turb_ρτxz, F_turb_ρτyz)
    F_turb_ρτxz_uv = sim.remapping.temp_uv_vec.components.data.:1
    F_turb_ρτyz_uv = sim.remapping.temp_uv_vec.components.data.:2

    # Remap momentum fluxes onto reduced 2D Center, Center fields using scratch arrays and fields
    CC.Remapping.interpolate!(
        sim.remapping.scratch_arr1,
        sim.remapping.remapper_cc,
        F_turb_ρτxz_uv,
    )
    OC.set!(sim.remapping.scratch_cc1, sim.remapping.scratch_arr1) # zonal momentum flux
    CC.Remapping.interpolate!(
        sim.remapping.scratch_arr2,
        sim.remapping.remapper_cc,
        F_turb_ρτyz_uv,
    )
    OC.set!(sim.remapping.scratch_cc2, sim.remapping.scratch_arr2) # meridional momentum flux

    # Rename for clarity; these are now Center, Center Oceananigans fields
    F_turb_ρτxz_cc = sim.remapping.scratch_cc1
    F_turb_ρτyz_cc = sim.remapping.scratch_cc2

    # Weight by (1 - sea ice concentration)
    OC.interior(F_turb_ρτxz_cc, :, :, 1) .=
        OC.interior(F_turb_ρτxz_cc, :, :, 1) .* (1.0 .- ice_concentration) ./
        reference_density
    OC.interior(F_turb_ρτyz_cc, :, :, 1) .=
        OC.interior(F_turb_ρτyz_cc, :, :, 1) .* (1.0 .- ice_concentration) ./
        reference_density

    # Set the momentum flux BCs at the correct locations using the remapped scratch fields
    oc_flux_u = surface_flux(sim.ocean.model.velocities.u)
    oc_flux_v = surface_flux(sim.ocean.model.velocities.v)
    set_from_extrinsic_vector!(
        (; u = oc_flux_u, v = oc_flux_v),
        grid,
        F_turb_ρτxz_cc,
        F_turb_ρτyz_cc,
    )

    # Remap the latent and sensible heat fluxes using scratch arrays
    CC.Remapping.interpolate!(sim.remapping.scratch_arr1, sim.remapping.remapper_cc, F_lh) # latent heat flux
    CC.Remapping.interpolate!(sim.remapping.scratch_arr2, sim.remapping.remapper_cc, F_sh) # sensible heat flux

    # Rename for clarity; recall F_turb_energy = F_lh + F_sh
    remapped_F_lh = sim.remapping.scratch_arr1
    remapped_F_sh = sim.remapping.scratch_arr2

    # TODO: Note, SW radiation penetrates the surface. Right now, we just put
    # everything on the surface, but later we will need to account for this.
    # One way we can do this is using directly ClimaOcean
    oc_flux_T = surface_flux(sim.ocean.model.tracers.T)
    OC.interior(oc_flux_T, :, :, 1) .=
        OC.interior(oc_flux_T, :, :, 1) .+
        (1.0 .- ice_concentration) .* (remapped_F_lh .+ remapped_F_sh) ./
        (reference_density * heat_capacity)

    # Add the part of the salinity flux that comes from the moisture flux, we also need to
    # add the component due to precipitation (that was done with the radiative fluxes)
    CC.Remapping.interpolate!(
        sim.remapping.scratch_arr1,
        sim.remapping.remapper_cc,
        F_turb_moisture,
    )
    moisture_fresh_water_flux = sim.remapping.scratch_arr1 ./ reference_density
    oc_flux_S = surface_flux(sim.ocean.model.tracers.S)
    surface_salinity = OC.interior(sim.ocean.model.tracers.S, :, :, grid.Nz)
    OC.interior(oc_flux_S, :, :, 1) .=
        OC.interior(oc_flux_S, :, :, 1) .-
        (1.0 .- ice_concentration) .* surface_salinity .* moisture_fresh_water_flux
    return nothing
end

function Interfacer.update_field!(sim::OceananigansSimulation, ::Val{:area_fraction}, field)
    sim.area_fraction .= field
    return nothing
end

"""
    FieldExchanger.update_sim!(sim::OceananigansSimulation, csf)

Update the ocean simulation with the provided fields, which have been filled in
by the coupler.

Update the portion of the surface_fluxes for T and S that is due to radiation and
precipitation. The rest will be updated in `update_turbulent_fluxes!`.

This function sets the surface fluxes directly, overwriting any previous values.
Additional contributions will be made in `update_turbulent_fluxes!` and `ocean_seaice_fluxes!`.

A note on sign conventions:
ClimaAtmos and Oceananigans both use the convention that a positive flux is an upward flux.
No sign change is needed during the exchange, except for precipitation/salinity fluxes.
ClimaAtmos provides precipitation as a negative flux at the surface, and
Oceananigans represents precipitation as a positive salinity flux,
so a sign change is needed when we convert from precipitation to salinity flux.
"""
function FieldExchanger.update_sim!(sim::OceananigansSimulation, csf)
    (; reference_density, heat_capacity) = sim.ocean_properties
    ice_concentration = sim.ice_concentration
    Nz = sim.ocean.model.grid.Nz

    # Remap radiative flux onto scratch array; rename for clarity
    CC.Remapping.interpolate!(
        sim.remapping.scratch_arr1,
        sim.remapping.remapper_cc,
        csf.SW_d,
    )
    remapped_SW_d = sim.remapping.scratch_arr1

    CC.Remapping.interpolate!(
        sim.remapping.scratch_arr2,
        sim.remapping.remapper_cc,
        csf.LW_d,
    )
    remapped_LW_d = sim.remapping.scratch_arr2

    # Update only the part due to radiative fluxes. For the full update, the component due
    # to latent and sensible heat is missing and will be updated in update_turbulent_fluxes.
    oc_flux_T = surface_flux(sim.ocean.model.tracers.T)
    (; σ, C_to_K) = sim.ocean_properties
    α = Interfacer.get_field(sim, Val(:surface_direct_albedo)) # scalar
    ϵ = Interfacer.get_field(sim, Val(:emissivity)) # scalar
    OC.interior(oc_flux_T, :, :, 1) .=
        (1.0 .- ice_concentration) .* (
            -(1 - α) .* remapped_SW_d .-
            ϵ * (
                remapped_LW_d .-
                σ .* (C_to_K .+ OC.interior(sim.ocean.model.tracers.T, :, :, Nz)) .^ 4
            )
        ) ./ (reference_density * heat_capacity)

    # Remap precipitation fields onto scratch arrays; rename for clarity
    CC.Remapping.interpolate!(
        sim.remapping.scratch_arr1,
        sim.remapping.remapper_cc,
        csf.P_liq,
    )
    CC.Remapping.interpolate!(
        sim.remapping.scratch_arr2,
        sim.remapping.remapper_cc,
        csf.P_snow,
    )
    remapped_P_liq = sim.remapping.scratch_arr1
    remapped_P_snow = sim.remapping.scratch_arr2

    # Virtual salt flux
    oc_flux_S = surface_flux(sim.ocean.model.tracers.S)
    OC.interior(oc_flux_S, :, :, 1) .=
        OC.interior(oc_flux_S, :, :, 1) .-
        OC.interior(sim.ocean.model.tracers.S, :, :, Nz) .* (1.0 .- ice_concentration) .*
        (remapped_P_liq .+ remapped_P_snow) ./ reference_density
    return nothing
end

"""
    get_model_prog_state(sim::OceananigansSimulation)

Returns the model state of a simulation as a `ClimaCore.FieldVector`.
It's okay to leave this unimplemented for now, but we won't be able to use the
restart system.

TODO extend this for non-ClimaCore states.
"""
function Checkpointer.get_model_prog_state(sim::OceananigansSimulation)
    @warn "get_model_prog_state not implemented for OceananigansSimulation"
end
