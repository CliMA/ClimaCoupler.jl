import ClimaComms
import SurfaceFluxes as SF
import Thermodynamics as TD
import ClimaUtilities.TimeManager: ITime, date, counter, period
import Dates

# TODO move to Oceananigans
# Change the `set_to_array!` of Oceananigans to work on a `Center`
# field (in y) with abstract 2D and 3D arrays.
const TG = OC.OrthogonalSphericalShellGrids.TripolarGrid{
    FT,
    TX,
    <:OC.Grids.RightFaceFolded,
} where {FT, TX}
const TGRF = Union{<:OC.ImmersedBoundaryGrid{<:Any, <:Any, <:Any, <:Any, <:TG}, TG}

@inline resize_to_facefolded(a::AbstractArray{T, 2}) where {T} = a[:, 1:(end - 1)]
@inline resize_to_facefolded(a::AbstractArray{T, 3}) where {T} = a[:, 1:(end - 1), :]

function OC.Fields.set_to_array!(
    u::OC.Field{LX, <:OC.Grids.Center, LZ, O, <:TGRF},
    a,
) where {LX, LZ, O}
    a = OC.Architectures.on_architecture(OC.Architectures.CPU(), a)
    a = resize_to_facefolded(a)
    a = OC.Architectures.on_architecture(OC.Architectures.architecture(u), a)

    try
        copyto!(OC.interior(u), a)
    catch err
        if err isa DimensionMismatch
            Nx, Ny, Nz = size(u)
            u .= reshape(a, Nx, Ny, Nz)

            msg = string("Reshaped ", summary(a), " to set! its data to ", '\n', summary(u))
            @warn msg
        else
            throw(err)
        end
    end

    return u
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
- `simple_ocean`: Whether to use a simple ocean model setup

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
    simple_ocean = false,
    dt = 1800.0, # 30 minutes
    comms_ctx = ClimaComms.context(),
    coupled_param_dict = CP.create_toml_dict(FT),
    extra_kwargs...,
) where {FT}
    arch = comms_ctx.device isa ClimaComms.CUDADevice ? OC.GPU() : OC.CPU()
    OC.defaults.FloatType = FT

    # Compute stop_date for oceananigans (needed for EN4 data retrieval)
    stop_date = start_date + Dates.Second(float(tspan[2] - tspan[1]))

    # Retrieve EN4 data (monthly)
    # (It requires username and password)
    dates = range(start_date, step = Dates.Month(1), stop = stop_date)
    en4_temperature = CO.Metadata(:temperature; dates, dataset = CO.EN4.EN4Monthly())
    en4_salinity = CO.Metadata(:salinity; dates, dataset = CO.EN4.EN4Monthly())
    CO.EN4.download_dataset(en4_temperature)
    CO.EN4.download_dataset(en4_salinity)

    # Set up tripolar ocean grid (0.5 degree)
    Nx = 720
    Ny = 360
    Nz = 100
    depth = 6000 # meters
    z = OC.ExponentialDiscretization(Nz, -depth, 0; scale = 1800)

    underlying_grid = OC.TripolarGrid(
        arch;
        size = (Nx, Ny, Nz),
        halo = (7, 7, 7),
        z,
        fold_topology = OC.Grids.RightFaceFolded,
    )

    # TODO revert this
    bottom_height = CO.regrid_bathymetry(
        underlying_grid;
        minimum_depth = 20,
        interpolation_passes = 2,
        major_basins = 1,
    )

    grid = OC.ImmersedBoundaryGrid(
        underlying_grid,
        OC.GridFittedBottom(bottom_height);
        active_cells_map = true,
    )

    # Create ocean simulation
    if !simple_ocean
        free_surface = OC.SplitExplicitFreeSurface(grid; substeps = 150)
        momentum_advection = OC.WENOVectorInvariant(order = 5)
        tracer_advection = OC.WENO(order = 7)
        eddy_closure = OC.TurbulenceClosures.IsopycnalSkewSymmetricDiffusivity(
            κ_skew = 500,
            κ_symmetric = 100,
        )
        @inline νhb(i, j, k, grid, ℓx, ℓy, ℓz, clock, fields, λ) =
            OC.Operators.Az(i, j, k, grid, ℓx, ℓy, ℓz)^2 / λ

        horizontal_viscosity = OC.HorizontalScalarBiharmonicDiffusivity(
            ν = νhb,
            discrete_form = true,
            parameters = 40 * 24 * 60 * 60, # 40 days
        )
        vertical_closure = OC.VerticalScalarDiffusivity(ν = 1e-5, κ = 2e-6)
        catke_closure = CO.Oceans.default_ocean_closure()
        closure = (catke_closure, eddy_closure, horizontal_viscosity, vertical_closure)
    else
        # Simpler setup
        @info "Using simpler ocean setup; to be used for software testing only."
        free_surface = OC.SplitExplicitFreeSurface(grid; substeps = 70)
        momentum_advection = OC.WENOVectorInvariant(order = 5)
        horizontal_viscosity = OC.HorizontalScalarDiffusivity(ν = 1e4)
        tracer_advection = OC.WENO(order = 5)
        vertical_mixing = OC.ConvectiveAdjustmentVerticalDiffusivity(
            background_κz = 1e-5,
            convective_κz = 0.1,
            background_νz = 1e-4,
            convective_νz = 0.1,
        )

        closure = (horizontal_viscosity, vertical_mixing)
    end

    if tspan[1] isa ITime
        # create a model clock that uses DateTime, for compatibility with ITime.
        model_clock = OC.TimeSteppers.Clock(time = start_date)
        stop_time = Dates.DateTime(tspan[2])
    elseif tspan[1] isa Float64
        model_clock = OC.TimeSteppers.Clock{Float64}(time = tspan[1])
        stop_time = tspan[2]
    else
        error("Unsupported time type: $(typeof(tspan[1]))")
    end
    model_Δt = dt
    ocean = CO.ocean_simulation(
        grid;
        clock = model_clock,
        stop_time,
        Δt = float(model_Δt),
        timestepper = :SplitRungeKutta3,
        momentum_advection,
        tracer_advection,
        free_surface,
        closure,
    )

    # Set initial condition to EN4 state estimate at start_date
    OC.set!(ocean.model, T = en4_temperature[1], S = en4_salinity[1])

    # Construct the remapper object and allocate scratch space
    remapping = construct_remapper(grid, boundary_space)

    # COARE3 roughness params (allocated once, reused each timestep)
    coare3_roughness_params = CC.Fields.Field(SF.COARE3RoughnessParams{FT}, boundary_space)
    coare3_roughness_params .= SF.COARE3RoughnessParams{FT}()

    # Get some ocean properties and parameters (including COARE3 roughness params)
    ocean_properties = (;
        reference_density = 1020,
        heat_capacity = 3991,
        σ = coupled_param_dict["stefan_boltzmann_constant"],
        C_to_K = coupled_param_dict["temperature_water_freeze"],
        coare3_roughness_params,
    )

    # Before version 0.96.22, the NetCDFWriter was broken on GPU
    if arch isa OC.CPU
        # Save all tracers and velocities to a NetCDF file at daily frequency
        outputs = merge(ocean.model.tracers, ocean.model.velocities)
        surface_writer = OC.NetCDFWriter(
            ocean.model,
            outputs;
            schedule = OC.TimeInterval(3600), # hourly snapshots
            filename = joinpath(output_dir, "ocean_diagnostics.nc"),
            indices = (:, :, grid.Nz),
            overwrite_existing = true,
            array_type = Array{FT},
        )
        free_surface_writer = OC.NetCDFWriter(
            ocean.model,
            (; displacement = ocean.model.free_surface.displacement);
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
        model_Δt,
    )
end

"""
    construct_remapper(grid_oc, boundary_space)

Given an Oceananigans grid and a ClimaCore boundary space, construct the
remappers needed to remap between the two grids in both directions.

Returns a remapper from the Oceananigans grid to the ClimaCore boundary space.
To regrid from Oceananigans to ClimaCore, use `CR.regrid!(dest_vector, remapper_oc_to_cc, src_vector)`.
To regrid from ClimaCore to Oceananigans, use `CR.regrid!(dest_vector, transpose(remapper_oc_to_cc), src_vector)`.
"""
function construct_remapper(grid_oc, boundary_space)
    # Move grids to CPU since ConservativeRegridding doesn't support GPU grids yet
    grid_oc_underlying_cpu = OC.on_architecture(OC.CPU(), grid_oc.underlying_grid)
    boundary_space_cpu = CC.Adapt.adapt(Array, boundary_space)

    # Create the remapper from the Oceananigans grid to the ClimaCore boundary space
    remapper_oc_to_cc = CR.Regridder(
        boundary_space_cpu,
        grid_oc_underlying_cpu;
        normalize = false,
        threaded = false,
    )

    # Move remapper to GPU if needed
    remapper_oc_to_cc = OC.on_architecture(OC.architecture(grid_oc), remapper_oc_to_cc)

    # Create a field of ones on the boundary space so we can compute element areas
    field_ones_cc = CC.Fields.ones(boundary_space)

    # Allocate a vector with length equal to the number of elements in the target space
    # To be used as a temp field for remapping
    FT = CC.Spaces.undertype(boundary_space)
    ArrayType = ClimaComms.array_type(boundary_space)
    value_per_element_cc =
        ArrayType(zeros(FT, CC.Meshes.nelements(boundary_space.grid.topology.mesh)))

    # Construct 2D Oceananigans Center/Center fields as scratch space while remapping
    scratch_field_oc1 = OC.Field{OC.Center, OC.Center, Nothing}(grid_oc)
    scratch_field_oc2 = OC.Field{OC.Center, OC.Center, Nothing}(grid_oc)
    scratch_field_oc3 = OC.Field{OC.Center, OC.Center, Nothing}(grid_oc)

    # Allocate space for a Field of UVVectors, which we need for remapping momentum fluxes
    temp_uv_vec = CC.Fields.Field(CC.Geometry.UVVector{FT}, boundary_space)

    # Precompute mask for removing ocean surface fluxes at the poles for LatLonGrid
    # and near south pole for TripolarGrid
    # This mask lives on cell centers
    # (Will be removed when we switch to other grids)
    polar_mask = ocean_polar_mask(grid_oc.underlying_grid; location = (OC.Center(), OC.Center(), OC.Center()))

    return (;
        remapper_oc_to_cc,
        field_ones_cc,
        value_per_element_cc,
        scratch_field_oc1,
        scratch_field_oc2,
        scratch_field_oc3,
        temp_uv_vec,
        polar_mask,
    )
end

"""
    FieldExchanger.resolve_area_fractions!(ocean_sim, ice_sim, land_fraction)

Ensure the ocean and ice area fractions are consistent with each other.

The ocean's LatitudeLongitudeGrid is only defined between -80 and 80
degrees latitude. In this case, we set ice and ocean area fractions
to 0 and land to 1 where |lat| ≥ 78° (same band used to zero ocean
surface fluxes).

Similarly, the ocean's TripolarGrid is defined between -80 and 90 degrees latitude.
In this case we set the ice and ocean area fractions to 0 and the land fraction to 1
on [-78°S, -90°S], and leave the northern latitudes unmodified.
"""
function FieldExchanger.resolve_area_fractions!(
    ocean_sim::OceananigansSimulation,
    ice_sim,
    land_fraction,
)
    ocean_fraction = Interfacer.get_field(ocean_sim, Val(:area_fraction))
    ice_fraction = Interfacer.get_field(ice_sim, Val(:area_fraction))

    boundary_space = axes(ocean_fraction)
    FT = CC.Spaces.undertype(boundary_space)
    lat = CC.Fields.coordinate_field(boundary_space).lat
    polar_mask = CC.Fields.ones(boundary_space)

    if ocean_sim.ocean.model.grid.underlying_grid isa OC.LatitudeLongitudeGrid
        # Polar mask: 1 where |lat| < 78° (valid ocean), 0 where polar
        polar_mask .= abs.(lat) .< FT(78)
    elseif ocean_sim.ocean.model.grid.underlying_grid isa OC.TripolarGrid
        # TODO we should be able to remove this polar mask for Tripolar - double check by plotting land/sea mask
        # Create a mask that's 1 where valid (lat >= -78°) and 0 at the south pole
        polar_mask .= lat .>= FT(-78)
    end

    # Set land fraction to 1 and ice/ocean fraction to 0 where polar_mask is 0
    @. land_fraction = ifelse.(polar_mask == FT(0), FT(1), land_fraction)
    @. ice_fraction = ifelse.(polar_mask == FT(0), FT(0), ice_fraction)
    @. ocean_fraction = ifelse.(polar_mask == FT(0), FT(0), ocean_fraction)

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

# Timestep the simulation forward to time `t`. This may not actually do anything.
function Interfacer.step!(sim::OceananigansSimulation, t::Float64)
    # `round(Int, ...)` tolerates floating point drift less than `model_dt / 2`
    n_steps = round(Int, (t - sim.ocean.model.clock.time) / sim.model_Δt)
    for _ in 1:n_steps
        OC.time_step!(sim.ocean, sim.model_Δt)
    end
    return nothing
end

function Interfacer.step!(sim::OceananigansSimulation, t::ITime)
    Δt_msec = date(t) - sim.ocean.model.clock.time
    model_Δt_msec = counter(sim.model_Δt) * Dates.Millisecond(period(sim.model_Δt))
    n_steps = div(Δt_msec, model_Δt_msec) # integer division; exact for Millisecond periods
    for _ in 1:n_steps
        OC.time_step!(sim.ocean, float(sim.model_Δt))
    end
    return nothing
end

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

"""
    ocean_polar_mask(underlying_grid; location)

Build the ocean polar mask once at setup.
Currently we define ocean between 80°S to 80°N with 2 degree overlap in the coupler mask.
Returns a 2D mask (1.0 where |lat| < 78°, 0.0 elsewhere). This mask is on the ocean grid
(unlike the polar mask which is defined on the boundary_space)
"""
function ocean_polar_mask(underlying_grid::OC.LatitudeLongitudeGrid; location = (OC.Center(), OC.Center(), OC.Center()))
    polar_flux_lat_deg = 78.0  # zero fluxes where |lat| ≥ this (same band as polar_mask for atmosphere)

    # latitude nodes: a StepRangeLen of size grid.Ny *in degrees*
    φ = OC.φnodes(underlying_grid, location[1], location[2], location[3])

    # compute mask (1.0 where |lat| < 78°, 0.0 elsewhere)
    mask = ifelse.(abs.(φ) .< polar_flux_lat_deg, 1.0, 0.0)  # Vector of size grid.Ny
    mask = reshape(mask, 1, :)  # make mask a row vector (1 × grid.Ny)
    mask = repeat(mask, underlying_grid.Nx, 1)  # repeat across longitude to get a grid.Nx × grid.Ny Matrix

    # move to architecture
    arch = OC.Architectures.architecture(underlying_grid)
    return OC.Architectures.on_architecture(arch, mask)
end

"""
    ocean_polar_mask(underlying_grid::TripolarGrid; location)

Same convention as the LatLon method: `1.0` where surface fluxes are applied, `0.0` in the
south polar cap (latitudes below about -78°). See `resolve_area_fractions!` for the matching
coupler boundary mask.
"""
function ocean_polar_mask(
    underlying_grid::OC.TripolarGrid;
    location = (OC.Center(), OC.Center(), OC.Center()),
)
    # zero fluxes where lat ≥ this (same band as polar_mask for atmosphere)
    polar_flux_lat_deg = -78.0

    # TODO check if this is correct, esp radians vs degrees and dimensions of arrays (different from LatLonGrid)
    # compute mask (1.0 where lat >= -78°, 0.0 elsewhere)
    φ = OC.φnodes(underlying_grid, location[1], location[2], location[3])
    φ_2D = Array(φ[:, :, 1])
    mask = ifelse.(φ_2D .> polar_flux_lat_deg, 1.0, 0.0)  # Vector of size grid.Ny

    # move to architecture
    architecture = OC.Architectures.architecture(underlying_grid)
    return OC.Architectures.on_architecture(architecture, mask)
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
    ice_concentration_field = sim.ice_concentration
    # for masking out the poles
    polar_mask = sim.remapping.polar_mask

    # Convert the momentum fluxes from contravariant to Cartesian basis
    contravariant_to_cartesian!(sim.remapping.temp_uv_vec, F_turb_ρτxz, F_turb_ρτyz)
    F_turb_ρτxz_uv = sim.remapping.temp_uv_vec.components.data.:1
    F_turb_ρτyz_uv = sim.remapping.temp_uv_vec.components.data.:2

    # Remap momentum fluxes onto reduced 2D Center, Center fields
    Interfacer.remap!(sim.remapping.scratch_field_oc1, F_turb_ρτxz_uv, sim.remapping) # zonal momentum flux
    Interfacer.remap!(sim.remapping.scratch_field_oc2, F_turb_ρτyz_uv, sim.remapping) # meridional momentum flux
    # Rename for clarity; these are now cell-centered (Center, Center) Oceananigans fields
    F_turb_ρτxz_oc = sim.remapping.scratch_field_oc1
    F_turb_ρτyz_oc = sim.remapping.scratch_field_oc2

    # mask out the poles
    @. F_turb_ρτxz_oc = polar_mask * F_turb_ρτxz_oc
    @. F_turb_ρτyz_oc = polar_mask * F_turb_ρτyz_oc

    # Weight by (1 - sea ice concentration)
    ice_concentration = OC.interior(ice_concentration_field, :, :, 1)
    OC.interior(F_turb_ρτxz_oc, :, :, 1) .=
        OC.interior(F_turb_ρτxz_oc, :, :, 1) .* (1.0 .- ice_concentration) ./
        reference_density
    OC.interior(F_turb_ρτyz_oc, :, :, 1) .=
        OC.interior(F_turb_ρτyz_oc, :, :, 1) .* (1.0 .- ice_concentration) ./
        reference_density

    # Set the momentum flux BCs at the correct locations using the remapped scratch fields
    oc_flux_u = surface_flux(sim.ocean.model.velocities.u)
    oc_flux_v = surface_flux(sim.ocean.model.velocities.v)
    set_from_extrinsic_vector!(
        (; u = oc_flux_u, v = oc_flux_v),
        grid,
        F_turb_ρτxz_oc,
        F_turb_ρτyz_oc,
    )

    # Remap the latent and sensible heat fluxes using scratch arrays
    Interfacer.remap!(sim.remapping.scratch_field_oc1, F_lh, sim.remapping) # latent heat flux
    Interfacer.remap!(sim.remapping.scratch_field_oc2, F_sh, sim.remapping) # sensible heat flux

    # Rename for clarity; recall F_turb_energy = F_lh + F_sh
    remapped_F_lh = OC.interior(sim.remapping.scratch_field_oc1, :, :, 1)
    remapped_F_sh = OC.interior(sim.remapping.scratch_field_oc2, :, :, 1)

    # TODO: Note, SW radiation penetrates the surface. Right now, we just put
    # everything on the surface, but later we will need to account for this.
    # One way we can do this is using directly ClimaOcean
    oc_flux_T = surface_flux(sim.ocean.model.tracers.T)
    # mask out the poles
    @. remapped_F_lh = polar_mask * remapped_F_lh
    @. remapped_F_sh = polar_mask * remapped_F_sh
    OC.interior(oc_flux_T, :, :, 1) .+=
        (1.0 .- ice_concentration) .* (remapped_F_lh .+ remapped_F_sh) ./
        (reference_density * heat_capacity)

    # Add the part of the salinity flux that comes from the moisture flux, we also need to
    # add the component due to precipitation (that was done with the radiative fluxes)
    Interfacer.remap!(sim.remapping.scratch_field_oc1, F_turb_moisture, sim.remapping) # moisture flux
    moisture_fresh_water_flux =
        OC.interior(sim.remapping.scratch_field_oc1, :, :, 1) ./ reference_density

    oc_flux_S = surface_flux(sim.ocean.model.tracers.S)
    surface_salinity = OC.interior(sim.ocean.model.tracers.S, :, :, grid.Nz)
    # mask out the poles
    @. moisture_fresh_water_flux = polar_mask * moisture_fresh_water_flux
    OC.interior(oc_flux_S, :, :, 1) .-=
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
    grid = sim.ocean.model.grid
    Nz = grid.Nz
    ice_concentration = OC.interior(sim.ice_concentration, :, :, 1)

    # Reset fluxes to 0 at the start of the step
    oc_flux_T = surface_flux(sim.ocean.model.tracers.T)
    OC.interior(oc_flux_T, :, :, 1) .= 0
    oc_flux_S = surface_flux(sim.ocean.model.tracers.S)
    OC.interior(oc_flux_S, :, :, 1) .= 0

    # for masking out the poles
    polar_mask = sim.remapping.polar_mask

    # Remap shortwave and longwave radiation into separate scratch fields to avoid aliasing
    Interfacer.remap!(sim.remapping.scratch_field_oc1, csf.SW_d, sim.remapping) # shortwave radiation
    remapped_SW_d = OC.interior(sim.remapping.scratch_field_oc1, :, :, 1)

    Interfacer.remap!(sim.remapping.scratch_field_oc2, csf.LW_d, sim.remapping) # longwave radiation
    remapped_LW_d = OC.interior(sim.remapping.scratch_field_oc2, :, :, 1)

    # Update only the part due to radiative fluxes. For the full update, the component due
    # to latent and sensible heat is missing and will be updated in update_turbulent_fluxes.
    (; σ, C_to_K) = sim.ocean_properties
    α = Interfacer.get_field(sim, Val(:surface_direct_albedo)) # scalar
    ϵ = Interfacer.get_field(sim, Val(:emissivity)) # scalar

    # Compute radiative contribution (scratch_field_oc3: remapped_SW_d aliases scratch_field_oc1)
    rad_T_flux = OC.interior(sim.remapping.scratch_field_oc3, :, :, 1)
    rad_T_flux .=
        (1.0 .- ice_concentration) .* (
            -(1 - α) .* remapped_SW_d .-
            ϵ * (
                remapped_LW_d .-
                σ .* (C_to_K .+ OC.interior(sim.ocean.model.tracers.T, :, :, Nz)) .^ 4
            )
        ) ./ (reference_density * heat_capacity)
    # mask out the poles
    @. rad_T_flux = polar_mask * rad_T_flux
    OC.interior(oc_flux_T, :, :, 1) .+= rad_T_flux

    # Update salinity flux with precipitation contribution
    Interfacer.remap!(sim.remapping.scratch_field_oc1, csf.P_liq, sim.remapping) # liquid precipitation
    remapped_P_liq = OC.interior(sim.remapping.scratch_field_oc1, :, :, 1)
    # mask out the poles
    @. remapped_P_liq = polar_mask * remapped_P_liq

    Interfacer.remap!(sim.remapping.scratch_field_oc2, csf.P_snow, sim.remapping) # snow precipitation
    remapped_P_snow = OC.interior(sim.remapping.scratch_field_oc2, :, :, 1)
    @. remapped_P_snow = polar_mask * remapped_P_snow

    # Note the negative sign here to account for the sign change from precipitation to salinity flux
    # polar-exclusion mask
    OC.interior(oc_flux_S, :, :, 1) .-=
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

# Additional OceananigansSimulation getter methods for plotting debug fields
Interfacer.get_field(sim::OceananigansSimulation, ::Val{:salinity}) =
    sim.ocean.model.tracers.S
Interfacer.get_field(sim::OceananigansSimulation, ::Val{:u}) = sim.ocean.model.velocities.u
Interfacer.get_field(sim::OceananigansSimulation, ::Val{:v}) = sim.ocean.model.velocities.v
Interfacer.get_field(sim::OceananigansSimulation, ::Val{:free_surface_displacement}) =
    sim.ocean.model.free_surface.displacement

"""
    Plotting.debug_plot_fields(sim::OceananigansSimulation)

Return the fields to include in debug plots for an Oceananigans simulation.
This includes the area fraction, surface temperature, salinity, velocity, and
free surface displacement. These plots are not polished, and are intended for debugging.
"""
Plotting.debug_plot_fields(sim::OceananigansSimulation) =
    (:area_fraction, :surface_temperature, :salinity, :u, :v, :free_surface_displacement)
