import ClimaComms
import NVTX
import SurfaceFluxes as SF
import Thermodynamics as TD
import ClimaOcean.EN4: download_dataset
import ClimaUtilities.TimeManager: ITime, date, counter, period
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
struct OceananigansSimulation{SIM, A, OPROP, REMAP, SIC, MDT} <:
       Interfacer.AbstractOceanSimulation
    ocean::SIM
    area_fraction::A
    ocean_properties::OPROP
    remapping::REMAP
    ice_concentration::SIC
    model_Δt::MDT
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
- `progress_interval`: iteration interval for printing progress information (default: `nothing`)

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
    progress_interval = nothing,
    ocean_diagnostic_interval = "1days",
    ocean_diagnostic_mode = :average,
    extra_kwargs...,
) where {FT}
    arch = comms_ctx.device isa ClimaComms.CUDADevice ? OC.GPU() : OC.CPU()
    OC.Oceananigans.defaults.FloatType = FT

    # Compute stop_date for oceananigans (needed for EN4 data retrieval)
    stop_date = start_date + Dates.Second(float(tspan[2]))

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
        latitude = (-80, 80),   # NOTE: Don't go too high up when using LatLongGrid, or the cells will be too small
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
        tracer_advection = OC.WENO(order = 5)
        vertical_mixing = OC.ConvectiveAdjustmentVerticalDiffusivity(
            background_κz = 1e-5,
            convective_κz = 0.1,
            background_νz = 1e-4,
            convective_νz = 0.1,
        )
        momentum_advection = OC.WENOVectorInvariant(order = 5)
        horizontal_viscosity = OC.HorizontalScalarDiffusivity(ν = 1e4)
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

    wall_time = Ref(time_ns())

    """
        progress(sim)

    Output the extrema of some prognostic variables, which can be useful for debugging.
    The frequency with which this is output is determined by the interval passed to
    `OC.add_callback!` below.
    """
    function progress(sim)
        ocean = sim.model

        (Tmax, Tmin) = extrema(ocean.tracers.T)
        (Smax, Smin) = extrema(ocean.tracers.S)
        (ηmax, ηmin) = extrema(ocean.free_surface.displacement)
        umax = maximum(abs, ocean.velocities.u)
        vmax = maximum(abs, ocean.velocities.v)
        wmax = maximum(abs, ocean.velocities.w)
        step_time = 1e-9 * (time_ns() - wall_time[])
        @info "time: $(OC.Utils.prettytime(sim)), iteration: $(OC.iteration(sim)), Δt: $(OC.Utils.prettytime(sim.Δt)), " *
              "extrema(η): ($(round(ηmin, sigdigits=2)), $(round(ηmax, sigdigits=2))) " *
              "extrema(T, S): ($(round(Tmin, digits=2)), $(round(Tmax, digits=2))) ᵒC, " *
              "($(round(Smin, digits=2)), $(round(Smax, digits=2))) psu " *
              "maximum(u): ($(round(umax, sigdigits=2)), $(round(vmax, sigdigits=2)), $(round(wmax, sigdigits=2))) m/s, " *
              "wall time: $(OC.Utils.prettytime(step_time))"

        wall_time[] = time_ns()

        return nothing
    end

    # Attaching a progress function to the ocean
    if !isnothing(progress_interval)
        OC.add_callback!(ocean, progress, OC.IterationInterval(progress_interval))
    end

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

    # Initialize with 0 ice concentration; this will be updated in `resolve_area_fractions!`
    # if the ocean is coupled to a non-prescribed sea ice model.
    ice_concentration = OC.Field{OC.Center, OC.Center, Nothing}(grid)

    # Create a dummy area fraction that will get overwritten in `update_surface_fractions!`
    area_fraction = ones(boundary_space)

    sim = OceananigansSimulation(
        ocean,
        area_fraction,
        ocean_properties,
        remapping,
        ice_concentration,
        model_Δt,
    )

    add_ocean_diagnostics!(
        sim;
        output_dir,
        interval = TimeManager.time_to_period(ocean_diagnostic_interval),
        mode = ocean_diagnostic_mode,
    )

    return sim
end

"""
    convert_regridder_eltype(::Type{FT}, regridder) where {FT}

Convert the element type of a ConservativeRegridding.Regridder's internal arrays
to the specified float type `FT`. 
"""
function convert_regridder_eltype(::Type{FT}, regridder) where {FT}
    intersections = regridder.intersections
    dst_areas = regridder.dst_areas
    src_areas = regridder.src_areas
    dst_temp = regridder.dst_temp
    src_temp = regridder.src_temp

    new_intersections = SparseArrays.SparseMatrixCSC(
        intersections.m,
        intersections.n,
        intersections.colptr,
        intersections.rowval,
        FT.(intersections.nzval),
    )

    return CR.Regridder(
        new_intersections,
        FT.(dst_areas),
        FT.(src_areas),
        FT.(dst_temp),
        FT.(src_temp),
    )
end

"""
    construct_remapper(grid_oc, boundary_space)

Given an Oceananigans grid and a ClimaCore boundary space, construct the
remappers needed to remap between the two grids in both directions.

Returns a remapper from the Oceananigans grid to the ClimaCore boundary space.
For low-level use: `CR.regrid!(dst, remapper_oc_to_cc, src)` and 
`CR.regrid!(dst, transpose(remapper_oc_to_cc), src)` for the reverse map.
"""
function construct_remapper(grid_oc, boundary_space)
    grid_oc_underlying_cpu = OC.on_architecture(OC.CPU(), grid_oc.underlying_grid)
    boundary_space_cpu = CC.Adapt.adapt(Array, boundary_space)

    FT_cc = CC.Spaces.undertype(boundary_space_cpu)
    R = CC.Spaces.topology(boundary_space_cpu).mesh.domain.radius
    manifold = CR.Spherical(; radius = FT_cc(R))

    remapper_oc_to_cc = CR.Regridder(
        manifold,
        boundary_space_cpu,
        grid_oc_underlying_cpu;
        normalize = false,
        threaded = false,
    )

    # Convert sparse matrix element type to match the simulation's float type.
    remapper_oc_to_cc = convert_regridder_eltype(FT_cc, remapper_oc_to_cc)

    # Move remapper to GPU if needed
    remapper_oc_to_cc =
        OC.Architectures.on_architecture(OC.architecture(grid_oc), remapper_oc_to_cc)

    field_ones_cc = CC.Fields.ones(boundary_space)

    # Allocate a vector with length equal to the number of elements in the target space
    # To be used as a temp field for remapping
    FT = CC.Spaces.undertype(boundary_space)
    ArrayType = ClimaComms.array_type(boundary_space)
    value_per_element_cc =
        ArrayType(zeros(FT, CC.Meshes.nelements(boundary_space.grid.topology.mesh)))

    # Construct three 2D Oceananigans Center/Center fields to use as scratch space while remapping
    scratch_field_oc1 = OC.Field{OC.Center, OC.Center, Nothing}(grid_oc)
    scratch_field_oc2 = OC.Field{OC.Center, OC.Center, Nothing}(grid_oc)
    scratch_field_oc3 = OC.Field{OC.Center, OC.Center, Nothing}(grid_oc)

    # Allocate space for a Field of UVVectors, which we need for remapping momentum fluxes
    temp_uv_vec = CC.Fields.Field(CC.Geometry.UVVector{FT}, boundary_space)

    # Precompute mask (zero where |lat| ≥ 78°) for removing ocean surface fluxes at the poles
    # This mask lives on cell centers
    # (Will be removed when we switch to other grids)
    polar_mask =
        ocean_polar_mask(grid_oc; location = (OC.Center(), OC.Center(), OC.Center()))

    # Extract intersection grid for intersection-grid flux calculations
    # This uses the CPU regridder since intersection grid iteration is CPU-based
    remapper_cpu = OC.on_architecture(OC.CPU(), remapper_oc_to_cc)
    intersection_grid = extract_intersection_grid(remapper_cpu)

    # Allocate flux state for intersection-grid calculations
    intersection_flux_state = IntersectionFluxState(FT, intersection_grid.n_intersections)

    # Allocate temporary vectors for CC and OC state on intersection grid
    n_cc = intersection_grid.n_cc
    n_oc = intersection_grid.n_oc
    cc_atmos_temp = (
        T = zeros(FT, n_cc),
        q_tot = zeros(FT, n_cc),
        q_liq = zeros(FT, n_cc),
        q_ice = zeros(FT, n_cc),
        ρ = zeros(FT, n_cc),
        u = zeros(FT, n_cc),
        v = zeros(FT, n_cc),
        h = zeros(FT, n_cc),
        h_sfc = zeros(FT, n_cc),
    )
    oc_surface_temp = (
        T = zeros(FT, n_oc),
        z0m = zeros(FT, n_oc),
        z0b = zeros(FT, n_oc),
        h = zeros(FT, n_oc),
    )

    return (;
        remapper_oc_to_cc,
        field_ones_cc,
        value_per_element_cc,
        scratch_field_oc1,
        scratch_field_oc2,
        scratch_field_oc3,
        temp_uv_vec,
        polar_mask,
        intersection_grid,
        intersection_flux_state,
        cc_atmos_temp,
        oc_surface_temp,
    )
end

"""
    ocean_polar_mask(grid; location)

Build the ocean polar mask once at setup.
Currently we define ocean between 80°S to 80°N with 2 degree overlap in the coupler mask.
Returns a 2D mask (1.0 where |lat| < 78°, 0.0 elsewhere). This mask is on the ocean grid
(unlike the polar mask which is defined on the boundary_space)
"""
function ocean_polar_mask(grid; location = (OC.Center(), OC.Center(), OC.Center()))
    polar_flux_lat_deg = 78.0  # zero fluxes where |lat| ≥ this (same band as polar_mask for atmosphere)

    # latitude nodes: a StepRangeLen of size grid.Ny *in degrees*
    φ = OC.φnodes(grid, location[1], location[2], location[3])

    # compute mask (1.0 where |lat| < 78°, 0.0 elsewhere)
    mask = ifelse.(abs.(φ) .< polar_flux_lat_deg, 1.0, 0.0)  # Vector of size grid.Ny
    mask = reshape(mask, 1, :)  # make mask a row vector (1 × grid.Ny)
    mask = repeat(mask, grid.Nx, 1)  # repeat across longitude to get a grid.Nx × grid.Ny Matrix

    # move to architecture
    arch = OC.Architectures.architecture(grid)
    return OC.Architectures.on_architecture(arch, mask)
end

"""
    FieldExchanger.resolve_area_fractions!(ocean_sim, ice_sim, land_fraction)

Ensure the ocean and ice area fractions are consistent with each other.
This matters in the case of a LatitudeLongitudeGrid, which only extends to
±80° latitude. We set ice and ocean area fractions to 0 and land to 1 where
|lat| ≥ 78° (same band used to zero ocean surface fluxes).

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

        # Polar mask: 0 where |lat| ≥ 78° (same band used to zero ocean fluxes)
        boundary_space = axes(ocean_fraction)
        FT = CC.Spaces.undertype(boundary_space)
        lat = CC.Fields.coordinate_field(boundary_space).lat
        polar_mask = CC.Fields.zeros(boundary_space)
        polar_mask .= abs.(lat) .< FT(78)

        # Set land fraction to 1 and ice/ocean fraction to 0 where polar_mask is 0
        @. land_fraction = polar_mask * land_fraction + (1 - polar_mask)
        @. ice_fraction = polar_mask * ice_fraction
        @. ocean_fraction = polar_mask * ocean_fraction
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

# Timestep the simulation forward to time `t`. This may not actually do anything.
NVTX.@annotate function Interfacer.step!(sim::OceananigansSimulation, t::Float64)
    # `round(Int, ...)` tolerates floating point drift less than `model_dt / 2`
    n_steps = round(Int, (t - sim.ocean.model.clock.time) / sim.model_Δt)
    for _ in 1:n_steps
        OC.time_step!(sim.ocean, sim.model_Δt)
    end
    return nothing
end

NVTX.@annotate function Interfacer.step!(sim::OceananigansSimulation, t::ITime)
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
    construct_polar_mask(grid)

Construct the polar exclusion flux masks on the ocean grid.
We define three masks to be used for different fluxes:
- `polar_exclusion_flux_mask_centers`: cell-centered (heat and salinity fluxes)
- `polar_exclusion_flux_mask_u`: Face/Center (zonal momentum flux)
- `polar_exclusion_flux_mask_v`: Center/Face (meridional momentum flux)
"""
function construct_polar_mask(grid)
    polar_exclusion_flux_mask_centers = ocean_flux_highlat_mask(
        grid.underlying_grid;
        location = (OC.Center(), OC.Center(), OC.Center()),
    )
    polar_exclusion_flux_mask_u = ocean_flux_highlat_mask(
        grid.underlying_grid;
        location = (OC.Face(), OC.Center(), OC.Center()),
    )
    polar_exclusion_flux_mask_v = ocean_flux_highlat_mask(
        grid.underlying_grid;
        location = (OC.Center(), OC.Face(), OC.Center()),
    )

    return (;
        polar_exclusion_flux_mask_centers,
        polar_exclusion_flux_mask_u,
        polar_exclusion_flux_mask_v,
    )
end

"""
    ocean_flux_highlat_mask(underlying_grid; location)

Build the ocean polar mask once at setup.
Currently we define ocean between 80°S to 80°N with 2 degree overlap in the coupler mask.
Returns a 2D mask (1.0 where |lat| < 78°, 0.0 elsewhere). This mask is on the ocean grid
(unlike the polar mask which is defined on the boundary_space)
"""
function ocean_flux_highlat_mask(
    underlying_grid::OC.LatitudeLongitudeGrid;
    location = (OC.Center(), OC.Center(), OC.Center()),
)
    polar_flux_lat_deg = 78.0  # zero fluxes where |lat| ≥ this

    # latitude nodes: a StepRangeLen of size grid.Ny *in degrees*
    φ = OC.φnodes(underlying_grid, location[1], location[2], location[3])

    # compute mask (1.0 where |lat| < 78°, 0.0 elsewhere)
    mask = ifelse.(abs.(φ) .< polar_flux_lat_deg, 1.0, 0.0)  # Vector of size grid.Ny
    mask = reshape(mask, 1, :)  # make mask a row vector (1 × grid.Ny)
    mask = repeat(mask, underlying_grid.Nx, 1)  # repeat across longitude to get grid.Nx × grid.Ny

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
    # for masking out the poles (cell-centered fluxes)
    polar_mask = sim.remapping.polar_exclusion_flux_mask_centers

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
    OC.interior(F_turb_ρτxz_oc, :, :, 1) .=
        polar_mask .* OC.interior(F_turb_ρτxz_oc, :, :, 1)
    OC.interior(F_turb_ρτyz_oc, :, :, 1) .=
        polar_mask .* OC.interior(F_turb_ρτyz_oc, :, :, 1)

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

    # for masking out the poles (cell-centered fluxes)
    polar_mask = sim.remapping.polar_exclusion_flux_mask_centers

    # Remap shortwave and longwave onto separate scratch fields
    Interfacer.remap!(sim.remapping.scratch_field_oc3, csf.SW_d, sim.remapping)
    remapped_SW_d = OC.interior(sim.remapping.scratch_field_oc3, :, :, 1)

    Interfacer.remap!(sim.remapping.scratch_field_oc2, csf.LW_d, sim.remapping)
    remapped_LW_d = OC.interior(sim.remapping.scratch_field_oc2, :, :, 1)

    # Update only the part due to radiative fluxes. For the full update, the component due
    # to latent and sensible heat is missing and will be updated in update_turbulent_fluxes.
    (; σ, C_to_K) = sim.ocean_properties
    α = Interfacer.get_field(sim, Val(:surface_direct_albedo)) # scalar
    ϵ = Interfacer.get_field(sim, Val(:emissivity)) # scalar

    # Compute radiative contribution; polar-exclusion mask applied after assignment
    rad_T_flux = OC.interior(sim.remapping.scratch_field_oc1, :, :, 1)
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

    # Remap precipitation fields onto scratch arrays; rename for clarity
    Interfacer.remap!(sim.remapping.scratch_field_oc1, csf.P_liq, sim.remapping) # liquid precipitation
    remapped_P_liq = OC.interior(sim.remapping.scratch_field_oc1, :, :, 1)
    @. remapped_P_liq = ifelse(polar_mask .≈ 0, zero(remapped_P_liq), remapped_P_liq)

    Interfacer.remap!(sim.remapping.scratch_field_oc2, csf.P_snow, sim.remapping) # snow precipitation
    remapped_P_snow = OC.interior(sim.remapping.scratch_field_oc2, :, :, 1)
    @. remapped_P_snow = ifelse(polar_mask .≈ 0, zero(remapped_P_snow), remapped_P_snow)

    # Note the negative sign here to account for the sign change from precipitation to salinity flux
    OC.interior(oc_flux_S, :, :, 1) .-=
        OC.interior(sim.ocean.model.tracers.S, :, :, Nz) .* (1.0 .- ice_concentration) .*
        (remapped_P_liq .+ remapped_P_snow) ./ reference_density
    return nothing
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

# Intersection-grid flux calculation methods

"""
    extract_cc_atmos_state!(cc_atmos_temp, csf, boundary_space)

Extract atmosphere state from coupler fields into per-CC-element vectors.

The coupler fields are stored as ClimaCore Fields on the boundary space.
This function extracts the per-element mean values into vectors suitable
for intersection-grid flux calculations.
"""
function extract_cc_atmos_state!(cc_atmos_temp::NamedTuple, csf, boundary_space)
    CRExt = get_ConservativeRegriddingCCExt()
    field_ones = CC.Fields.ones(boundary_space)

    CRExt.get_value_per_element!(cc_atmos_temp.T, csf.T_atmos, field_ones)
    CRExt.get_value_per_element!(cc_atmos_temp.q_tot, csf.q_tot_atmos, field_ones)
    CRExt.get_value_per_element!(cc_atmos_temp.q_liq, csf.q_liq_atmos, field_ones)
    CRExt.get_value_per_element!(cc_atmos_temp.q_ice, csf.q_ice_atmos, field_ones)
    CRExt.get_value_per_element!(cc_atmos_temp.ρ, csf.ρ_atmos, field_ones)
    CRExt.get_value_per_element!(cc_atmos_temp.u, csf.u_int, field_ones)
    CRExt.get_value_per_element!(cc_atmos_temp.v, csf.v_int, field_ones)
    CRExt.get_value_per_element!(cc_atmos_temp.h, csf.height_int, field_ones)
    CRExt.get_value_per_element!(cc_atmos_temp.h_sfc, csf.height_sfc, field_ones)
    return nothing
end

"""
    extract_oc_surface_state!(oc_surface_temp, sim::OceananigansSimulation)

Extract ocean surface state into per-OC-cell vectors for intersection-grid flux calculations.
"""
function extract_oc_surface_state!(oc_surface_temp::NamedTuple, sim::OceananigansSimulation)
    grid = sim.ocean.model.grid
    Nz = size(grid, 3)

    # Surface temperature (convert from Celsius to Kelvin)
    C_to_K = sim.ocean_properties.C_to_K
    T_oc = OC.interior(sim.ocean.model.tracers.T, :, :, Nz)
    oc_surface_temp.T .= vec(T_oc) .+ C_to_K

    # Use constant roughness for ocean (will be replaced by COARE3 in future)
    FT = eltype(oc_surface_temp.T)
    fill!(oc_surface_temp.z0m, FT(1e-4))
    fill!(oc_surface_temp.z0b, FT(1e-4))
    fill!(oc_surface_temp.h, FT(0))

    return nothing
end

"""
    compute_intersection_grid_fluxes!(sim::OceananigansSimulation, csf, surface_fluxes_params, thermo_params)

Compute turbulent fluxes on the intersection grid for the ocean simulation.

This provides better coastline representation by computing fluxes directly on the
intersection polygons between CC elements and OC cells, rather than remapping
surface properties to the CC grid before flux computation.

# Arguments
- `sim`: OceananigansSimulation
- `csf`: Coupler fields containing atmosphere state on the boundary space
- `surface_fluxes_params`: SurfaceFluxes parameters
- `thermo_params`: Thermodynamics parameters

After calling this function, the intersection-grid fluxes are stored in
`sim.remapping.intersection_flux_state` and can be scattered to CC or OC grids
using `scatter_flux_to_cc!` or `scatter_flux_to_oc!`.
"""
function compute_intersection_grid_fluxes!(
    sim::OceananigansSimulation,
    csf,
    surface_fluxes_params,
    thermo_params,
)
    (; intersection_grid, intersection_flux_state, cc_atmos_temp, oc_surface_temp) = sim.remapping
    boundary_space = axes(csf)

    # Extract atmosphere state from coupler fields to per-element vectors
    extract_cc_atmos_state!(cc_atmos_temp, csf, boundary_space)

    # Prepare cc_atmos_state NamedTuple (using h_sfc for surface height)
    cc_atmos_state = (
        T = cc_atmos_temp.T,
        q_tot = cc_atmos_temp.q_tot,
        q_liq = cc_atmos_temp.q_liq,
        q_ice = cc_atmos_temp.q_ice,
        ρ = cc_atmos_temp.ρ,
        u = cc_atmos_temp.u,
        v = cc_atmos_temp.v,
        h = cc_atmos_temp.h,
    )

    # Extract ocean surface state to per-cell vectors
    extract_oc_surface_state!(oc_surface_temp, sim)

    oc_surface_state = (
        T = oc_surface_temp.T,
        z0m = oc_surface_temp.z0m,
        z0b = oc_surface_temp.z0b,
        h = oc_surface_temp.h,
    )

    # Compute fluxes on intersection grid
    compute_surface_fluxes_on_intersection!(
        intersection_flux_state,
        intersection_grid,
        cc_atmos_state,
        oc_surface_state,
        surface_fluxes_params,
        thermo_params,
    )

    return nothing
end

"""
    get_intersection_grid(sim::OceananigansSimulation)

Return the intersection grid for an OceananigansSimulation.
"""
get_intersection_grid(sim::OceananigansSimulation) = sim.remapping.intersection_grid

"""
    get_intersection_flux_state(sim::OceananigansSimulation)

Return the intersection flux state for an OceananigansSimulation.
"""
get_intersection_flux_state(sim::OceananigansSimulation) = sim.remapping.intersection_flux_state
