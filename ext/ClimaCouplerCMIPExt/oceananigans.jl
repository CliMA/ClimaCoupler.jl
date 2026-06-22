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
    tripolar_ocean_simulation(arch; active_cells_map, kwargs...)

One-degree tripolar ocean setup matching
`ClimaOcean.OceanConfigurations.one_degree_tripolar_ocean`, but exposing
`active_cells_map`. Tripolar immersed grids with `active_cells_map = true`
exceed the 4 KiB CUDA kernel parameter limit on sm_60 (P100) even with a
minimal closure; disabling the map keeps the tripolar grid while staying
within the limit.
"""
function tripolar_ocean_simulation(arch;
    active_cells_map = true,
    zstar = true,
    Nz = 32,
    depth = 5500,
    momentum_advection = OC.WENOVectorInvariant(order = 5),
    tracer_advection = OC.WENO(order = 5),
    closure = nothing,
    halo = (5, 5, 4),
    minimum_depth = 10,
    interpolation_passes = 10,
    substeps = 70,
    kwargs...)
    z = OC.ExponentialDiscretization(Nz, -depth, 0; mutable = zstar)

    grid = OC.TripolarGrid(arch; size = (360, 180, Nz), z, halo)

    bottom_height = CO.regrid_bathymetry(grid;
        minimum_depth,
        major_basins = 2,
        interpolation_passes)

    grid = OC.ImmersedBoundaryGrid(grid, OC.GridFittedBottom(bottom_height);
        active_cells_map)

    free_surface = OC.SplitExplicitFreeSurface(grid; substeps)

    return CO.ocean_simulation(grid;
        momentum_advection,
        tracer_advection,
        free_surface,
        closure,
        kwargs...)
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
    en4_temperature = CO.Metadata(:temperature; dates, dataset = CO.EN4Monthly())
    en4_salinity = CO.Metadata(:salinity; dates, dataset = CO.EN4Monthly())
    CO.download_with_fallback(en4_temperature)
    CO.download_with_fallback(en4_salinity)

    # Create ocean simulation
    closure = if simple_ocean
        @info "Using simpler ocean setup; to be used for software testing only."
        CO.OceanConfigurations.simplified_ocean_closure()
    else
        CO.OceanConfigurations.default_one_degree_closure()
    end

    # Tripolar active_cells_map kernels exceed the 4 KiB sm_60 limit on P100.
    active_cells_map = !simple_ocean

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

    ocean = tripolar_ocean_simulation(
        arch;
        active_cells_map,
        clock = model_clock,
        depth = 5500,
        Nz = 32,
        closure,
        substeps = simple_ocean ? 70 : 150,
    )
    ocean.stop_time = stop_time
    ocean.Δt = float(dt)

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
    grid = ocean.model.grid
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
        ocean.Δt,
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
    new_intersections = SparseArrays.SparseMatrixCSC{FT}(regridder.intersections)
    return CR.Regridder(
        new_intersections,
        FT.(regridder.dst_areas),
        FT.(regridder.src_areas),
        FT.(regridder.dst_temp),
        FT.(regridder.src_temp),
    )
end

"""
    underlying_grid(grid)

Strip an `ImmersedBoundaryGrid` wrapper if present; `ConservativeRegridding`
operates on the underlying lat/lon (or tripolar / orthogonal-shell) grid and
is unaware of the immersed mask. Raw grids pass through unchanged so that
`construct_remapper` is callable for any grid type accepted by the
`Interfacer.remap` method dispatch in `climaocean_helpers.jl`.
"""
underlying_grid(grid::OC.ImmersedBoundaryGrid) = grid.underlying_grid
underlying_grid(grid) = grid


"""
    construct_remapper(grid_oc, boundary_space)

Given an Oceananigans grid and a ClimaCore boundary space, construct the two
independent sparse regridders needed to remap between them in both directions.

* `remapper_oc_to_cc` — FV → SE, built via the per-element L2 projection
  (`fv_to_se_l2_projection` in `ConservativeRegriddingClimaCoreExt`). The
  inverse element mass matrix `Mᵉ⁻¹` is baked into the sparse matrix and the
  SE finalizer applies `Spaces.weighted_dss!` to reconcile shared nodes.
* `remapper_cc_to_oc` — SE → FV, built via the principled polygon-intersection
  operator (`se_to_fv_principled`). Each row integrates the SE basis over the
  FV intersection polygon and divides by the FV cell area, giving a mean-
  preserving cell average that preserves constants exactly.

For low-level use: `CR.regrid!(dst, remapper_oc_to_cc, src)` and
`CR.regrid!(dst, remapper_cc_to_oc, src)`. The `Interfacer.remap!` methods in
`climaocean_helpers.jl` accept `CC.Fields.Field` directly as source or
destination and route through `ConservativeRegriddingClimaCoreExt`'s nodal
extract / finalize overrides — no per-element scratch buffer is required.

In addition to the two regridders this function also allocates the
intersection-grid flux pipeline used by `compute_intersection_grid_fluxes!`:

* `intersection_grid::IntersectionGrid` — CC-element × OC-cell sparse
  intersection (constructed via [`extract_intersection_grid`](@ref), which
  uses CR's tree-level `intersection_areas`). On a `TripolarGrid` the fold
  row is handled by CR's `PaddedTreeWrapper` automatically.
* `intersection_flux_state::IntersectionFluxState` — per-polygon scratch
  for atmosphere/surface state plus the five computed fluxes
  (`flux_sh`, `flux_lh`, `flux_τx`, `flux_τy`, `flux_evap`).
* `cc_atmos_temp::NamedTuple` — per-SE-node scratch written by
  `extract_cc_atmos_state!` from the coupler fields (flat GLL layout).
* `oc_surface_temp::NamedTuple` — per-OC-cell scratch written by
  `extract_oc_surface_state!` from the live ocean state.

All four are allocated on the architecture of `grid_oc` so they can be
written in place from device-resident `OC.interior(...)` views and
participate in GPU kernels without a host bounce.
"""
function construct_remapper(grid_oc, boundary_space)
    grid_oc_underlying_cpu = OC.on_architecture(OC.CPU(), underlying_grid(grid_oc))
    boundary_space_cpu = CC.Adapt.adapt(Array, boundary_space)

    FT_cc = CC.Spaces.undertype(boundary_space_cpu)
    R = CC.Spaces.topology(boundary_space_cpu).mesh.domain.radius
    manifold = CR.Spherical(; radius = FT_cc(R))

    # FV → SE per-element L2 projection
    remapper_oc_to_cc = CR.Regridder(
        manifold,
        boundary_space_cpu,
        grid_oc_underlying_cpu;
        normalize = false,
        threaded = false,
    )

    # SE → FV principled polygon-intersection operator
    remapper_cc_to_oc = CR.Regridder(
        manifold,
        grid_oc_underlying_cpu,
        boundary_space_cpu;
        normalize = false,
        threaded = false,
    )

    # Convert sparse matrix element types to match the simulation's float type,
    # then move both remappers to the architecture of `grid_oc`.
    remapper_oc_to_cc = convert_regridder_eltype(FT_cc, remapper_oc_to_cc)
    remapper_cc_to_oc = convert_regridder_eltype(FT_cc, remapper_cc_to_oc)

    arch = OC.architecture(grid_oc)
    remapper_oc_to_cc = OC.Architectures.on_architecture(arch, remapper_oc_to_cc)
    remapper_cc_to_oc = OC.Architectures.on_architecture(arch, remapper_cc_to_oc)

    FT = CC.Spaces.undertype(boundary_space)

    # Construct three 2D Oceananigans Center/Center fields to use as scratch space while remapping
    scratch_field_oc1 = OC.Field{OC.Center, OC.Center, Nothing}(grid_oc)
    scratch_field_oc2 = OC.Field{OC.Center, OC.Center, Nothing}(grid_oc)
    scratch_field_oc3 = OC.Field{OC.Center, OC.Center, Nothing}(grid_oc)

    # Allocate space for a Field of UVVectors, which we need for remapping momentum fluxes
    temp_uv_vec = CC.Fields.Field(CC.Geometry.UVVector{FT}, boundary_space)

    # ----- Intersection-grid flux pipeline -----
    #
    # Build the CC-element × OC-cell intersection grid + scratch buffers
    # used by `compute_intersection_grid_fluxes!`. The intersection matrix
    # is computed once here via the tree-level `CR.intersection_areas`
    # (correctly element-indexed on CC and cell-indexed on FV — distinct
    # from `remapper_*` whose sparse matrices are SE-*node*-indexed), and
    # on a `TripolarGrid` the tree builder dispatches to CR's
    # `PaddedTreeWrapper` so the fold row is handled implicitly. The
    # backing arrays live on the OC architecture so `extract_oc_surface_state!`
    # can do an in-place `vec(OC.interior(T_oc))`-style copy.
    intersection_grid = extract_intersection_grid(boundary_space, grid_oc)
    if !(arch isa OC.CPU)
        AT = arr -> OC.on_architecture(arch, arr)
        intersection_grid = IntersectionGrid(
            AT(intersection_grid.cc_indices),
            AT(intersection_grid.oc_indices),
            AT(intersection_grid.areas),
            AT(intersection_grid.cc_areas),
            AT(intersection_grid.oc_areas),
            intersection_grid.n_cc,
            intersection_grid.n_oc,
            intersection_grid.n_intersections,
            intersection_grid.n_nodes,
            AT(intersection_grid.node_gather_polygon),
            AT(intersection_grid.node_gather_node),
            AT(intersection_grid.node_gather_weight),
        )
    end

    _alloc(n) = OC.on_architecture(arch, zeros(FT, n))
    n_cc = intersection_grid.n_cc
    n_int = intersection_grid.n_intersections
    n_nodes = intersection_grid.n_nodes

    Nx_oc, Ny_oc, _ = size(grid_oc)
    n_oc_layout = Nx_oc * Ny_oc

    intersection_flux_state = IntersectionFluxState(
        _alloc(n_int), _alloc(n_int), _alloc(n_int), _alloc(n_int),
        _alloc(n_int), _alloc(n_int), _alloc(n_int), _alloc(n_int),
        _alloc(n_int), _alloc(n_int), _alloc(n_int), _alloc(n_int),
        _alloc(n_int), _alloc(n_int), _alloc(n_int), _alloc(n_int),
        _alloc(n_int), _alloc(n_int),
    )

    # Per-SE-node scratch for atmosphere state. `extract_cc_atmos_state!`
    # writes the nine fields below via `se_field_to_vec`.
    cc_atmos_temp = (
        T     = _alloc(n_nodes),
        q_tot = _alloc(n_nodes),
        q_liq = _alloc(n_nodes),
        q_ice = _alloc(n_nodes),
        ρ     = _alloc(n_nodes),
        u     = _alloc(n_nodes),
        v     = _alloc(n_nodes),
        h     = _alloc(n_nodes),
        h_sfc = _alloc(n_nodes),
    )

    # Per-OC-cell scratch for ocean surface state. `extract_oc_surface_state!`
    oc_surface_temp = (
        T   = _alloc(n_oc_layout),
        z0m = _alloc(n_oc_layout),
        z0b = _alloc(n_oc_layout),
        h   = _alloc(n_oc_layout),
    )

    return (;
        remapper_oc_to_cc,
        remapper_cc_to_oc,
        scratch_field_oc1,
        scratch_field_oc2,
        scratch_field_oc3,
        temp_uv_vec,
        intersection_grid,
        intersection_flux_state,
        cc_atmos_temp,
        oc_surface_temp,
    )
end

"""
    FieldExchanger.resolve_area_fractions!(ocean_sim, ice_sim, land_fraction)

Sync the ice concentration field on the ocean simulation with the ice
simulation's concentration so that it can be used for weighting flux updates.

With a `TripolarGrid` covering the full sphere there is no polar band to
exclude from physical contribution; previously this routine forced land
fraction to 1 (and ice / ocean fractions to 0) for `|lat| ≥ 78°` on a
`LatitudeLongitudeGrid`, which is no longer needed.
"""
function FieldExchanger.resolve_area_fractions!(
    ocean_sim::OceananigansSimulation,
    ice_sim,
    land_fraction,
)
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

    # Weight by (1 - sea ice concentration)
    ice_concentration = OC.interior(ice_concentration_field,:,:,1)
    OC.interior(F_turb_ρτxz_oc,:,:,1) .=
        OC.interior(F_turb_ρτxz_oc,:,:,1) .* (1.0 .- ice_concentration) ./ reference_density
    OC.interior(F_turb_ρτyz_oc,:,:,1) .=
        OC.interior(F_turb_ρτyz_oc,:,:,1) .* (1.0 .- ice_concentration) ./ reference_density

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
    remapped_F_lh = OC.interior(sim.remapping.scratch_field_oc1,:,:,1)
    remapped_F_sh = OC.interior(sim.remapping.scratch_field_oc2,:,:,1)

    # TODO: Note, SW radiation penetrates the surface. Right now, we just put
    # everything on the surface, but later we will need to account for this.
    oc_flux_T = surface_flux(sim.ocean.model.tracers.T)
    OC.interior(oc_flux_T, :, :, 1) .+=
        (1.0 .- ice_concentration) .* (remapped_F_lh .+ remapped_F_sh) ./
        (reference_density * heat_capacity)

    # Add the part of the salinity flux that comes from the moisture flux, we also need to
    # add the component due to precipitation (that was done with the radiative fluxes)
    Interfacer.remap!(sim.remapping.scratch_field_oc1, F_turb_moisture, sim.remapping) # moisture flux
    moisture_fresh_water_flux =
        OC.interior(sim.remapping.scratch_field_oc1,:,:,1) ./ reference_density

    oc_flux_S = surface_flux(sim.ocean.model.tracers.S)
    surface_salinity = OC.interior(sim.ocean.model.tracers.S, :, :, grid.Nz)
    OC.interior(oc_flux_S, :, :, 1) .-=
        (1.0 .- ice_concentration) .* surface_salinity .* moisture_fresh_water_flux
    return nothing
end

function Interfacer.update_field!(sim::OceananigansSimulation, ::Val{:area_fraction}, field)
    sim.area_fraction .= field
    return nothing
end

"""
    update_turbulent_fluxes_intersection_grid!(sim::OceananigansSimulation, csf,
                                               surface_fluxes_params, thermo_params)


Caution: Claude assisted refactor - not yet reviewed!!

Opt-in counterpart to
[`FluxCalculator.update_turbulent_fluxes!(::OceananigansSimulation, fields)`](@ref)
that performs the full *intersection-grid* flux exchange end-to-end:

1. `compute_intersection_grid_fluxes!` computes one SurfaceFluxes call per
   CC↔OC intersection polygon, against per-CC-element atmosphere state
   (taken from `csf`) and per-OC-cell surface state (taken from the live
   ocean tracer fields). Dry polygons are already absent from
   `sim.remapping.intersection_grid` because `construct_remapper` honors
   the immersed mask via `extract_intersection_grid(...,
   respect_immersed_mask = true)`.
2. `scatter_to_oc!` area-averages the polygon fluxes back to the live
   `(Nx, Ny)` OC layout. This is the *exact* area-weighted mean over
   the polygon partition of each cell — it conserves
   `Σ_poly flux × area` and is gradient-preserving at coastlines, in
   contrast to the SE → FV principled regrid which blends across
   element boundaries unconditionally.
3. The OC-aggregated fluxes are then pushed into the ocean's surface
   BCs exactly as the production path does, including:
   * sea-ice-concentration weighting `(1 - α)`,
   * extrinsic-vector rotation of `(F_τx, F_τy)` (in lon/lat / UV basis,
     because `SurfaceFluxes` returns `ρτxz, ρτyz` in the same frame as
     its input wind, and `csf.u_int / v_int` are stored in extrinsic
     atmos basis) into the OC's intrinsic curvilinear basis via
     `set_from_extrinsic_vector!`,
   * unit-conversion divisions by `reference_density` and
     `reference_density * heat_capacity`,
   * salinity sign flip from upward evaporation to downward salinity
     tendency.

# Frame note

Unlike the standard path, *no* `contravariant_to_cartesian!` step is
required here. The standard path receives `F_turb_ρτxz / F_turb_ρτyz`
from ClimaAtmos as *contravariant* components on the cubed sphere; the
intersection-grid path consumes `csf.u_int / v_int` (which are already
stored as `UVVector` lon/lat components on the boundary space) and
passes them through `SurfaceFluxes`, so the polygon momentum fluxes
come back in lon/lat too.

# Vector-rotation caveat

For non-`LatitudeLongitudeGrid` ocean grids (e.g. `TripolarGrid`),
`set_from_extrinsic_vector!` handles the lon/lat → curvilinear-intrinsic
rotation correctly on a per-cell basis. The intra-polygon variation of
the rotation matrix is *not* applied — i.e. all polygons under OC cell
`j` use cell `j`'s rotation. At realistic OC resolutions this is far
below SF discretization noise, but if it ever became material the fix
would be to rotate `(flux_τx, flux_τy)` from lon/lat → cell intrinsic
inside `compute_surface_fluxes_on_intersection!` before aggregating.
"""
function update_turbulent_fluxes_intersection_grid!(
    sim::OceananigansSimulation,
    csf,
    surface_fluxes_params,
    thermo_params,
)
    # Step 1: compute per-polygon fluxes (writes into
    # `sim.remapping.intersection_flux_state`).
    compute_intersection_grid_fluxes!(sim, csf, surface_fluxes_params, thermo_params)

    (; intersection_grid, intersection_flux_state) = sim.remapping
    (; reference_density, heat_capacity) = sim.ocean_properties
    grid = sim.ocean.model.grid
    ice_concentration_field = sim.ice_concentration

    Nx_oc, Ny_oc, _ = size(grid)
    n_oc_layout = Nx_oc * Ny_oc
    FT = eltype(intersection_flux_state.flux_sh)

    # Step 2: aggregate polygon fluxes to the live `(Nx, Ny)` OC layout.
    # `intersection_grid.oc_indices` ranges over `1:n_oc_layout` (we
    # sized `oc_surface_temp` to `Nx*Ny` for exactly this reason — see
    # `construct_remapper`).
    F_sh_oc   = OC.on_architecture(OC.architecture(grid), zeros(FT, n_oc_layout))
    F_lh_oc   = OC.on_architecture(OC.architecture(grid), zeros(FT, n_oc_layout))
    F_τx_oc   = OC.on_architecture(OC.architecture(grid), zeros(FT, n_oc_layout))
    F_τy_oc   = OC.on_architecture(OC.architecture(grid), zeros(FT, n_oc_layout))
    F_evap_oc = OC.on_architecture(OC.architecture(grid), zeros(FT, n_oc_layout))
    scatter_to_oc!(F_sh_oc,   intersection_grid, intersection_flux_state.flux_sh)
    scatter_to_oc!(F_lh_oc,   intersection_grid, intersection_flux_state.flux_lh)
    scatter_to_oc!(F_τx_oc,   intersection_grid, intersection_flux_state.flux_τx)
    scatter_to_oc!(F_τy_oc,   intersection_grid, intersection_flux_state.flux_τy)
    scatter_to_oc!(F_evap_oc, intersection_grid, intersection_flux_state.flux_evap)

    # Reshape to the 2D layout the OC surface fields expect.
    F_sh_oc_2d   = reshape(F_sh_oc,   Nx_oc, Ny_oc)
    F_lh_oc_2d   = reshape(F_lh_oc,   Nx_oc, Ny_oc)
    F_τx_oc_2d   = reshape(F_τx_oc,   Nx_oc, Ny_oc)
    F_τy_oc_2d   = reshape(F_τy_oc,   Nx_oc, Ny_oc)
    F_evap_oc_2d = reshape(F_evap_oc, Nx_oc, Ny_oc)

    ice_concentration = OC.interior(ice_concentration_field, :, :, 1)
    one_minus_ice = 1.0 .- ice_concentration

    # Step 3a: momentum BCs.
    # Stage the lon/lat-frame, ice-weighted, density-normalized fluxes
    # into the same scratch Center,Center fields the standard path uses,
    # then rotate to OC intrinsic via `set_from_extrinsic_vector!`.
    OC.interior(sim.remapping.scratch_field_oc1, :, :, 1) .=
        F_τx_oc_2d .* one_minus_ice ./ reference_density
    OC.interior(sim.remapping.scratch_field_oc2, :, :, 1) .=
        F_τy_oc_2d .* one_minus_ice ./ reference_density
    oc_flux_u = surface_flux(sim.ocean.model.velocities.u)
    oc_flux_v = surface_flux(sim.ocean.model.velocities.v)
    set_from_extrinsic_vector!(
        (; u = oc_flux_u, v = oc_flux_v),
        grid,
        sim.remapping.scratch_field_oc1,
        sim.remapping.scratch_field_oc2,
    )

    # Step 3b: T surface BC from sensible + latent heat. Sign convention
    # matches the standard path: SurfaceFluxes returns SH / LH positive
    # upward and the ocean sees an upward heat loss as a *positive*
    # temperature-tendency contribution. (Radiation and precipitation
    # already filled `oc_flux_T` upstream via `update_sim!`; we add.)
    oc_flux_T = surface_flux(sim.ocean.model.tracers.T)
    OC.interior(oc_flux_T, :, :, 1) .+=
        one_minus_ice .* (F_lh_oc_2d .+ F_sh_oc_2d) ./
        (reference_density * heat_capacity)

    # Step 3c: salinity BC from moisture (evaporation).
    # The OC convention is that moisture leaving the surface (positive
    # `F_evap`) raises near-surface salinity, which appears as a
    # *negative* salinity-flux contribution at the bottom of the surface
    # tendency: `-S_sfc · (E / ρ₀)`.
    moisture_fresh_water_flux = F_evap_oc_2d ./ reference_density
    oc_flux_S = surface_flux(sim.ocean.model.tracers.S)
    surface_salinity = OC.interior(sim.ocean.model.tracers.S, :, :, grid.Nz)
    OC.interior(oc_flux_S, :, :, 1) .-=
        one_minus_ice .* surface_salinity .* moisture_fresh_water_flux

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
    ice_concentration = OC.interior(sim.ice_concentration,:,:,1)

    # Reset fluxes to 0 at the start of the step
    oc_flux_T = surface_flux(sim.ocean.model.tracers.T)
    OC.interior(oc_flux_T,:,:,1) .= 0
    oc_flux_S = surface_flux(sim.ocean.model.tracers.S)
    OC.interior(oc_flux_S, :, :, 1) .= 0

    # Remap shortwave and longwave onto separate scratch fields
    Interfacer.remap!(sim.remapping.scratch_field_oc3, csf.SW_d, sim.remapping)
    remapped_SW_d = OC.interior(sim.remapping.scratch_field_oc3,:,:,1)

    Interfacer.remap!(sim.remapping.scratch_field_oc2, csf.LW_d, sim.remapping)
    remapped_LW_d = OC.interior(sim.remapping.scratch_field_oc2,:,:,1)

    # Update only the part due to radiative fluxes. For the full update, the component due
    # to latent and sensible heat is missing and will be updated in update_turbulent_fluxes.
    (; σ, C_to_K) = sim.ocean_properties
    α = Interfacer.get_field(sim, Val(:surface_direct_albedo)) # scalar
    ϵ = Interfacer.get_field(sim, Val(:emissivity)) # scalar

    # Compute radiative contribution.
    rad_T_flux = OC.interior(sim.remapping.scratch_field_oc1, :, :, 1)
    rad_T_flux .=
        (1.0 .- ice_concentration) .* (
            -(1 - α) .* remapped_SW_d .-
            ϵ * (
                remapped_LW_d .-
                σ .* (C_to_K .+ OC.interior(sim.ocean.model.tracers.T,:,:,Nz)) .^ 4
            )
        ) ./ (reference_density * heat_capacity)
    OC.interior(oc_flux_T, :, :, 1) .+= rad_T_flux

    # Remap precipitation fields onto scratch arrays; rename for clarity
    Interfacer.remap!(sim.remapping.scratch_field_oc1, csf.P_liq, sim.remapping) # liquid precipitation
    remapped_P_liq = OC.interior(sim.remapping.scratch_field_oc1, :, :, 1)

    Interfacer.remap!(sim.remapping.scratch_field_oc2, csf.P_snow, sim.remapping) # snow precipitation
    remapped_P_snow = OC.interior(sim.remapping.scratch_field_oc2, :, :, 1)

    # Note the negative sign here to account for the sign change from precipitation to salinity flux
    OC.interior(oc_flux_S,:,:,1) .-=
        OC.interior(sim.ocean.model.tracers.S,:,:,Nz) .* (1.0 .- ice_concentration) .*
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

Extract atmosphere state from coupler fields into flat nodal vectors suitable
for polygon-averaged intersection-grid flux calculations.
"""
function extract_cc_atmos_state!(cc_atmos_temp::NamedTuple, csf, boundary_space)
    CRExt = get_ConservativeRegriddingCCExt()

    cc_atmos_temp.T     .= CRExt.se_field_to_vec(csf.T_atmos)
    cc_atmos_temp.q_tot .= CRExt.se_field_to_vec(csf.q_tot_atmos)
    cc_atmos_temp.q_liq .= CRExt.se_field_to_vec(csf.q_liq_atmos)
    cc_atmos_temp.q_ice .= CRExt.se_field_to_vec(csf.q_ice_atmos)
    cc_atmos_temp.ρ     .= CRExt.se_field_to_vec(csf.ρ_atmos)
    cc_atmos_temp.u     .= CRExt.se_field_to_vec(csf.u_int)
    cc_atmos_temp.v     .= CRExt.se_field_to_vec(csf.v_int)
    cc_atmos_temp.h     .= CRExt.se_field_to_vec(csf.height_int)
    cc_atmos_temp.h_sfc .= CRExt.se_field_to_vec(csf.height_sfc)
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

    # Extract atmosphere state from coupler fields to per-node vectors
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
