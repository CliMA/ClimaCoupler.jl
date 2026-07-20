import ClimaComms
import NVTX
import SurfaceFluxes as SF
import Thermodynamics as TD
import ClimaUtilities.TimeManager: ITime, date, counter, period
import ClimaUtilities.ClimaArtifacts: @clima_artifact
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
    tripolar_ocean_simulation(arch; clock, stop_time, active_cells_map, kwargs...)

One-degree tripolar ocean setup matching
`ClimaOcean.OceanConfigurations.one_degree_tripolar_ocean`, but exposing
`active_cells_map`. Tripolar immersed grids with `active_cells_map = true`
exceed the 4 KiB CUDA kernel parameter limit on sm_60 (P100) even with a
minimal closure; disabling the map keeps the tripolar grid while staying
within the limit.

`clock` and `stop_time` are required and forwarded to `ClimaOcean.ocean_simulation`;
the coupler builds them from `tspan` so the ocean steps in calendar time when coupled.
"""
function tripolar_ocean_simulation(
    arch;
    clock,
    stop_time,
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
    kwargs...,
)
    z = OC.ExponentialDiscretization(Nz, -depth, 0; mutable = zstar)

    grid = OC.TripolarGrid(arch; size = (360, 180, Nz), z, halo)

    bottom_height =
        CO.regrid_bathymetry(grid; minimum_depth, major_basins = 2, interpolation_passes)

    grid =
        OC.ImmersedBoundaryGrid(grid, OC.GridFittedBottom(bottom_height); active_cells_map)

    free_surface = OC.SplitExplicitFreeSurface(grid; substeps)

    return CO.ocean_simulation(
        grid;
        clock,
        stop_time,
        momentum_advection,
        tracer_advection,
        free_surface,
        closure,
        kwargs...,
    )
end

"""
    ocean_simulation(arch, ocean_grid; simple_ocean, closure, clock, depth, kwargs...)

Build an Oceananigans `Simulation` on either the standard one-degree tripolar grid
or the NEMO eORCA1 mesh (`orca`).

`depth` defaults to 5500 m, matching the EN4 vertical extent.
"""
function ocean_simulation(
    arch,
    ocean_grid::Symbol;
    simple_ocean,
    closure,
    clock,
    depth = 5500,
    kwargs...,
)
    substeps = simple_ocean ? 70 : 150

    if ocean_grid == :orca
        if simple_ocean
            @warn "`simple_ocean=true` not supported for ORCA1 grid,  using one_deg_tripolar grid instead"
            ocean_grid = :one_deg_tripolar
        else
            @info "Using ORCA1 ocean grid"
            return CO.OceanConfigurations.orca_ocean(
                arch;
                closure,
                clock,
                depth,
                substeps,
                kwargs...,
            )
        end
    end

    @info "Using one-degree tripolar ocean grid"
    active_cells_map = !simple_ocean
    zstar = !simple_ocean

    return tripolar_ocean_simulation(
        arch;
        zstar,
        active_cells_map,
        clock,
        depth,
        Nz = 32,
        closure,
        substeps,
        kwargs...,
    )
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
- `output_dir`: Directory for output files

# Optional keyword arguments
- `simple_ocean`: Whether to use a simple ocean model setup (default: `false`)
- `ocean_grid`: Horizontal grid for Oceananigans (`:one_deg_tripolar` or `:orca`, default: `:one_deg_tripolar`)
- `depth`: Maximum ocean depth in metres (default: 5500 m for EN4 compatibility)
- `dt`: Time step (default: `1800.0` seconds, or 30 minutes)
- `comms_ctx`: Communication context (default: `ClimaComms.context()`)
- `coupled_param_dict`: Coupled parameter dictionary (default: created from `ClimaParams.create_toml_dict(FT)`)
- `ocean_diagnostic_interval`: Interval for ocean diagnostics (default: `"1days"`)
- `ocean_diagnostic_mode`: Mode for ocean diagnostics (default: `:average`)

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
    ocean_grid = :one_deg_tripolar,
    use_intersection_grid = true,
    topography_damping_factor = 5,
    depth = 5500,
    dt = 1800.0, # 30 minutes
    comms_ctx = ClimaComms.context(),
    coupled_param_dict = CP.create_toml_dict(FT),
    ocean_diagnostic_interval = "1days",
    ocean_diagnostic_mode = :average,
    extra_kwargs...,
) where {FT}
    arch = comms_ctx.device isa ClimaComms.CUDADevice ? OC.GPU() : OC.CPU()
    OC.Oceananigans.defaults.FloatType = FT

    # Use Float64 for the ocean to avoid precision issues
    FT_ocean = Float64
    OC.Oceananigans.defaults.FloatType = FT_ocean

    # Create ocean simulation
    closure = if simple_ocean
        @info "Using simpler ocean setup; to be used for software testing only."
        CO.OceanConfigurations.simplified_ocean_closure()
    else
        CO.OceanConfigurations.default_one_degree_closure()
    end

    if tspan[1] isa ITime
        # create a model clock that uses DateTime, for compatibility with ITime.
        clock = OC.TimeSteppers.Clock(time = start_date)
        stop_time = Dates.DateTime(tspan[2])
    elseif tspan[1] isa Float64
        clock = OC.TimeSteppers.Clock{Float64}(time = tspan[1])
        stop_time = tspan[2]
    else
        error("Unsupported time type: $(typeof(tspan[1]))")
    end

    ocean =
        ocean_simulation(arch, ocean_grid; clock, stop_time, simple_ocean, closure, depth)
    ocean.Δt = float(dt)

    # Set initial condition to EN4 state estimate at start_date (monthly)
    date = start_date
    # set up the `dir` keyword argument for `Metadatum`
    if date == Dates.Date(2010, 1, 1)
        # we have a ClimaArtifact saved for January 1, 2010 (so that CI can always run)
        dir_kw = (; dir = @clima_artifact("en4_temperature_salinity_2010_01"))
        @info "Using $(dir_kw.dir) ClimaArtifact for ocean initialization on $(date)"
    else
        # otherwise, download the data
        # (or load from scratchspace; ClimaOcean will automatically handle this)
        dir_kw = (;)
    end

    en4_temperature = CO.Metadatum(:temperature; date, dataset = CO.EN4Monthly(), dir_kw...)
    en4_salinity = CO.Metadatum(:salinity; date, dataset = CO.EN4Monthly(), dir_kw...)

    @info "EN4 temperature data path: $(CO.DataWrangling.metadata_path(en4_temperature))"
    @info "EN4 salinity data path: $(CO.DataWrangling.metadata_path(en4_salinity))"

    OC.set!(ocean.model, T = en4_temperature, S = en4_salinity)

    # Construct the remapper object and allocate scratch space
    grid = ocean.model.grid
    remapping = construct_remapper(
        grid,
        boundary_space;
        use_intersection_grid,
        topography_damping_factor,
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
        dt,
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
    Interfacer.progress(ocean_sim::OceananigansSimulation, cs)

Extension of `Interfacer.progress` for Oceananigans.

Print some statistics with a frequency determined by the `ocean_progress_interval` config option.
"""
function Interfacer.progress(ocean_sim::OceananigansSimulation, cs)
    ocean = ocean_sim.ocean
    model = ocean.model

    (Tmin, Tmax) = extrema(model.tracers.T)
    (Smin, Smax) = extrema(model.tracers.S)
    (ηmin, ηmax) = extrema(model.free_surface.displacement)
    umax = maximum(abs, model.velocities.u)
    vmax = maximum(abs, model.velocities.v)
    wmax = maximum(abs, model.velocities.w)
    if ClimaComms.iamroot(ClimaComms.context(cs))
        @info "Ocean | time: $(Interfacer.current_date(cs, model.clock.time)), iteration: $(OC.iteration(ocean)), " *
              "extrema(η): ($(round(ηmin, sigdigits=2)), $(round(ηmax, sigdigits=2))) " *
              "extrema(T, S): ($(round(Tmin, digits=2)), $(round(Tmax, digits=2))) ᵒC, " *
              "($(round(Smin, digits=2)), $(round(Smax, digits=2))) psu " *
              "maximum(u): ($(round(umax, sigdigits=2)), $(round(vmax, sigdigits=2)), $(round(wmax, sigdigits=2))) m/s"
    end
    return nothing
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
    construct_remapper(grid_oc, boundary_space;
                       use_intersection_grid = true,
                       topography_damping_factor = 5)

Given an Oceananigans grid and a ClimaCore boundary space, construct the two
independent sparse regridders needed to remap between them in both directions,
and (when supported) the exchange (intersection) grid used for surface
fractions and per-polygon fluxes.

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

When `use_intersection_grid = true` and the setup supports it (a
`SpectralElementSpace2D` boundary space on a single process), the returned
NamedTuple additionally carries:
- `exchange_grid`: device-resident [`ExchangeGrid`](@ref);
- `wet_ocean_fraction`: static boundary-space `CC.Field` from
  [`wet_ocean_fraction_field`](@ref), filtered consistently with the
  atmosphere's orography smoothing (`topography_damping_factor`);
- `use_exchange_grid::Bool`: whether the exchange-grid path is active.
"""
function construct_remapper(
    grid_oc,
    boundary_space;
    use_intersection_grid = true,
    topography_damping_factor = 5,
)
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

    # Exchange (intersection) grid: geometry built on CPU, applied on device.
    # Only supported for a process-local spectral-element boundary space; the
    # column (PointSpace) and distributed setups fall back to the regridder-
    # only path.
    use_exchange_grid =
        use_intersection_grid &&
        boundary_space isa CC.Spaces.SpectralElementSpace2D &&
        ClimaComms.context(boundary_space) isa ClimaComms.SingletonCommsContext
    if use_exchange_grid
        exchange_grid_cpu = build_exchange_grid(boundary_space, grid_oc)
        wet_ocean_fraction = wet_ocean_fraction_field(
            boundary_space,
            exchange_grid_cpu;
            topography_damping_factor,
        )
        exchange_grid = on_device(arch, exchange_grid_cpu)

        # Per-polygon flux scratch, boundary-space flux scratch fields (in the
        # layout `update_flux_fields!` expects), a shared DSS buffer, and the
        # DSS'd nodal wet coverage used to normalize the SE-side flux scatter.
        ocean_flux_state = ExchangeFluxState{FT}(arch, exchange_grid_cpu.n_poly)
        momentum_basis = momentum_basis_fields(boundary_space)
        ice_flux_state = IceExchangeState{FT}(arch, exchange_grid_cpu.n_poly)
        weight_cov_scratch = CC.Fields.zeros(boundary_space)
        flux_scratch = (;
            F_turb_ρτxz = CC.Fields.zeros(boundary_space),
            F_turb_ρτyz = CC.Fields.zeros(boundary_space),
            F_sh = CC.Fields.zeros(boundary_space),
            F_lh = CC.Fields.zeros(boundary_space),
            F_turb_moisture = CC.Fields.zeros(boundary_space),
        )
        flux_dss_buffer = Utilities.init_dss_buffer(flux_scratch.F_sh)
        node_cov_dss = CC.Fields.zeros(boundary_space)
        CRExt = get_ConservativeRegriddingCCExt()
        device_array_type = ClimaComms.array_type(ClimaComms.device(boundary_space))
        CRExt.vec_to_se_field!(
            node_cov_dss,
            Adapt.adapt(device_array_type, exchange_grid_cpu.node_cov),
        )
        Utilities.apply_dss!(node_cov_dss, flux_dss_buffer)
    else
        exchange_grid = nothing
        wet_ocean_fraction = nothing
        ocean_flux_state = nothing
        momentum_basis = nothing
        ice_flux_state = nothing
        weight_cov_scratch = nothing
        flux_scratch = nothing
        flux_dss_buffer = nothing
        node_cov_dss = nothing
    end

    # `TripolarGrid` covers the full sphere by construction, so no polar
    # masking is needed on either the OC or CC side. The FV → SE projection
    # has no structural-zero "no-data" nodes to repair, and `weighted_dss!`
    # operates on a fully-physical SE field.
    return (;
        remapper_oc_to_cc,
        remapper_cc_to_oc,
        scratch_field_oc1,
        scratch_field_oc2,
        scratch_field_oc3,
        temp_uv_vec,
        exchange_grid,
        wet_ocean_fraction,
        ocean_flux_state,
        momentum_basis,
        ice_flux_state,
        weight_cov_scratch,
        flux_scratch,
        flux_dss_buffer,
        node_cov_dss,
        use_exchange_grid,
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

"""
    FieldExchanger.align_surface_fractions!(ocean_sim::OceananigansSimulation,
                                            cs::Interfacer.CoupledSimulation) -> Bool

Ocean-bathymetry-authoritative surface fractions on the exchange grid.

The static wet-ocean fraction (`remapping.wet_ocean_fraction`, derived from
the intersection areas with the ocean's immersed wet mask and filtered
consistently with the atmosphere's orography smoothing) partitions each
boundary node into wet and land parts. Sea ice and open ocean subdivide the
wet part; land fills the remainder:

    ice   = clamp(ice_concentration, 0, wet)
    ocean = wet - ice
    land  = 1 - wet

so the three fractions sum to 1 identically and the flux weights are, by
construction, consistent with where the ocean model actually has wet cells
(issue #1838).

Returns `false` (falling back to the legacy ETOPO-based update) when the
exchange grid is not active.
"""
function FieldExchanger.align_surface_fractions!(
    ocean_sim::OceananigansSimulation,
    cs::Interfacer.CoupledSimulation,
)
    ocean_sim.remapping.use_exchange_grid || return false
    # Without a land model nothing can absorb the `1 - wet` remainder, so the
    # legacy residual update (ocean = 1 - ice) is the only consistent choice.
    haskey(cs.model_sims, :land_sim) || return false

    FT = CC.Spaces.undertype(Interfacer.boundary_space(cs))
    wet_fraction = ocean_sim.remapping.wet_ocean_fraction

    if haskey(cs.model_sims, :ice_sim)
        ice_sim = cs.model_sims.ice_sim
        Interfacer.get_field!(cs.fields.scalar_temp2, ice_sim, Val(:ice_concentration))
        ice_concentration = cs.fields.scalar_temp2
        @. cs.fields.scalar_temp3 = clamp(ice_concentration, FT(0), wet_fraction)
        ice_fraction = cs.fields.scalar_temp3
        Interfacer.update_field!(ice_sim, Val(:area_fraction), ice_fraction)
    else
        cs.fields.scalar_temp3 .= FT(0)
        ice_fraction = cs.fields.scalar_temp3
    end

    @. cs.fields.scalar_temp2 = max(wet_fraction - ice_fraction, FT(0))
    ocean_fraction = cs.fields.scalar_temp2
    Interfacer.update_field!(ocean_sim, Val(:area_fraction), ocean_fraction)

    @. cs.fields.scalar_temp1 = max(FT(1) - wet_fraction, FT(0))
    land_fraction = cs.fields.scalar_temp1
    Interfacer.update_field!(cs.model_sims.land_sim, Val(:area_fraction), land_fraction)
    cs.fields.land_area_fraction .= land_fraction

    if haskey(cs.model_sims, :ice_sim)
        FieldExchanger.resolve_area_fractions!(
            ocean_sim,
            cs.model_sims.ice_sim,
            land_fraction,
        )
    end

    cs.fields.ice_area_fraction .= ice_fraction
    cs.fields.ocean_area_fraction .= ocean_fraction

    @assert minimum(ice_fraction .+ land_fraction .+ ocean_fraction) ≈ FT(1)
    @assert maximum(ice_fraction .+ land_fraction .+ ocean_fraction) ≈ FT(1)
    return true
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
    if sim.remapping.use_exchange_grid
        # The exchange-grid path pushes the per-polygon fluxes currently held
        # in `sim.remapping.ocean_flux_state` (`fields` is ignored: the
        # boundary-space fields are a coarser view of the same fluxes).
        push_exchange_fluxes_to_ocean!(sim)
    else
        _update_turbulent_fluxes_boundary!(sim, fields)
    end
    return nothing
end

function _update_turbulent_fluxes_boundary!(sim::OceananigansSimulation, fields)
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
    OC.interior(oc_flux_S, :, :, 1) .-=
        (1.0 .- ice_concentration) .* surface_salinity .* moisture_fresh_water_flux
    return nothing
end

"""
    FluxCalculator.compute_surface_fluxes!(csf, sim::OceananigansSimulation,
                                           atmos_sim, thermo_params,
                                           accumulator = nothing)

Compute atmosphere-ocean turbulent fluxes on the exchange (intersection)
grid: one SurfaceFluxes evaluation per polygon, with the atmospheric state
gathered from the SE nodes and the SST read directly from the owning ocean
cell. The per-polygon fluxes are aggregated to the boundary space (see
`scatter_poly_fluxes_to_boundary!`) and handed to
`FluxCalculator.update_flux_fields!`, which applies the (exchange-grid
derived) area-fraction weighting for the coupler sums and triggers the
ocean-side push or the accumulator.

Falls back to the generic boundary-space computation when the exchange grid
is not active.
"""
NVTX.@annotate function FluxCalculator.compute_surface_fluxes!(
    csf,
    sim::OceananigansSimulation,
    atmos_sim::Interfacer.AbstractAtmosSimulation,
    thermo_params,
    accumulator = nothing,
)
    remapping = sim.remapping
    if !remapping.use_exchange_grid
        return invoke(
            FluxCalculator.compute_surface_fluxes!,
            Tuple{
                Any,
                Interfacer.AbstractSurfaceSimulation,
                Interfacer.AbstractAtmosSimulation,
                Any,
                Any,
            },
            csf,
            sim,
            atmos_sim,
            thermo_params,
            accumulator,
        )
    end

    FT = CC.Spaces.undertype(axes(csf))
    eg = remapping.exchange_grid
    fs = remapping.ocean_flux_state
    surface_fluxes_params = FluxCalculator.get_surface_params(atmos_sim)

    # Gather the atmospheric and ocean-surface state onto the polygons.
    gather_atmos_state_to_polys!(fs, eg, csf, remapping.temp_uv_vec, remapping.momentum_basis)
    Nz = size(sim.ocean.model.grid, 3)
    gather_cells_to_polys!(
        fs.T_sfc,
        eg,
        vec(OC.interior(sim.ocean.model.tracers.T, :, :, Nz)),
    )
    fs.T_sfc .+= FT(sim.ocean_properties.C_to_K)
    gather_cells_to_polys!(fs.sic, eg, vec(OC.interior(sim.ice_concentration, :, :, 1)))

    # COARE3 roughness is spatially uniform, so the whole flux configuration
    # is a scalar kernel argument.
    config = SF.SurfaceFluxConfig(
        SF.COARE3RoughnessParams{FT}(),
        SF.ConstantGustinessSpec(FT(1)),
    )
    compute_ocean_polygon_fluxes!(fs, surface_fluxes_params, thermo_params, config)
    fs.n_acc[] += 1

    # The ocean fluxes apply to the open-water part of each polygon.
    @. fs.scratch2 = 1 - fs.sic
    scatter_poly_fluxes_to_boundary!(remapping, eg, fs, fs.scratch2)
    FluxCalculator.update_flux_fields!(csf, sim, remapping.flux_scratch, accumulator)
    return nothing
end

"""
    push_exchange_fluxes_to_ocean!(sim::OceananigansSimulation)

Push the per-polygon turbulent fluxes currently held in
`sim.remapping.ocean_flux_state` into the ocean boundary conditions, with the
sea-ice-concentration weighting applied *per polygon* (rather than per ocean
cell after remapping, as in the boundary-space path). Sign and scaling
conventions match `_update_turbulent_fluxes_boundary!`: momentum sets the
velocity flux BCs, heat and salinity accumulate onto the tracer flux BCs.
"""
NVTX.@annotate function push_exchange_fluxes_to_ocean!(sim::OceananigansSimulation)
    remapping = sim.remapping
    eg = remapping.exchange_grid
    fs = remapping.ocean_flux_state
    (; reference_density, heat_capacity) = sim.ocean_properties
    grid = sim.ocean.model.grid

    # Momentum (UV basis): weight by open-ocean fraction per polygon, scatter
    # to cells, mirror the tripolar fold, then rotate/stagger onto the C-grid.
    @. fs.scratch1 = fs.F_τu * (1 - fs.sic) / reference_density
    @. fs.scratch2 = fs.F_τv * (1 - fs.sic) / reference_density
    τu_cells = vec(OC.interior(remapping.scratch_field_oc1, :, :, 1))
    τv_cells = vec(OC.interior(remapping.scratch_field_oc2, :, :, 1))
    scatter_polys_to_cells!(τu_cells, eg, fs.scratch1)
    scatter_polys_to_cells!(τv_cells, eg, fs.scratch2)
    mirror_fold_partners!(τu_cells, grid)
    mirror_fold_partners!(τv_cells, grid)
    oc_flux_u = surface_flux(sim.ocean.model.velocities.u)
    oc_flux_v = surface_flux(sim.ocean.model.velocities.v)
    set_from_extrinsic_vector!(
        (; u = oc_flux_u, v = oc_flux_v),
        grid,
        remapping.scratch_field_oc1,
        remapping.scratch_field_oc2,
    )

    # Heat: (1 - SIC)-weighted turbulent heat flux per polygon.
    @. fs.scratch1 =
        (1 - fs.sic) * (fs.F_lh + fs.F_sh) / (reference_density * heat_capacity)
    heat_cells = vec(OC.interior(remapping.scratch_field_oc3, :, :, 1))
    scatter_polys_to_cells!(heat_cells, eg, fs.scratch1)
    mirror_fold_partners!(heat_cells, grid)
    oc_flux_T = surface_flux(sim.ocean.model.tracers.T)
    OC.interior(oc_flux_T, :, :, 1) .+= OC.interior(remapping.scratch_field_oc3, :, :, 1)

    # Salinity: moisture flux (upward positive) per polygon; multiplied by the
    # local surface salinity at the cell level.
    @. fs.scratch1 = (1 - fs.sic) * fs.F_moisture / reference_density
    moisture_cells = vec(OC.interior(remapping.scratch_field_oc3, :, :, 1))
    scatter_polys_to_cells!(moisture_cells, eg, fs.scratch1)
    mirror_fold_partners!(moisture_cells, grid)
    oc_flux_S = surface_flux(sim.ocean.model.tracers.S)
    surface_salinity = OC.interior(sim.ocean.model.tracers.S, :, :, grid.Nz)
    OC.interior(oc_flux_S, :, :, 1) .-=
        surface_salinity .* OC.interior(remapping.scratch_field_oc3, :, :, 1)
    return nothing
end

"""
    FluxCalculator.push_and_reset!(sim::OceananigansSimulation, acc)

Slow-surface push for the exchange-grid path: time-average the *per-polygon*
flux accumulators into the flux-state outputs, push them to the ocean
boundary conditions, and reset both the per-polygon and the boundary-space
accumulators. Falls back to the generic boundary-space behavior when the
exchange grid is not active.
"""
function FluxCalculator.push_and_reset!(
    sim::OceananigansSimulation,
    acc::FluxCalculator.FluxAccumulator,
)
    if !sim.remapping.use_exchange_grid
        return invoke(
            FluxCalculator.push_and_reset!,
            Tuple{Any, FluxCalculator.FluxAccumulator},
            sim,
            acc,
        )
    end
    fs = sim.remapping.ocean_flux_state
    average_and_reset_exchange_accumulators!(fs) || return nothing
    push_exchange_fluxes_to_ocean!(sim)
    # Keep the (unused but still accumulated) boundary-space accumulator in
    # sync so its averages stay well-defined.
    FluxCalculator.reset!(acc)
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

    # Compute radiative contribution.
    rad_T_flux = OC.interior(sim.remapping.scratch_field_oc1, :, :, 1)
    rad_T_flux .=
        (1.0 .- ice_concentration) .* (
            -(1 - α) .* remapped_SW_d .-
            ϵ * (
                remapped_LW_d .-
                σ .* (C_to_K .+ OC.interior(sim.ocean.model.tracers.T, :, :, Nz)) .^ 4
            )
        ) ./ (reference_density * heat_capacity)
    OC.interior(oc_flux_T, :, :, 1) .+= rad_T_flux

    # Remap precipitation fields onto scratch arrays; rename for clarity
    Interfacer.remap!(sim.remapping.scratch_field_oc1, csf.P_liq, sim.remapping) # liquid precipitation
    remapped_P_liq = OC.interior(sim.remapping.scratch_field_oc1, :, :, 1)

    Interfacer.remap!(sim.remapping.scratch_field_oc2, csf.P_snow, sim.remapping) # snow precipitation
    remapped_P_snow = OC.interior(sim.remapping.scratch_field_oc2, :, :, 1)

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
