"""
    Utilities

This module contains functions, objects, and constants used by various
modules in the coupler.
"""
module Utilities

import ClimaComms
import ClimaCore as CC
import Logging
import ClimaUtilities.OutputPathGenerator: generate_output_path
import ..Interfacer

export get_device,
    get_comms_context,
    show_memory_usage,
    setup_output_dirs,
    time_to_seconds,
    integral,
    create_boundary_space,
    diagnostics_global_attribs,
    report_nans

"""
    get_device(config_dict)

Returns the device on which the model is being run

# Arguments
- `config_dict`: dictionary containing a "device" flag which decides which device to run on
"""
function get_device(config_dict)
    if config_dict["device"] == "auto"
        return ClimaComms.device()
    elseif config_dict["device"] == "CUDADevice"
        return ClimaComms.CUDADevice()
    elseif config_dict["device"] == "CPUMultiThreaded" || Threads.nthreads() > 1
        return ClimaComms.CPUMultiThreaded()
    else
        return ClimaComms.CPUSingleThreaded()
    end
end


"""
    get_comms_context(config_dict)

Sets up the appropriate ClimaComms context for the device the model is to be run on,
choosing from the following options:
    - CPU single threaded
    - CPU with MPI
    - GPU

If no device is passed to `ClimaComms.context()` then `ClimaComms` automatically
selects the device from which this code is called.

# Arguments
`config_dict`: dictionary containing a "device" flag which decides which device context is needed
"""
function get_comms_context(config_dict)
    device = get_device(config_dict)
    comms_ctx = ClimaComms.context(device)
    ClimaComms.init(comms_ctx)

    if pkgversion(ClimaComms) < v"0.6.6"
        # For older versions of ClimaComms, we have to manually ensure that logging only
        # happens on the root process when using multiple processes
        ClimaComms.iamroot(comms_ctx) || Logging.disable_logging(Logging.AboveMaxLevel)
    end

    if comms_ctx isa ClimaComms.SingletonCommsContext
        @info "Setting up single-process ClimaCoupler run on device: $(nameof(typeof(device)))."
    else
        @info "Setting up distributed ClimaCoupler run on " nprocs =
            ClimaComms.nprocs(comms_ctx) device = "$(nameof(typeof(device)))"
    end

    return comms_ctx
end

"""
    show_memory_usage()

Display and return the maximum resident set size (RSS) memory footprint on the
CPU of this process since it began.
"""
function show_memory_usage()
    cpu_max_rss_GB = ""
    cpu_max_rss_GB = string(round(Sys.maxrss() / 1e9, digits = 3)) * " GiB"
    @info "Memory in use: $(cpu_max_rss_GB)"
    return cpu_max_rss_GB
end

"""
    setup_output_dirs(output_dir_root = pwd(),
        artifacts_dir = joinpath(output_dir, "artifacts"),
        checkpoints_dir = joinpath(output_dir, "checkpoints"),
        regrid_dir = nothing,
        comms_ctx,
    )

Create output directories for the experiment. If `comms_ctx` is provided,
only the root process will create the directories.
By default, the artifacts and checkpoints directories are created inside the
root output directory with the names `artifacts/` and `checkpoints/`.
The regrid directory is by default created as a temporary directory inside the
root output directory and is automatically deleted when the process exits.

`ClimaUtilities.OutputPathGenerator` is used so that simulations can be re-run and re-started.
The output path looks like:
```
coupler_output_dir_amip/
├── checkpoints
│       └── checkpoints for the various models
├── artifacts
│       └── plots produced by the postporcessing step
├── regrid_tmp_<random_tempdir>/
│       └── temporary files used for regridding
├── output_0000/
│   ├── atmos/
│   │   └── output of the atmos model
│   └── ocean/
│       └── output of the ocean model
├── output_0001/
│   └── ... component model outputs in their folders ...
├── output_0002/
│   └── ... component model outputs in their folders ...
└── output_active -> output_0002/
```

# Arguments
- `output_dir::String`: The directory where the output files will be stored.
        Default is the current directory.
- `artifacts_dir::String`: The directory where plots (from postprocessing and conservation checks) will be stored.
        Default is `output_dir/artifacts/`.
- `checkpoint_dir::String`: The directory where the checkpoint files will be stored.
        Default is `output_dir/checkpoints/`.
- `regrid_dir::String`: The directory where the regridded files will be stored.
        Default is `output_dir/regrid_tmp_<random_tempdir>/`.
- `comms_ctx::Union{Nothing, ClimaComms.AbstractCommsContext}`: The communicator context.
        If provided, only the root process will create the directories.

# Returns
- A NamedTuple with the paths to output directories for each component and the coupler,
as well as paths toartifacts, regrid, and checkpoints directories.
"""
function setup_output_dirs(;
    output_dir_root = pwd(),
    artifacts_dir = joinpath(output_dir_root, "artifacts"),
    checkpoints_dir = joinpath(output_dir_root, "checkpoints"),
    regrid_dir = nothing,
    comms_ctx,
)
    output_dir_root = generate_output_path(output_dir_root, context = comms_ctx)

    # Make component-specific output directories, and one for the coupler
    atmos_output_dir = joinpath(output_dir_root, "clima_atmos")
    land_output_dir = joinpath(output_dir_root, "clima_land")
    ocean_output_dir = joinpath(output_dir_root, "clima_ocean")
    ice_output_dir = joinpath(output_dir_root, "clima_seaice")
    coupler_output_dir = joinpath(output_dir_root, "clima_coupler")

    if ClimaComms.iamroot(comms_ctx)
        mkpath(atmos_output_dir)
        mkpath(land_output_dir)
        mkpath(ocean_output_dir)
        mkpath(ice_output_dir)
        mkpath(coupler_output_dir)

        mkpath(artifacts_dir)
        mkpath(checkpoints_dir)
        # If no regrid_dir is provided, create a temporary directory
        # Note this must be done on the root process because `mktempdir` chooses a random name
        regrid_dir =
            isnothing(regrid_dir) ? mktempdir(; prefix = "regrid_tmp_") : mkpath(regrid_dir)
    end
    regrid_dir = ClimaComms.bcast(comms_ctx, regrid_dir)

    @info "Output directories $(output_dir_root)"
    @info "Coupler artifacts directory $(artifacts_dir)"
    @info "Coupler checkpoint directory $(checkpoints_dir)"
    return (;
        output_dir_root,
        atmos_output_dir,
        land_output_dir,
        ocean_output_dir,
        ice_output_dir,
        coupler_output_dir,
        artifacts_dir,
        regrid_dir,
        checkpoints_dir,
    )
end

"""
    time_to_seconds(s::String)

Convert a string to seconds. The string should be in the format `numberunit`, where `unit` is one of `secs`, `mins`, `hours`, or `days`.

# Arguments
- `s::String`: The string to convert to seconds.

# Returns
- The number of seconds represented by the string.
"""
function time_to_seconds(s::String)
    factor = Dict("secs" => 1, "mins" => 60, "hours" => 60 * 60, "days" => 60 * 60 * 24)
    s == "Inf" && return Inf
    if count(occursin.(keys(factor), Ref(s))) != 1
        error("Bad format for flag $s. Examples: [`10secs`, `20mins`, `30hours`, `40days`]")
    end
    for match in keys(factor)
        occursin(match, s) || continue
        return parse(Float64, first(split(s, match))) * factor[match]
    end
    error("Uncaught case in computing time from given string.")
end

"""
    diagnostics_global_attribs(start_date)

Build a dictionary of global (file-level) attributes to attach to the diagnostic
NetCDF files produced by the coupler. These are passed to the `global_attribs`
keyword of `ClimaDiagnostics.Writers.NetCDFWriter` and stored as metadata in
every NetCDF file, which is useful when visualizing or sharing output.

Always includes the simulation `start_date`. When running on Buildkite, also
includes the build number (the "model run number") read from the
`BUILDKITE_BUILD_NUMBER` environment variable; this is omitted on local runs
where that variable is not set.

Also includes the git commit hash of the code being run, read from the
`BUILDKITE_COMMIT` environment variable on Buildkite or, on local runs, from
`git rev-parse HEAD`. It is omitted if neither source is available (e.g. when
the working directory is not a git repository).

The returned dictionary maps `String` keys to `String` values, as required by
`NetCDFWriter`.
"""
function diagnostics_global_attribs(start_date)
    attribs = Dict{String, String}("start_date" => string(start_date))
    build_number = get(ENV, "BUILDKITE_BUILD_NUMBER", "")
    isempty(build_number) || (attribs["buildkite_build_number"] = build_number)
    commit_hash = get(ENV, "BUILDKITE_COMMIT", "")
    if isempty(commit_hash)
        # Fall back to querying git directly on local runs
        commit_hash = try
            readchomp(`git rev-parse HEAD`)
        catch
            ""
        end
    end
    isempty(commit_hash) || (attribs["commit_hash"] = commit_hash)
    return attribs
end

"""
    integral(field)

Return the integral (a scalar) for `field` along its spatial dimensions.
"""
function integral(field::CC.Fields.Field)
    if axes(field) isa CC.Spaces.SpectralElementSpace2D
        # ClimaCore #1578, if field comes from Fields.level
        # TODO: This should be fixed in ClimaCore
        rad_grid = CC.Spaces.grid(axes(field))
        if rad_grid isa CC.Grids.LevelGrid
            if rad_grid.level isa CC.Utilities.PlusHalf
                # FaceSpace
                return sum(field ./ CC.Fields.Δz_field(field) .* 2)
            else
                # Center
                return sum(field ./ CC.Fields.Δz_field(field))
            end
        else
            return sum(field)
        end
    else
        return sum(field)
    end
end

function integral(field::AbstractArray)
    return sum(field)
end

"""
    _column_boundary_space(::Type{FT}, latlon, comms_ctx) where {FT}

Create a `PointSpace` with lat/long coordinates for use as a single-column
boundary space. Constructs a minimal extruded domain with `LatPoint`/`LongPoint`
intervals centered at the given coordinates, then extracts the top face level
as a `PointSpace`.

The resulting `PointSpace` has `LatLongZPoint` coordinates (with `.lat`, `.long`,
and `.z`), so `SpaceVaryingInput`, `TimeVaryingInput`, insolation calculations,
and coordinate-dependent physics all work correctly.

# Arguments
- `FT`: floating-point type
- `latlon`: tuple of `(latitude, longitude)` in degrees
- `comms_ctx`: ClimaComms context
"""
function _column_boundary_space(::Type{FT}, latlon, comms_ctx) where {FT}
    lat, long = FT.(latlon)
    device = ClimaComms.device(comms_ctx)

    # Build a small non-degenerate horizontal domain centered at (lat, long).
    # GL{1} gives 1 quadrature point per element, located at the interval center.
    offset = FT(0.2)
    domain_x = CC.Domains.IntervalDomain(
        CC.Geometry.LatPoint(lat - offset),
        CC.Geometry.LatPoint(lat + offset);
        boundary_names = (:north, :south),
    )
    domain_y = CC.Domains.IntervalDomain(
        CC.Geometry.LongPoint(long - offset),
        CC.Geometry.LongPoint(long + offset);
        boundary_names = (:west, :east),
    )
    plane = CC.Domains.RectangleDomain(domain_x, domain_y)
    mesh = CC.Meshes.RectilinearMesh(plane, 1, 1)
    topology = CC.Topologies.Topology2D(comms_ctx, mesh)
    quad = CC.Spaces.Quadratures.GL{1}()
    hspace = CC.Spaces.SpectralElementSpace2D(topology, quad)

    # Minimal vertical domain
    z_sfc = FT(0)
    vertdomain = CC.Domains.IntervalDomain(
        CC.Geometry.ZPoint(z_sfc - FT(1)),
        CC.Geometry.ZPoint(z_sfc);
        boundary_names = (:bottom, :top),
    )
    vertmesh = CC.Meshes.IntervalMesh(vertdomain; nelems = 2)
    vert_center_space = CC.Spaces.CenterFiniteDifferenceSpace(device, vertmesh)
    center_extruded = CC.Spaces.ExtrudedFiniteDifferenceSpace(hspace, vert_center_space)

    # Extract the top face level as a PointSpace.
    # level(FaceFiniteDifferenceSpace, PlusHalf) → PointSpace
    face_extruded = CC.Spaces.face_space(center_extruded)
    face_column = CC.Spaces.column(face_extruded, 1, 1, 1)
    point_space = CC.Spaces.level(face_column, CC.Utilities.half)

    return point_space
end

"""
    create_boundary_space(::Type{FT}, domain_type, atmos_sim, share_surface_space, comms_ctx; kwargs...)

Construct the 2D boundary space used for coupler field exchange.

For `domain_type == "column"`, returns a `PointSpace` with lat/long coordinates.
For global simulations, returns either the atmosphere's horizontal surface space
(when `share_surface_space` is true) or an independent `CubedSphereSpace`.

# Arguments
- `FT`: floating-point type
- `domain_type`: `"global"` or `"column"`
- `atmos_sim`: atmosphere simulation (used when sharing surface space)
- `share_surface_space`: whether to reuse the atmosphere's horizontal space
- `comms_ctx`: ClimaComms context
- `column_latlon`: `(lat, lon)` tuple, required when `domain_type == "column"`
- `nh_poly`: polynomial order, required when not sharing surface space
- `h_elem`: number of horizontal elements, required when not sharing surface space
- `coupled_param_dict`: parameter dictionary, required when not sharing surface space
"""
function create_boundary_space(
    ::Type{FT},
    domain_type,
    atmos_sim,
    share_surface_space,
    comms_ctx;
    column_latlon = nothing,
    nh_poly_coupler = nothing,
    h_elem_coupler = nothing,
    coupled_param_dict = nothing,
) where {FT}
    if domain_type == "column"
        return _column_boundary_space(FT, column_latlon, comms_ctx)
    elseif share_surface_space
        return CC.Spaces.horizontal_space(atmos_sim.domain.face_space)
    else
        n_quad_points = nh_poly_coupler + 1
        radius = coupled_param_dict["planet_radius"]
        return CC.CommonSpaces.CubedSphereSpace(
            FT;
            radius,
            n_quad_points,
            h_elem = h_elem_coupler,
        )
    end
end

"""
    init_dss_buffer(field::CC.Fields.Field)
    init_dss_buffer(field_object::Union{CC.Fields.FieldVector, NamedTuple})

Initialize a DSS buffer for the given Field or FieldVector object. DSS only
applies to spectral-element (horizontal) spaces. If the object is defined on a
space with no horizontal spectral-element direction (a `PointSpace` or a column
`FiniteDifferenceSpace`), return `nothing`.
"""
function init_dss_buffer(field::CC.Fields.Field)
    space = axes(field)
    if space isa CC.Spaces.PointSpace || space isa CC.Spaces.AbstractFiniteDifferenceSpace
        return nothing
    else
        return CC.Spaces.create_dss_buffer(field)
    end
end
function init_dss_buffer(field_object::Union{CC.Fields.FieldVector, NamedTuple})
    buffers = map(
        key -> init_dss_buffer(getproperty(field_object, key)),
        propertynames(field_object),
    )
    all(isnothing, buffers) && return nothing
    return NamedTuple{propertynames(field_object)}(buffers)
end

"""
    apply_dss!(field::Union{CC.Fields.Field, CC.Fields.FieldVector}, buffer)

Apply DSS to the given Field or FieldVector object.
If the object is defined on a `PointSpace`, return `nothing`.
"""
function apply_dss!(field::Union{CC.Fields.Field, CC.Fields.FieldVector}, buffer)
    if isnothing(buffer)
        return nothing
    else
        return CC.Spaces.weighted_dss!(field, buffer)
    end
end

"""
    report_nans(obj; name = "state") → Int

Scan every scalar leaf field of `obj` (a `ClimaCore.Fields.Field` or
`FieldVector`, e.g. an atmos integrator state) for NaNs. For each contaminated
leaf, log an error with the NaN count and the coordinate bounding box
(lat/long/z, when the space has them) of the NaN locations. Return the total
NaN count, so a caller can halt with `report_nans(obj) == 0 || error(...)`.

Unlike the ClimaAtmos `check_nans` callback (a bare `any(isnan, parent(u))`),
this localizes which variable went NaN and roughly where, e.g.

    NaNs in atmos.c.ρe_tot: 12 of 55296 values, lat ∈ [67.3, 71.2], long ∈ [-162.4, -158.0], z ∈ [30.0, 30.0]

Uses only broadcasts and array reductions, so it is GPU-safe; reductions are
per-process (no MPI reduction). Intended for debugging, not hot loops: the
bounding box allocates a few temporary fields per contaminated leaf.
"""
function report_nans(obj; name = "state")
    total = 0
    for chain in CC.Fields.property_chains(obj)
        # `identity`: index covariant components directly instead of
        # converting whole vector fields to UVVector/WVector per leaf
        leaf = CC.Fields.single_field(obj, chain, identity)
        eltype(leaf) <: AbstractFloat || continue
        n = count(isnan, parent(leaf))
        n == 0 && continue
        total += n
        label = string(name, join(map(p -> p isa Symbol ? ".$p" : "[$p]", chain)))
        @error "NaNs in $label: $n of $(length(parent(leaf))) values$(nan_bounding_box(leaf))"
    end
    return total
end

"""
    report_nans(cs::Interfacer.CoupledSimulation) → Int

Scan the coupler exchange fields and the prognostic state of every component
simulation with a ClimaCore state (`sim.integrator.u`, i.e. atmos and land;
Oceananigans-based components are skipped) and log the location of any NaNs.
Return the total NaN count. Typical use in a manual stepping loop:

    ClimaCoupler.step!(cs)
    Utilities.report_nans(cs) == 0 || error("NaN at \$(Interfacer.current_date(cs))")
"""
function report_nans(cs::Interfacer.CoupledSimulation)
    total = report_nans(cs.fields; name = "coupler_fields")
    for (sim_name, sim) in pairs(cs.model_sims)
        hasproperty(sim, :integrator) || continue
        u = sim.integrator.u
        u isa Union{CC.Fields.Field, CC.Fields.FieldVector} || continue
        total += report_nans(u; name = string(sim_name))
    end
    return total
end

# Coordinate bounding box of the NaN entries of scalar field `leaf`, as a
# string; empty for spaces without lat/long/z coordinates.
function nan_bounding_box(leaf)
    coords = CC.Fields.coordinate_field(axes(leaf))
    mask = isnan.(leaf)
    io = IOBuffer()
    for p in (:lat, :long, :z)
        p in propertynames(coords) || continue
        c = getproperty(coords, p)
        FT = eltype(parent(c))
        lo = minimum(parent(ifelse.(mask, c, FT(Inf))))
        hi = maximum(parent(ifelse.(mask, c, FT(-Inf))))
        print(io, ", $p ∈ [$(round(lo; sigdigits = 6)), $(round(hi; sigdigits = 6))]")
    end
    return String(take!(io))
end

end # module
