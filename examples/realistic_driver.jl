# Import necessary packages
import ClimaComms
import ClimaCore as CC
import ClimaAtmos as CA
import ClimaLand as CL

# Include adapters for component models
## helpers for component models
include("../experiments/ClimaEarth/components/atmosphere/climaatmos.jl")
include("../experiments/ClimaEarth/components/land/climaland_bucket.jl")
include("../experiments/ClimaEarth/components/ocean/slab_ocean.jl")
include("../experiments/ClimaEarth/components/ocean/prescr_seaice.jl")

# Choose device to run on (CPU or GPU)
context = ClimaComms.context()
device = ClimaComms.device(context)

# Choose float type
FT = Float64

# Set up simulation space
space = create_sphere(FT, context)

# Set up timestepping information
tspan = (; t_start = Float64(0), t_end = Float64(60 * 60 * 24)) # simulation length 1 day
dt_cpl = (400)

# Initialize component models
atmos = atmos_init()
land = land_init() # or stub_init()
ocean = ocean_init() # or stub_init()
ice = ice_init() # or stub_init()

# Initialize coupled simulation
coupled_simulation = CoupledSimulation(atmos, ocean, land, seaice, tspan, dt_cpl)

# Run simulation
solve_coupler!(coupled_simulation)

# Postprocessing
postprocess(coupled_simulation)


"""
    create_sphere(FT, context)

Create a spherical shell space. Use a default radius, depth, and number
of horizontal elements. For now we just create a spherical shell with depth,
representing the subsurface, and take the top of the shell to be the surface.
We should figure out how to construct subsurface, surface, and above-surface
spaces if we want all component models to inherit all 3 from the coupler.
"""
function create_sphere(FT, context)
    # # TODO how can we specify number of atmos vert elements and number of land vert elements separately?
    nelements = (16, 15) # 16 horizontal elements; 15 vertical elements
    # h_elem = 4 # number of horizontal elements

    radius = FT(6378.1e3)
    depth = FT(50)

    vertdomain = ClimaCore.Domains.IntervalDomain(
        ClimaCore.Geometry.ZPoint(FT(-depth)),
        ClimaCore.Geometry.ZPoint(FT(0));
        boundary_names = (:bottom, :top),
    )
    if dz_tuple isa Nothing
        vertmesh =
            ClimaCore.Meshes.IntervalMesh(vertdomain; nelems = nelements[2])
    else
        vertmesh = ClimaCore.Meshes.IntervalMesh(
            vertdomain,
            ClimaCore.Meshes.GeneralizedExponentialStretching{FT}(
                dz_tuple[1],
                dz_tuple[2],
            );
            nelems = nelements[2],
            reverse_mode = true,
        )
    end
    device = ClimaComms.device(context)
    if pkgversion(ClimaCore) >= v"0.14.10"
        vert_center_space =
            ClimaCore.Spaces.CenterFiniteDifferenceSpace(device, vertmesh)
    else
        vert_center_space =
            ClimaCore.Spaces.CenterFiniteDifferenceSpace(vertmesh)
    end

    horzdomain = ClimaCore.Domains.SphereDomain(radius)
    horzmesh = ClimaCore.Meshes.EquiangularCubedSphere(horzdomain, nelements[1])
    horztopology = ClimaCore.Topologies.Topology2D(comms_ctx, horzmesh)
    quad = ClimaCore.Spaces.Quadratures.GLL{npolynomial + 1}()
    horzspace = ClimaCore.Spaces.SpectralElementSpace2D(horztopology, quad)

    subsurface_space = ClimaCore.Spaces.ExtrudedFiniteDifferenceSpace(
        horzspace,
        vert_center_space,
    )
    boundary_space = CC.Spaces.horizontal_space(face_space)
    space = (; surface = boundary_space, subsurface = subsurface_space)
    return boundary_space
end
