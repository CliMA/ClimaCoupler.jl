#
# driverdriver.jl::
#
# trivial coupled system mock up
#
# very simple top level driver that show proof of concept implementation for 
# simple couplng API sufficient for v0 atmos-land-ocean DG tools.
# Components
# _A - atmos
# _L - land
# _O - ocean
# _C - coupler
#
# driverdriver multiplexes land ocean atmos and coupler components.
#
# all exist on aligned horizontal grids decomposed on ranks in aligned way.
#
# components in this example run sequentially on shared set/subset of 
# underlying resources without spinwait contention.
#
# components have their own MPI communicators to support MPI rank disaggregation onto
# distinct resources or to interleave on same resources with some lightweight
# resource contention coordination.
#
# disaggregated ranks require coupler to handle asynchronous send/receive posting
# 
#
abstract type CLIMAcomponent end
abstract type  DGModelComponent <: CLIMAcomponent end
struct CouplerComponent         <: CLIMAcomponent end
struct AtmosModelComponent      <: DGModelComponent end
struct LandModelComponent       <: DGModelComponent end
struct OceanModelComponent      <: DGModelComponent end

comp_config(args...; kargs...) = nothing;
#
using MPI
using ClimateMachine
using ClimateMachine.Mesh.Topologies
using ClimateMachine.Mesh.Grids
using ClimateMachine.DGMethods
using ClimateMachine.DGMethods.NumericalFluxes

 #####
 # Basic initialization
 #####
 ClimateMachine.init()
 ArrayType = ClimateMachine.array_type()
 mpicomm = MPI.COMM_WORLD
 FT = Float64

 # In general we might dup communicator and/or create intercommunicators/
 # intracommicators depending on resource use and parallelism mode. Here we 
 # assume sequential, all componenta on all same ranks, aligned. So each component 
 # just uses default mpicomm - (also not sure if CLiMA as written always uses 
 # component communicator). We carry different names so we can make them different 
 # in the future.
 mpicomm_A=mpicomm
 mpicomm_O=mpicomm
 mpicomm_L=mpicomm
 mpicomm_C=mpicomm

 #####
 # Create mesh
 #####

 # Horiz mesh
 # Same horizontal mesh size, extents and polynomical order for all components
 # We still create distinct grids so that they can use different communicators when
 # needed.
 const N = 4
 const Nˣ = 10
 const Nʸ = 10
 const Lˣ = 4e6  # m
 const Lʸ = 4e6  # m
 xrange = range(FT(0); length = Nˣ + 1, stop = Lˣ)
 yrange = range(FT(0); length = Nʸ + 1, stop = Lʸ)
 brickrange_2D = (xrange, yrange)
 topol_2D( comm ) = (
   BrickTopology( comm, brickrange_2D, periodicity = (false, false) ) );
 grid_2D( topol ) = ( DiscontinuousSpectralElementGrid(
   topol,
   FloatType = FT,
   DeviceArray = ArrayType,
   polynomialorder = N,
 )
 )
 topo_2D_A  = topol_2D( mpicomm_A )
 grid_2D_A  = grid_2D(  topo_2D_A )
 topo_2D_O  = topol_2D( mpicomm_O )
 grid_2D_O  = grid_2D(  topo_2D_O )
 topo_2D_L  = topol_2D( mpicomm_L )
 grid_2D_L  = grid_2D(  topo_2D_L )
 topo_2D_C  = topol_2D( mpicomm_C )
 grid_2D_C  = grid_2D(  topo_2D_C )


 # Vert mesh (stacked topology)
 # Different components can have different vertical extents and element counts
 const Nᶻ = 20
 # Atmos vert
 const Nᶻ_A = 20
 # Ocean vert
 const Nᶻ_O = 50
 # Land vert - lets makes this 2d as an example
 const Nᶻ_L =  0
 # Coupler vert - lets makes this 2d too for now
 const Nᶻ_C =  0

 const H_A = 100000  # m
 const H_O = 5000    # m
 zrange_A  = range(FT(0   ); length = Nᶻ_A + 1, stop = H_A)
 zrange_O  = range(FT(-H_O); length = Nᶻ_O + 1, stop = 0  )
 brickrange_3D_A = (xrange, yrange, zrange_A)
 brickrange_3D_O = (xrange, yrange, zrange_O)
 topol_3D( comm, br ) = (
  StackedBrickTopology(
    mpicomm,
    br;
    periodicity = (false, false, false),
    boundary = ((1, 1), (1, 1), (2, 3)),
  )
 );
 topol_3D_A = topol_3D( mpicomm_A, brickrange_3D_A  )
 topol_3D_O = topol_3D( mpicomm_A, brickrange_3D_O  )
 grid_3D( topol ) = ( 
  DiscontinuousSpectralElementGrid(
   topol,
   FloatType = FT,
   DeviceArray = ArrayType,
   polynomialorder = N,
  )
 );
 grid_3D_A=grid_3D( topol_3D_A )
 grid_3D_O=grid_3D( topol_3D_O )

 # Config and init coupler
 comp_C=comp_config( CouplerComponent(), ( (), ) )

 # Config (including mask) and init drivers of individual components
 comp_A=comp_config( AtmosModelComponent() ; grids=( (grid=grid_3D_A, label="g3da",), ) )
 comp_L=comp_config( LandModelComponent()  ; grid2D=grid_2D_L )
 comp_O=comp_config( OceanModelComponent() ; grid3D=grid_3D_O )
