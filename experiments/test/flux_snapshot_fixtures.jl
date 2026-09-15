#=
Fixtures for the flux snapshot tests and manual checks. The including file must
already have imported ClimaComms, ClimaCore as CC, Oceananigans as OC and
ClimaCoupler, and set `CMIPExt`.
=#

exchange_arch() = ClimaComms.device() isa ClimaComms.CUDADevice ? OC.GPU() : OC.CPU()

"""A 4 × 4-element cubed sphere with 4 GLL points: tiny, but Nq ≥ 3 like production."""
tiny_boundary_space(FT = Float64) = CC.CommonSpaces.CubedSphereSpace(
    FT;
    radius = FT(6.371e6),
    n_quad_points = 4,
    h_elem = 4,
    context = ClimaComms.context(),
)

"""A 36 × 18 tripolar ocean, dry for longitudes in [0°, 180°), so it has coastlines."""
function tiny_coastal_grid(arch)
    underlying = OC.TripolarGrid(
        arch;
        size = (36, 18, 1),
        southernmost_latitude = -80,
        north_poles_latitude = 55,
        first_pole_longitude = 70,
        fold_topology = OC.RightCenterFolded,
        z = (-100.0, 0.0),
        halo = (4, 4, 4),
    )
    bottom(x, y) = mod(x, 360) < 180 ? 100.0 : -200.0
    return OC.ImmersedBoundaryGrid(
        underlying,
        OC.GridFittedBottom(bottom);
        active_cells_map = false,
    )
end

# Stand-in for an Oceananigans model: capture reads only these fields.
struct SnapshotFakeOceanModel{G, T, V}
    grid::G
    tracers::T
    velocities::V
end

"""
    fake_ocean_setup(boundary_space, grid)

An exchange grid between `boundary_space` and `grid`, an ocean simulation stub
carrying it (with surface-flux BCs on T, S, u, v), and coupler fields on
`boundary_space`, all zeroed. Returns `(; sim, csf, eg_cpu)`.
"""
function fake_ocean_setup(boundary_space, grid)
    FT = CC.Spaces.undertype(boundary_space)
    arch = exchange_arch()
    eg_cpu = CMIPExt.build_exchange_grid(boundary_space, grid)
    eg = CMIPExt.on_device(arch, eg_cpu)

    function with_top_flux(FieldType, LX, LY)
        bcs = OC.FieldBoundaryConditions(
            grid,
            (LX(), LY(), OC.Center());
            top = OC.FluxBoundaryCondition(OC.Field{LX, LY, Nothing}(grid)),
        )
        return FieldType(grid; boundary_conditions = bcs)
    end
    model = SnapshotFakeOceanModel(
        grid,
        (;
            T = with_top_flux(OC.CenterField, OC.Center, OC.Center),
            S = with_top_flux(OC.CenterField, OC.Center, OC.Center),
        ),
        (;
            u = with_top_flux(OC.XFaceField, OC.Face, OC.Center),
            v = with_top_flux(OC.YFaceField, OC.Center, OC.Face),
        ),
    )

    boundary_zeros() = CC.Fields.zeros(boundary_space)
    flux_scratch = (;
        F_turb_ρτxz = boundary_zeros(),
        F_turb_ρτyz = boundary_zeros(),
        F_sh = boundary_zeros(),
        F_lh = boundary_zeros(),
        F_turb_moisture = boundary_zeros(),
    )
    remapping = (;
        use_exchange_grid = true,
        exchange_grid = eg,
        ocean_flux_state = CMIPExt.ExchangeFluxState{FT}(arch, eg_cpu.n_poly),
        flux_scratch,
        flux_dss_buffer = ClimaCoupler.Utilities.init_dss_buffer(flux_scratch.F_sh),
        weight_cov_scratch = boundary_zeros(),
        temp_uv_vec = CC.Fields.Field(CC.Geometry.UVVector{FT}, boundary_space),
        uv_basis = CMIPExt.uv_basis_coefficients(boundary_space),
    )
    ocean_properties = (;
        reference_density = FT(1020),
        heat_capacity = FT(3995),
        σ = FT(5.67e-8),
        C_to_K = FT(273.15),
    )
    ice_concentration = OC.Field{OC.Center, OC.Center, Nothing}(grid)
    sim = CMIPExt.OceananigansSimulation(
        (; model),
        nothing,
        ocean_properties,
        remapping,
        ice_concentration,
        nothing,
    )

    names = (
        :F_sh,
        :F_lh,
        :F_turb_moisture,
        :F_turb_ρτxz,
        :F_turb_ρτyz,
        :SW_d,
        :LW_d,
        :P_liq,
        :P_snow,
        :ocean_area_fraction,
        :ice_area_fraction,
        :land_area_fraction,
    )
    csf = NamedTuple{names}(ntuple(_ -> boundary_zeros(), length(names)))
    return (; sim, csf, eg_cpu)
end
