## Plot VTK output on the 2D cubed-sphere # from ClimaCoreWorkshops
# building on the ClimaCore Workshop
using ClimaCoreVTK


# Plot the time series
times = 0:saveat:t_end

function R_ne_Nq_from_space(h_space)
    R = h_space.topology.mesh.domain.radius
    ne = h_space.topology.mesh.nev
    Nq = Spaces.Quadratures.polynomial_degree(h_space.quadrature_style) + 1
    return R, ne, Nq
end
function convert_to_type(TT, h_space, Sol_qt)
    _, ne, Nq = R_ne_Nq_from_space(h_space)
    horizontal_mesh64 = cubed_sphere_mesh(; radius = TT(Planet.planet_radius(params)), h_elem = ne)
    h_space64 = make_horizontal_space(horizontal_mesh64, quad, nothing)

    z_stretch = Meshes.GeneralizedExponentialStretching(TT(500), TT(5000))
    z_max = TT(30e3)
    z_elem = 10
    center_space64, _ = make_hybrid_spaces(h_space64, z_max, z_elem, z_stretch)

    Sol_qt64 = ones(center_space64)
    Sol_qt64_data = convert.(TT, parent(Sol_qt))
    parent(Sol_qt64) .= Sol_qt64_data
    return Sol_qt64
end

Sol_qt = map(u -> convert_to_type(Float64, h_space, u.c.ρq_tot), sol_atm.u)
Sol_qt_5 = map(x -> Fields.level(x, 1), Sol_qt)


function convert_to_type_T_sfc(TT, h_space, slab_u, slab_ocean_u)
    _, ne, Nq = R_ne_Nq_from_space(h_space)
    horizontal_mesh64 = cubed_sphere_mesh(; radius = TT(Planet.planet_radius(params)), h_elem = ne)
    h_space64 = make_horizontal_space(horizontal_mesh64, quad, nothing)

    Sol_qt64 = ones(h_space64)
    Sol_T_sfc_comb_data = combine_surface.(parent(mask), parent(slab_u.T_sfc), parent(slab_ocean_u.T_sfc))

    Sol_qt64_data = convert.(TT, parent(Sol_T_sfc_comb_data))
    parent(Sol_qt64) .= Sol_qt64_data
    return Sol_qt64
end

Sol_Tsfc = map(i -> convert_to_type_T_sfc(Float64, h_space, sol_slab.u[i], sol_slab_ocean.u[i]), 1:length(times))

function convert_to_type_mask(TT, h_space, mask)
    _, ne, Nq = R_ne_Nq_from_space(h_space)
    horizontal_mesh64 = cubed_sphere_mesh(; radius = TT(Planet.planet_radius(params)), h_elem = ne)
    h_space64 = make_horizontal_space(horizontal_mesh64, quad, nothing)

    Sol_qt64 = ones(h_space64)
    Sol_T_sfc_comb_data = parent(mask)

    Sol_qt64_data = convert.(TT, parent(Sol_T_sfc_comb_data))
    parent(Sol_qt64) .= Sol_qt64_data
    return Sol_qt64
end

con_nans(mask64) = (mask64 < 0.3 ? NaN : mask64)

mask64 = clean_mask.(Float64, convert_to_type_mask(Float64, h_space, mask))
# parent(mask64) .= con_nans.(parent(mask64)) # TODO: writepvd doesnt deal well with nans or Float32
mask_tsrs = map(i -> mask64, 1:length(times))

ClimaCoreVTK.writepvd(
    joinpath("data", "tt_5"),
    times,
    (qt = Sol_qt_5, T_sfc = Sol_Tsfc, mask = mask_tsrs),
    basis = :lagrange,
)




domain = Domains.SphereDomain()
mesh = Meshes.EquiangularCubedSphere(domain, ne)
grid_topology = Topologies.Topology2D(mesh)
quad = Spaces.Quadratures.GLL{Nq}()
space = Spaces.SpectralElementSpace2D(grid_topology, quad)






## Plot on a lat-long grid
ClimaCoreVTK.writepvd(joinpath(path, "humidity_lat_long"), times, (ω = Sol_vort,); latlong = true, basis = :point)
