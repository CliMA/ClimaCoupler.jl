#=
# Flux snapshot tests

Capture, file writing, NaN dump, config and plotting for flux snapshots
(`SimOutput.write_flux_snapshot`, `Plotting.plot_flux_snapshot`).

    julia --project=experiments/CMIP experiments/test/flux_snapshots_test.jl
=#
import Test: @test, @testset
import Dates
import JLD2
import ClimaComms
ClimaComms.@import_required_backends
import ClimaCore as CC
import Oceananigans as OC
import ClimaSeaIce
import ClimaAtmos
import KernelAbstractions
import ConservativeRegridding
import Adapt
import CairoMakie, GeoMakie, Makie, Poppler_jll, Printf
import ClimaCoupler
import ClimaCoupler: Input, Plotting, SimOutput, TimeManager

CMIPExt = Base.get_extension(ClimaCoupler, :ClimaCouplerCMIPExt)
@assert !isnothing(CMIPExt)
@assert !isnothing(Base.get_extension(ClimaCoupler, :ClimaCouplerMakieExt))
ClimaComms.init(ClimaComms.context())
FT = Float64 # not `const`, so the file can be re-included in a REPL
include("flux_snapshot_fixtures.jl")
boundary_space = tiny_boundary_space(FT)

@testset "Flux snapshot schedule" begin
    start_date = Dates.DateTime(2010, 1, 1)
    @test isnothing(
        TimeManager.flux_snapshot_schedule("never", 0.0, 86400.0, start_date, 0.0),
    )

    # Hourly within [2 h, 4 h]. The calendar schedule keeps ticking outside the
    # window, so it is still on the hour inside it.
    schedule =
        TimeManager.flux_snapshot_schedule("1hours", 7200.0, 14400.0, start_date, 0.0)
    fired = [t for t in 360.0:360.0:21600.0 if schedule((; t, step = round(Int, t / 360)))]
    @test fired == [7200.0, 10800.0, 14400.0]
end

@testset "Exchange-grid geometry" begin
    grid = tiny_coastal_grid(OC.CPU())
    eg = CMIPExt.build_exchange_grid(boundary_space, grid)
    geo = CMIPExt.exchange_grid_geometry(boundary_space, grid, eg)
    (; poly, nodes, cells) = geo

    # One real polygon ring per exchange-grid polygon.
    @test length(poly.vert_ptr) == eg.n_poly + 1
    @test all(>=(3), diff(poly.vert_ptr))
    @test length(poly.vert_lon) == length(poly.vert_lat) == poly.vert_ptr[end] - 1
    @test all(lat -> -90 <= lat <= 90, poly.vert_lat)
    @test poly.area == eg.area

    # Same numbering as the exchange grid, checked exactly: the spherical area of
    # each vertex ring (fan of Van Oosterom–Strackee triangles) matches the
    # exchange grid's area for that polygon. A ring from a different polygon,
    # even a neighbor, would not.
    dot3(a, b) = a[1] * b[1] + a[2] * b[2] + a[3] * b[3]
    cross3(a, b) =
        (a[2] * b[3] - a[3] * b[2], a[3] * b[1] - a[1] * b[3], a[1] * b[2] - a[2] * b[1])
    unit_vector(lon, lat) = (cosd(lat) * cosd(lon), cosd(lat) * sind(lon), sind(lat))
    function ring_area(k)
        r = poly.vert_ptr[k]:(poly.vert_ptr[k + 1] - 1)
        p = [unit_vector(poly.vert_lon[v], poly.vert_lat[v]) for v in r]
        excess = sum(2:(length(p) - 1)) do v
            a, b, c = p[1], p[v], p[v + 1]
            2 * atan(abs(dot3(a, cross3(b, c))), 1 + dot3(a, b) + dot3(b, c) + dot3(c, a))
        end
        return excess * 6.371e6^2
    end
    atol = 1e-8 * maximum(eg.area) # degenerate slivers carry zero area in `eg`
    @test all(k -> isapprox(ring_area(k), eg.area[k]; rtol = 1e-8, atol), 1:(eg.n_poly))

    # Nodes: flat SE order, weights that integrate the space.
    @test length(nodes.lon) == eg.n_nodes
    @test sum(nodes.Jw) ≈ sum(CC.Fields.ones(boundary_space)) rtol = 1e-12
    @test nodes.poly_ptr == eg.snode_ptr

    # Cells: one column per FV cell, wet area matching the polygons.
    @test size(cells.corner_lon) == (4, 36 * 18)
    @test cells.wet_area == eg.oc_wet_area
    @test sum(cells.wet_area) ≈ sum(eg.area)
end

@testset "Ocean flux snapshot capture" begin
    (; sim, csf, eg_cpu) =
        fake_ocean_setup(boundary_space, tiny_coastal_grid(exchange_arch()))
    remapping = sim.remapping
    fs = remapping.ocean_flux_state
    uniform = (;
        F_sh = 10.0,
        F_lh = 20.0,
        F_moisture = 1e-5,
        F_τu = 0.1,
        F_τv = -0.2,
        SW_d = 300.0,
        LW_d = 350.0,
        P_liq = -1e-5,
        P_snow = -2e-6,
        T_sfc = 290.0,
        sic = 0.25,
    )
    for (name, value) in pairs(uniform)
        fill!(getproperty(fs, name), value)
    end
    csf.F_sh .= 7
    csf.ocean_area_fraction .= 0.5
    uv = CC.Fields.Field(CC.Geometry.UVVector{FT}, boundary_space)
    parent(uv.components.data.:1) .= 0.3
    parent(uv.components.data.:2) .= -0.4
    CMIPExt.cartesian_to_contravariant!(
        csf.F_turb_ρτxz,
        csf.F_turb_ρτyz,
        uv,
        remapping.uv_basis,
    )
    fill!(parent(CMIPExt.surface_flux(sim.ocean.model.tracers.T)), 1e-6)

    snap = SimOutput.capture_flux_snapshot(sim, csf)

    # Exchange grid: this step's per-polygon values.
    @test all(==(10.0), snap.ocean.poly.F_sh)
    @test all(==(300.0), snap.ocean.poly.SW_d)
    @test all(==(0.75), snap.ocean.poly.weight)

    # SE nodes, ocean only: a weighted average of a constant is that constant
    # wherever the node received flux.
    covered = snap.ocean.nodes.cov .> 1e-3
    @test count(covered) > 0
    @test all(isapprox.(snap.ocean.nodes.F_sh[covered], 10.0; rtol = 1e-10))
    @test all(isapprox.(snap.ocean.nodes.F_τu[covered], 0.1; rtol = 1e-10))

    # SE nodes, combined: as the atmosphere receives them, momentum east/north.
    @test all(==(7.0), snap.atmos.F_sh)
    @test all(==(0.5), snap.atmos.ocean_fraction)
    @test all(isapprox.(snap.atmos.F_τu, 0.3; rtol = 1e-10))
    @test all(isapprox.(snap.atmos.F_τv, -0.4; rtol = 1e-10))

    # Ocean cells: open-water weighted, and the scatter conserves the integral.
    wet = eg_cpu.oc_wet_area .> 0
    @test all(isapprox.(snap.ocean.cells.F_sh[wet], 0.75 * 10.0; rtol = 1e-10))
    @test sum(eg_cpu.area .* snap.ocean.poly.weight .* snap.ocean.poly.F_sh) ≈
          sum(eg_cpu.oc_wet_area .* snap.ocean.cells.F_sh)

    # What the ocean model receives, surface block only.
    @test all(==(1e-6), snap.ocean.totals.T)
    @test length(snap.ocean.context.SST) == 36 * 18

    # Components without an exchange grid capture nothing.
    @test isnothing(SimOutput.capture_flux_snapshot(nothing, csf))
end

@testset "Flux snapshot files and NaN dump" begin
    (; sim, csf, eg_cpu) =
        fake_ocean_setup(boundary_space, tiny_coastal_grid(exchange_arch()))
    fill!(sim.remapping.ocean_flux_state.F_sh, 10.0)
    model_sims = (; ocean_sim = sim)
    dir = mktempdir()
    date = Dates.DateTime(2010, 1, 1, 6)

    @test SimOutput.flux_snapshot_dir("out/clima_coupler") ==
          joinpath("out/clima_coupler", "fluxes")

    path = SimOutput.write_flux_snapshot(model_sims, csf, dir, date)
    @test basename(path) == "fluxes_2010-01-01T060000.jld2"
    JLD2.jldopen(path, "r") do file
        @test file["date"] == "2010-01-01T06:00:00"
        @test all(==(10.0), file["ocean/poly/F_sh"])
        @test length(file["atmos/F_sh"]) == eg_cpu.n_nodes
        @test !haskey(file, "ice")
    end
    geometry_path = joinpath(dir, "exchange_grid_geometry.jld2")
    JLD2.jldopen(geometry_path, "r") do file
        @test length(file["poly/vert_ptr"]) == eg_cpu.n_poly + 1
    end
    geometry_mtime = mtime(geometry_path)
    SimOutput.write_flux_snapshot(model_sims, csf, dir, date + Dates.Hour(1))
    @test mtime(geometry_path) == geometry_mtime # written once

    # Nothing to capture: no file.
    @test isnothing(SimOutput.write_flux_snapshot((;), csf, dir, date))

    # NaN dump: one file on the first NaN, then nothing.
    nan_dump = SimOutput.NaNFluxSnapshot()
    @test isnothing(
        SimOutput.maybe_write_nan_snapshot!(nan_dump, model_sims, csf, dir, date),
    )
    fill!(parent(csf.F_lh), NaN)
    @test SimOutput.turbulent_fluxes_have_nan(csf)
    nan_path = SimOutput.maybe_write_nan_snapshot!(nan_dump, model_sims, csf, dir, date)
    @test endswith(nan_path, "_nan.jld2") && isfile(nan_path)
    later = date + Dates.Hour(2)
    @test isnothing(
        SimOutput.maybe_write_nan_snapshot!(nan_dump, model_sims, csf, dir, later),
    )
    @test count(endswith("_nan.jld2"), readdir(dir)) == 1
end

@testset "Flux snapshot config" begin
    config_file = joinpath(
        pkgdir(ClimaCoupler),
        "config",
        "ci_configs",
        "cmip_oceananigans_climaseaice.yml",
    )
    config_dict = Input.get_coupler_config_dict(config_file)
    args = Input.get_coupler_args(config_dict)
    @test args.flux_snapshot_interval == "never"
    @test args.flux_snapshot_start == float(args.t_start)
    @test args.flux_snapshot_end == float(args.t_end)
    @test args.flux_snapshot_on_nan

    config_dict["flux_snapshot_interval"] = "6hours"
    config_dict["flux_snapshot_start"] = "2days"
    config_dict["flux_snapshot_end"] = "3days"
    config_dict["flux_snapshot_on_nan"] = false
    args = Input.get_coupler_args(config_dict)
    @test args.flux_snapshot_interval == "6hours"
    @test args.flux_snapshot_start == 2 * 86400.0
    @test args.flux_snapshot_end == 3 * 86400.0
    @test !args.flux_snapshot_on_nan
end

@testset "Flux snapshot plots" begin
    (; sim, csf) = fake_ocean_setup(boundary_space, tiny_coastal_grid(exchange_arch()))
    fs = sim.remapping.ocean_flux_state
    fill!(fs.F_sh, 10.0)
    fill!(fs.F_τu, 0.1)
    dir = mktempdir()
    path = SimOutput.write_flux_snapshot(
        (; ocean_sim = sim),
        csf,
        dir,
        Dates.DateTime(2010, 1, 1, 6),
    )
    pngs = Plotting.plot_flux_snapshot(
        path;
        fluxes = ("F_sh", "F_τu"),
        regions = (nothing, (0.0, 90.0, 3000.0)),
    )
    @test length(pngs) == 4
    @test all(isfile, pngs)
    @test any(endswith("F_sh_global.png"), pngs)
end
