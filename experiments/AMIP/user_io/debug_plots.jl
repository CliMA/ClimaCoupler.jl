using Plots
using ClimaCorePlots
using Printf
using ClimaCoupler.Interfacer: ComponentModelSimulation, SurfaceModelSimulation
using ClimaCoupler: TestHelper
using ClimaCore

# plotting functions for the coupled simulation
"""
    debug(cs::CoupledSimulation, dir = "debug", cs_fields_ref = nothing)

Plot the fields of a coupled simulation and save plots to a directory.
"""
function debug(cs::CoupledSimulation, dir = "debug", cs_fields_ref = nothing)
    mkpath(dir)
    # @info "plotting debug in " * dir
    for sim in cs.model_sims
        debug(sim, dir)
    end
    debug(cs.fields, dir, cs_fields_ref)
end

"""
    debug(cs_fields::NamedTuple, dir, cs_fields_ref = nothing)

Plot useful coupler fields (in `field_names`) and save plots to a directory.

If `cs_fields_ref` is provided (e.g., using a copy of cs.fields from the initialization),
plot the anomalies of the fields with respect to `cs_fields_ref`.
"""
function debug(cs_fields::NamedTuple, dir, cs_fields_ref = nothing)
    field_names = (
        :surface_albedo,
        :F_radiative,
        :F_turb_energy,
        :F_turb_moisture,
        :P_liq,
        :T_S,
        :ρ_sfc,
        :q_sfc,
        :beta,
        :z0b_S,
        :z0m_S,
    )
    all_plots = []
    cpu_comms_ctx = ClimaComms.SingletonCommsContext(ClimaComms.CPUSingleThreaded())

    for field_name in field_names
        field = getproperty(cs_fields, field_name)

        # Copy field onto cpu space
        space = axes(field)
        FT = ClimaCore.Spaces.undertype(space)
        R = get_R(space.grid)
        ne = get_ne(space.grid)
        polynomial_degree = ClimaCore.Quadratures.polynomial_degree(ClimaCore.Spaces.quadrature_style(space.grid))
        nz = ClimaCore.Spaces.nlevels(space)
        height = get_height(space.grid)

        cpu_space = TestHelper.create_space(
            FT,
            comms_ctx = cpu_comms_ctx,
            R = R,
            ne = ne,
            polynomial_degree = polynomial_degree,
            nz = nz,
            height = height,
        )
        cpu_field = ClimaCore.Fields.ones(cpu_space)

        parent(cpu_field) .= Array(parent(field))

        push!(all_plots, Plots.plot(cpu_field, title = string(field_name) * print_extrema(field)))
        if (field_name == :T_S) && (@isdefined debug_csf0)
            push!(
                all_plots,
                Plots.plot(
                    cpu_field .- debug_csf0.T_S,
                    title = string(field_name) * "_anom" * print_extrema(field),
                    color = :viridis,
                ),
            )
        end
    end
    Plots.plot(all_plots..., size = (1500, 800))
    Plots.png(dir * "debug_coupler")

    # plot anomalies if a reference cs.fields, `cs_fields_ref`, are provided
    if !isnothing(cs_fields_ref)
        all_plots = []
        for field_name in field_names
            field = getproperty(cs_fields, field_name)
            # Copy field onto cpu space
            parent(cpu_field) .= parent(field)

            push!(
                all_plots,
                Plots.plot(
                    cpu_field .- getproperty(cs_fields_ref, field_name),
                    title = string(field_name) * print_extrema(field),
                ),
            )
        end
        Plots.plot(all_plots..., size = (1500, 800))
        Plots.png(dir * "debug_coupler_amomalies")
    end
end

function get_ne(grid::ClimaCore.Grids.SpectralElementGrid2D)
    return grid.topology.mesh.ne
end
function get_ne(grid::ClimaCore.Grids.LevelGrid)
    return get_ne(grid.full_grid.horizontal_grid)
end
function get_ne(grid::ClimaCore.Grids.ExtrudedFiniteDifferenceGrid)
    return get_ne(grid.horizontal_grid)
end

function get_height(grid::ClimaCore.Grids.ExtrudedFiniteDifferenceGrid)
    return grid.vertical_grid.topology.mesh.domain.coord_max.z
end
function get_height(grid)
    return nothing # 2d case
end

function get_R(grid)#::ClimaCore.Grids.SpectralElementGrid2D)
    return ClimaCore.Grids.global_geometry(grid).radius
end
"""
    debug(sim::ComponentModelSimulation, dir)

Plot the fields of a component model simulation and save plots to a directory.
"""
function debug(sim::ComponentModelSimulation, dir)

    @show name(sim)
    field_names = plot_field_names(sim)
    cpu_comms_ctx = ClimaComms.SingletonCommsContext(ClimaComms.CPUSingleThreaded())

    all_plots = []
    for field_name in field_names
        field = get_field(sim, Val(field_name))
        # Copy field onto cpu space
        space = axes(field)
        FT = ClimaCore.Spaces.undertype(space)
        R = get_R(space.grid)
        ne = get_ne(space.grid)
        polynomial_degree = ClimaCore.Quadratures.polynomial_degree(ClimaCore.Spaces.quadrature_style(space.grid))
        nz = ClimaCore.Spaces.nlevels(space)
        height = get_height(space.grid)
        cpu_space = TestHelper.create_space(
            FT,
            comms_ctx = cpu_comms_ctx,
            R = R,
            ne = ne,
            polynomial_degree = polynomial_degree,
            nz = nz,
            height = height,
        )
        cpu_field = ClimaCore.Fields.ones(cpu_space)
        parent(cpu_field) .= Array(parent(field))

        push!(all_plots, Plots.plot(cpu_field, title = string(field_name) * print_extrema(field)))
    end
    fig = Plots.plot(all_plots..., size = (1500, 800))
    Plots.png(dir * "debug_$(name(sim))")

end

"""
    print_extrema(field::ClimaCore.Fields.Field)

Return the minimum and maximum values of a field as a string.
"""
function print_extrema(field::ClimaCore.Fields.Field)
    ext_vals = extrema(field)
    min = @sprintf("%.2E", ext_vals[1])
    max = @sprintf("%.2E", ext_vals[2])
    return " [$min, $max]"
end

# below are additional fields specific to this experiment (ourside of the required coupler fields) that we are interested in plotting for debugging purposes

# additional ClimaAtmos model debug fields
function get_field(sim::ClimaAtmosSimulation, ::Val{:w})
    w_c = ones(ClimaCore.Spaces.horizontal_space(sim.domain.face_space))
    parent(w_c) .= parent(ClimaCore.Fields.level(ClimaCore.Geometry.WVector.(sim.integrator.u.f.u₃), 5 .+ half))
    return w_c
end
get_field(sim::ClimaAtmosSimulation, ::Val{:ρq_tot}) = sim.integrator.u.c.ρq_tot
get_field(sim::ClimaAtmosSimulation, ::Val{:ρe_tot}) = sim.integrator.u.c.ρe_tot

# additional BucketSimulation debug fields
get_field(sim::BucketSimulation, ::Val{:σS}) = sim.integrator.u.bucket.σS
get_field(sim::BucketSimulation, ::Val{:Ws}) = sim.integrator.u.bucket.Ws
get_field(sim::BucketSimulation, ::Val{:W}) = sim.integrator.u.bucket.W

# currently selected plot fields
plot_field_names(sim::SurfaceModelSimulation) = (:area_fraction, :surface_temperature, :surface_humidity)
plot_field_names(sim::BucketSimulation) =
    (:area_fraction, :surface_temperature, :surface_humidity, :air_density, :σS, :Ws, :W)
plot_field_names(sim::ClimaAtmosSimulation) = (:w, :ρq_tot, :ρe_tot, :liquid_precipitation, :snow_precipitation)
