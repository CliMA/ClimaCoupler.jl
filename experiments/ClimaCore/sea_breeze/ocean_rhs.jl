# slab ocean ODE
function ocn_rhs!(du, u, (parameters, F_accumulated), t)
    """
    Slab layer equation
        d(T_sfc)/dt = - (F_accumulated) / (h_ocn * ρ_ocn * c_ocn)
        where 
            F_accumulated = F_integrated / Δt_coupler
    """
    @unpack ocn_h, ocn_ρ, ocn_c = parameters
    @unpack T_sfc = du

    @. T_sfc = (-F_accumulated) / (ocn_h * ocn_ρ * ocn_c)
end

# set up domain
function hspace_1D(xlim = (-π, π), npoly = 0, helem = 10)
    FT = Float64

    domain = Domains.IntervalDomain(Geometry.XPoint{FT}(xlim[1]) .. Geometry.XPoint{FT}(xlim[2]), periodic = true)
    mesh = Meshes.IntervalMesh(domain; nelems = helem)
    topology = Topologies.IntervalTopology(mesh)

    # Finite Volume Approximation: Gauss-Lobatto with 1pt per element
    quad = Spaces.Quadratures.GL{npoly + 1}()
    space = Spaces.SpectralElementSpace1D(topology, quad)

    return space
end

# init simulation
function ocn_init(; xmin = -1000, xmax = 1000, helem = 20, npoly = 0)

    # construct domain spaces - get only surface layer (NB: z should be zero, not z = first central height)
    space = hspace_1D((xmin, xmax), npoly, helem)
    coords = Fields.coordinate_field(space)
    domain = space

    # initial condition
    T_sfc = map(coords) do coord
        T_sfc = 267.0
    end

    # prognostic variable
    Y = Fields.FieldVector(T_sfc = T_sfc)

    return Y, domain
end
