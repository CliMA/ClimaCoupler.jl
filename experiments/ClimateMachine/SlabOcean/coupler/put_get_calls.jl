# Ocean put and get calls

"""
function preOcean(csolver)
    - saves and regrids the `BoundaryEnergyFlux` couplerfield (i.e. regridded `mAtmos.state.F_ρθ_accum[mAtmos.boundary]`) to `F_ρθ_prescribed[mOcean.boundary]` on the `mOcean` grid
    csolver::CplSolver
"""
function preOcean(csolver)
    mOcean = csolver.component_list.domainOcean.component_model
    mAtmos = csolver.component_list.domainAtmos.component_model

    mOcean.odesolver.rhs!.state_auxiliary.F_ρe_prescribed[mOcean.boundary] .=
        coupler_get(csolver.coupler, :BoundaryEnergyFlux, mOcean.grid, DateTime(0), u"J")
end

"""
function postOcean(csolver)
    - updates couplerfield `OceanSurfaceTemerature` with `mOcean.state.T_sfc[mOcean.boundary]` regridded to the coupler grid, and updates the coupler time
    csolver::CplSolver
"""
function postOcean(csolver)
    mOcean = csolver.component_list.domainOcean.component_model
    mAtmos = csolver.component_list.domainAtmos.component_model
    coupler_put!(
        csolver.coupler,
        :SeaSurfaceTemerature,
        mOcean.state.T_sfc[mOcean.boundary],
        mOcean.grid.numerical,
        DateTime(0),
        u"K",
    )
end

# Atmos put and get calls
"""
function preAtmos(csolver)
    - saves and regrids the `SeaSurfaceTemerature` couplerfield (i.e. regridded `mOcean.state.T_sfc[mOcean.boundary]``) on coupler grid to `mAtmos.aux.T_sfc[mAtmos.boundary]` on the `mAtmos` grid
    - resets the accumulated flux `F_ρe_accum` in `mAtmos` to 0
    csolver::CplSolver
"""
function preAtmos(csolver)
    mOcean = csolver.component_list.domainOcean.component_model
    mAtmos = csolver.component_list.domainAtmos.component_model
    # Set boundary T_sfc used in atmos at the start of the coupling cycle.
    mAtmos.odesolver.rhs!.state_auxiliary.T_sfc[mAtmos.boundary] .=
        coupler_get(csolver.coupler, :SeaSurfaceTemerature, mAtmos.grid.numerical, DateTime(0), u"K")
    # Set atmos boundary flux accumulator to 0.
    mAtmos.state.F_ρe_accum .= 0

    # get indicies of the required variables from the `mAtmos` MPIStateArray
    idx = varsindex(vars(mAtmos.state), :ρe)[1]
    idx_rho = varsindex(vars(mAtmos.state), :ρ)[1]

    # Calculate, print and log domain-averaged energy
    FT = eltype(mAtmos.state)

    model_type = cpl_solver.component_list.domainAtmos.component_model.model[1] # TODO: bl should be more accessible
    p = model_type.model.physics.parameters # maybe the parameter list could be saved in the coupler? (this would requre all the compute kernels below to import the coupler)
    nel = mAtmos.grid.resolution.elements.vertical
    po = mAtmos.grid.resolution.polynomial_order.vertical

    Earth_sfc_area = FT(4π * 6371000.0^2)
    ocean_layer_volume = sum((mOcean.state.weights[:, :, mOcean.state.realelems])[mOcean.boundary])
    E_Ocean = weightedsum(mOcean.state, 1) .* p.ρ_o .* p.c_o .* p.h_o ./ ocean_layer_volume # J / m^2
    E_Atmos = weightedsum(mAtmos.state, idx) ./ Earth_sfc_area  # J / m^2 

    @info(
        "preatmos",
        time = string(csolver.t) * "/" * string(mAtmos.time.finish),
        total_energyA = E_Ocean,
        total_energyB = E_Atmos,
        total_energy = E_Atmos + E_Ocean,
        Ocean_T_sfc = maximum(mOcean.state.T_sfc[mOcean.boundary]),
        atmos_ρe_sfc = maximum(mAtmos.state.ρe[mAtmos.boundary]),
    )

    isnothing(csolver.fluxlog) ? nothing : csolver.fluxlog.A[csolver.steps] = E_Ocean
    isnothing(csolver.fluxlog) ? nothing : csolver.fluxlog.B[csolver.steps] = E_Atmos
end

"""
function postAtmos(csolver)
    - updates couplerfield `BoundaryEnergyFlux` with mAtmos.state.F_ρe_accum[mAtmos.boundary] regridded to the coupler grid, and updates the coupler time
    csolver::CplSolver
"""
function postAtmos(csolver)
    mOcean = csolver.component_list.domainOcean.component_model
    mAtmos = csolver.component_list.domainAtmos.component_model

    # Pass atmos exports to "coupler" namespace
    # 1. Save mean θ flux at the Atmos boundary during the coupling period

    @info(
        "postatmos",
        time = string(csolver.t) * "/" * string(mAtmos.time.finish),
        Ocean_T_sfc = maximum(mOcean.state.T_sfc[mOcean.boundary]),
        atmos_ρe_sfc = maximum(mAtmos.state.ρe[mAtmos.boundary]),
    )

    coupler_put!(
        csolver.coupler,
        :BoundaryEnergyFlux,
        mAtmos.state.F_ρe_accum[mAtmos.boundary] ./ csolver.dt,
        mAtmos.grid.numerical,
        DateTime(0),
        u"J",
    )
end
