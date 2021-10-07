# Land put and get calls
"""
    preLand(csolver::CplSolver)

Updates the Land model's prescribed boundary flux with `BoundaryEnergyFlux`.
"""
function preLand(csolver::CplSolver)
    mLand = csolver.component_list.domainLand.component_model
    mAtmos = csolver.component_list.domainAtmos.component_model

    mLand.odesolver.rhs!.state_auxiliary.F_ρθ_prescribed[mLand.boundary] .= 
        coupler_get(csolver.coupler, :BoundaryEnergyFlux, mLand.grid, DateTime(0), u"J")

end

"""
    postLand(csolver::CplSolver)

Updates coupler field `LandSurfaceTemerature`.
"""
function postLand(csolver::CplSolver)
    mLand = csolver.component_list.domainLand.component_model
    mAtmos = csolver.component_list.domainAtmos.component_model
    coupler_put!(csolver.coupler, :LandSurfaceTemerature, mLand.state.T_sfc[mLand.boundary], mLand.grid.numerical, DateTime(0), u"K")
end

# Atmos put and get calls
"""
    preAtmos(csolver::CplSolver)

Gets the `LandSurfaceTemerature` field for the Atmos model and resets the accumulated flux.
"""
function preAtmos(csolver::CplSolver)
    mLand = csolver.component_list.domainLand.component_model
    mAtmos = csolver.component_list.domainAtmos.component_model
    # Set boundary T_sfc used in atmos at the start of the coupling cycle.
    mAtmos.odesolver.rhs!.state_auxiliary.T_sfc[mAtmos.boundary] .= 
        coupler_get(csolver.coupler, :LandSurfaceTemerature, mAtmos.grid.numerical, DateTime(0), u"K")
    # Set atmos boundary flux accumulator to 0.
    mAtmos.state.F_ρθ_accum .= 0

    # get indicies of the required variables from the `mAtmos` MPIStateArray
    idx = varsindex(vars(mAtmos.state), :ρθ)[1]
    idx_rho = varsindex(vars(mAtmos.state), :ρ)[1]

    # Calculate, print and log domain-averaged energy
    FT = eltype(mAtmos.state)
    p = cpl_solver.component_list.domainAtmos.component_model.model.parameters # maybe the parameter list could be saved in the coupler? (this would requre all the compute kernels below to import the coupler)
    nel = mAtmos.grid.resolution.elements.vertical
    po = mAtmos.grid.resolution.polynomial_order.vertical

    horiz_sfc_area = p.xmax * p.ymax

    E_Land = weightedsum(mLand.state, 1) .* p.ρ_s .* p.h_s .* p.c_s ./ p.zmax ./ horiz_sfc_area # J / m^2
    E_Atmos = weightedsum(mAtmos.state, idx) .* p.cp_d ./ horiz_sfc_area ./  nel ./ (po+1) ./ (po+2)  # J / m^2 

    @info(
        "preatmos",
        time = string(csolver.t) * "/" * string(mAtmos.time.finish),
        total_energyA = E_Land,
        total_energyB = E_Atmos,
        total_energy = E_Atmos + E_Land,
        land_T_sfc = maximum(mLand.state.T_sfc[mLand.boundary]),
        atmos_ρθ_sfc = maximum(mAtmos.state.ρθ[mAtmos.boundary]),
    )

    isnothing(csolver.fluxlog) ? nothing : csolver.fluxlog.A[csolver.steps] = E_Land
    isnothing(csolver.fluxlog) ? nothing : csolver.fluxlog.B[csolver.steps] = E_Atmos

end

"""
    postAtmos(csolver::CplSolver)

Updates coupler field `BoundaryEnergyFlux` with the average flux at the atmosphere's boundary.
"""
function postAtmos(csolver::CplSolver)
    mLand = csolver.component_list.domainLand.component_model
    mAtmos = csolver.component_list.domainAtmos.component_model
    # Pass atmos exports to "coupler" namespace
    # 1. Save mean θ flux at the Atmos boundary during the coupling period
    coupler_put!(csolver.coupler, :BoundaryEnergyFlux, mAtmos.state.F_ρθ_accum[mAtmos.boundary] ./ csolver.dt,
        mAtmos.grid.numerical, DateTime(0), u"J")
end
