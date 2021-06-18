# Land put and get calls
"""
function preA(csolver)
    - saves and regrids the `EnergyFluxAtmos` couplerfield (i.e. regridded `mAtmos.state.F_ρθ_accum[mAtmos.boundary]`) to `F_ρθ_prescribed[mLand.boundary]` on the `mLand` grid
    csolver::CplSolver
"""
function preLand(csolver)
    mLand = csolver.component_list.domainLand.component_model
    mAtmos = csolver.component_list.domainAtmos.component_model

    mLand.odesolver.rhs!.state_auxiliary.F_ρθ_prescribed[mLand.boundary] .= 
        coupler_get(csolver.coupler, :EnergyFluxAtmos, mLand.grid, DateTime(0), u"J")


    # For conservation checks (TODO: find a better way to do 2D domain integrals):
    # temp =  coupler_get(csolver.coupler, :EnergyFluxAtmos, mLand.grid, DateTime(0), u"J")    
    # abs_temp = abs.(temp)
    # if sum(abs_temp) == 0
    #     sign = temp .* 0.0 .+ 1.0
    # else
    #     sign = temp ./ abs_temp
    # end

    # mLand.odesolver.rhs!.state_auxiliary.F_ρθ_prescribed[:,:,:] .= sign[1,1,1] .* maximum(abs_temp)
end

"""
function postA(csolver)
    - updates couplerfield `LandSurfaceTemerature` with `mLand.state.T_sfc[mLand.boundary]` regridded to the coupler grid, and updates the coupler time
    csolver::CplSolver
"""
function postLand(csolver)
    mLand = csolver.component_list.domainLand.component_model
    mAtmos = csolver.component_list.domainAtmos.component_model
    coupler_put!(csolver.coupler, :LandSurfaceTemerature, mLand.state.T_sfc[mLand.boundary], mLand.grid.numerical, DateTime(0), u"K")
end

# Atmos put and get calls
"""
function preAtmos(csolver)
    - saves and regrids the `LandSurfaceTemerature` couplerfield (i.e. regridded `mLand.state.T_sfc[mLand.boundary]``) on coupler grid to `mAtmos.aux.T_sfc[mAtmos.boundary]` on the `mAtmos` grid
    - resets the accumulated flux `F_ρθ_accum` in `mAtmos` to 0
    csolver::CplSolver
"""
function preAtmos(csolver)
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

    E_Land = weightedsum(mLand.state, 1) .* p.ρ_s .* p.h_s .* p.c_s ./ p.zmax ./ horiz_sfc_area ./ csolver.dt # W / m^2
    E_Atmos = weightedsum(mAtmos.state, idx) .* p.cp_d ./ horiz_sfc_area ./ csolver.dt  ./  nel ./ (po+1) ./ (po+2)  # W / m^2 
    #E_Land = weightedsum(mLand.state, 1) .* p.ρ_s .* p.h_s .* p.c_s  # J / m^2
    #E_Atmos = weightedsum(mAtmos.state, idx) .* p.cp_d .* p.zmax ./  nel ./ (po+1) ./ (po+2) # J / m^2 

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
function postAtmos(csolver)
    - updates couplerfield `EnergyFluxAtmos` with mAtmos.state.F_ρθ_accum[mAtmos.boundary] regridded to the coupler grid, and updates the coupler time
    csolver::CplSolver
"""
function postAtmos(csolver)
    mLand = csolver.component_list.domainLand.component_model
    mAtmos = csolver.component_list.domainAtmos.component_model
    # Pass atmos exports to "coupler" namespace
    # 1. Save mean θ flux at the Atmos boundary during the coupling period
    coupler_put!(csolver.coupler, :EnergyFluxAtmos, mAtmos.state.F_ρθ_accum[mAtmos.boundary] ./ csolver.dt,
        mAtmos.grid.numerical, DateTime(0), u"J")
end