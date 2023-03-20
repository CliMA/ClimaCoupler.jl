# list of exchanged fields

"""
    CoupledFieldInfo(name::Val{Symbol})

Set a CoupledFieldInfo
"""
# function CoupledFieldInfo(name::Val{:T_S}) # ; space = false) # option for exchange grid

#     target_model = ThermalSlab()
#     source_model = (DiffusiveColumn(), DiffusiveColumn())

#     # set regridding type
#     regrid_map = :dummy_regrid

#     # how to obtain a variable
#     get_source_path(cs) =  # how coupler views this field
#         (cs.model_sims.slab.integrator.u.T, cs.model_sims.slab.integrator.u.T)

#     # how to write a variable
#     get_target_path(cs) = cs.model_sims.slab.integrator.u.T

#     # set masks
#     get_source_mask(cs) = (0.5, 0.5)
#     get_target_mask(cs) = nothing

#     get_coupler_path = get_target_path # where the coupler stores the variable

#     _CoupledFieldInfo(name, get_coupler_path, target_model, get_target_path, source_model, get_source_path, regrid_map, mask)
# end

# function CoupledFieldInfo(name::Val{:F_A}) # ; space = false) # option for exchange grid

#     target_model = DiffusiveColumn()
#     source_model = ThermalSlab()

#     regrid_map = :dummy_regrid

#     # how to obtain a variable
#     function get_source_path(cs)
#         cs.model_sims.col.integrator.u.flux_accumulated / cs.dt
#     end

#     # how to write a variable
#     function get_target_path(cs)
#         regrid(cs, cs.regrid_map)
#         cs.model_sims.slab.integrator.u.T
#     end

#     get_coupler_path(cs) = cs.fields.:name

#     _CoupledFieldInfo(name, get_coupler_path, target_model, get_target_path, source_model, get_source_path, regrid_map)
# end





