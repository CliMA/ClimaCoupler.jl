# list of exchanged fields

"""
    CplFieldInfo(name::Val{Symbol})

Set a CplFieldInfo
"""
function CplFieldInfo(name::Val{:T_S}) # ; space = false) # option for exchange grid

    writer_model = ThermalSlab()
    reader_model = DiffusiveColumn()

    regrid_map = :dummy_regrid

    # how to obtain a variable
    function get_reader_path(cs)
        f1 = cs.model_sims.slab.integrator.u.T
        f2 = cs.model_sims.slab.integrator.u.T
        combine_surfaces!(cs.combined_field, cs.surface_masks, (f1, f2))
        regrid(cs, getproperty(cs.regrid_maps, name))
    end

    # how to write a variable
    function get_writer_path(cs)
        cs.model_sims.slab.integrator.u.T
    end

    get_coupler_path = get_writer_path # where the coupler stores the variable
    _CplFieldInfo(name, get_coupler_path, writer_model, get_writer_path, reader_model, get_reader_path, regrid_map)
end

function CplFieldInfo(name::Val{:F_A}) # ; space = false) # option for exchange grid

    writer_model = DiffusiveColumn()
    reader_model = ThermalSlab()

    regrid_map = :dummy_regrid

    # how to obtain a variable
    function get_reader_path(cs)
        cs.model_sims.col.integrator.u.flux_accumulated / cs.dt
    end

    # how to write a variable
    function get_writer_path(cs)
        regrid(cs, cs.regrid_map)
        cs.model_sims.slab.integrator.u.T
    end

    get_coupler_path(cs) = cs.fields.:name

    _CplFieldInfo(name, get_coupler_path, writer_model, get_writer_path, reader_model, get_reader_path, regrid_map)
end





