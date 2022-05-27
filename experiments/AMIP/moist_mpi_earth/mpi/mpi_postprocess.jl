
# function mpi_collapse_solu()# (comms_ctx, horizontal_mesh, quad, z_max, z_elem, z_stretch, Y) # TODO: it's using the global scope here
if is_distributed # replace sol.u on the root processor with the global sol.u

    if ClimaComms.iamroot(comms_ctx)
        global_h_space = make_horizontal_space(horizontal_mesh, quad, nothing)
        global_center_space, global_face_space = make_hybrid_spaces(global_h_space, z_max, z_elem, z_stretch)
        global_Y_c_type = Fields.Field{typeof(Fields.field_values(Y.c)), typeof(global_center_space)}
        global_Y_f_type = Fields.Field{typeof(Fields.field_values(Y.f)), typeof(global_face_space)}
        global_Y_type = Fields.FieldVector{FT, NamedTuple{(:c, :f), Tuple{global_Y_c_type, global_Y_f_type}}}
        global_Y_slab_type =
            Fields.Field{typeof(Fields.field_values(bucket_sim.integrator.u.bucket.T_sfc)), typeof(global_h_space)}
        global_sol_u_atmos = similar(sol_atm.u, global_Y_type)
        global_sol_u_slab = similar(sol_slab.u, global_Y_slab_type)
    end
    for i in 1:length(sol_atm.u)
        global_Y_c = DataLayouts.gather(comms_ctx, Fields.field_values(sol_atm.u[i].c))
        global_Y_f = DataLayouts.gather(comms_ctx, Fields.field_values(sol_atm.u[i].f))
        global_Y_slab = DataLayouts.gather(comms_ctx, Fields.field_values(sol_slab.u[i].T_sfc))
        if ClimaComms.iamroot(comms_ctx)
            global_sol_u_atmos[i] = Fields.FieldVector(
                c = Fields.Field(global_Y_c, global_center_space),
                f = Fields.Field(global_Y_f, global_face_space),
            )
            global_sol_u_slab[i] = Fields.Field(global_Y_slab, global_h_space)
        end
    end
    if ClimaComms.iamroot(comms_ctx)
        solu_atm = global_sol_u_atmos
        solu_slab = global_sol_u_slab
        return solu_atm, solu_slab
    end
else
    solu_atm = sol_atm.u
    h_space = make_horizontal_space(horizontal_mesh, quad, nothing) #TODO move this to the beginning (once same the instance error sorted)
    solu_slab = Fields.FieldVector(T_sfc = [Fields.Field(Fields.field_values(u.bucket.T_sfc), h_space) for u in sol_slab.u])
end
