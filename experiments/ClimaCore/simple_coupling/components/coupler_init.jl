#### move to general:
# abstract type AbstractComponentModel end # tag upon which to dispatch

# abstract type AbstractSimulation end  # Simulation struct containing all sim info

# struct ComponentSimulation{T, Y, I, A} <: AbstractSimulation
#     tag::T
#     Y_init::Y
#     integrator::I
#     area_fraction::A
# end

# function ComponentSimulation(tag; Y_init = nothing, integrator = nothing, area_fraction = nothing)
#     ComponentSimulation(tag, Y_init, integrator, area_fraction)
# end

# init(::AbstractComponentModel, timestepping, domain) = nothing
# ####

# struct _CplFieldInfo
#     # coupler-specific name
#     name::Val
#     # where the coupler reads/stores info
#     coupler_path::Function
#     # writer model
#     writer_model::Union{AbstractComponentModel, Vector}
#     # writer path
#     writer_path::Function
#     # reader model
#     reader_model::Union{AbstractComponentModel, Vector}
#     # reader model
#     reader_path::Function
#     # # if using an exchange grid
#     # exchange_grid::Bool
#     # map to use for regridding or Nothing (if not needed)
#     regrid_map::Symbol
# end

# mutable struct __CoupledSimulation{FT}
#     comms_context # communication (MPI/GPU) context
#     boundary_space
#     parsed_args
#     model_sims
#     exchange_masks # (on coupler grid)
#     exchange_field_info
#     exchange_field_data
#     clock
#     calendar
#     diagnostics # callback
#     conservation_checks # callback
#     flux_calculator
#     regrid_maps
#     mode_specifics # to revamp with BCReader
#     # callbacks # TBD
#     # running
#     # initialized
# end
# function _CoupledSimulation(comms_context, boundary_space; parsed_args = (;), model_sims=(;), coupler_field_names = (;), clock = (;), calendar = (;), diagnostics = (;), conservation_checks = (;), flux_calculator = (;), regrid_maps = (;), exchange_masks = (;), mode_specifics = (;))

#     # set up coupler fields
#     exchange_field_info= NamedTuple{coupler_field_names}(map(i -> CouplerFieldInfo(Val(i)), coupler_field_names))
#     exchange_field_data = (;)

#     # set the coupled simulation
#     cs = __CoupledSimulation{FT}(comms_context, boundary_space, parsed_args, model_sims, exchange_masks, exchange_field_info, exchange_field_data, clock, calendar, diagnostics, conservation_checks, flux_calculator, regrid_maps, mode_specifics)

#     return cs
# end


# # function coupler_init!(cs)

# #     for i in cs.exchange_field_info
# #         name = (Symbol(string(i.name)[6:end-3]), ) # TODO
# #         data = cs.exchange_field_data

# #         # build cached fields
# #         try i.coupler_path(cs)
# #             nothing
# #         catch
# #             merge(cs.exchange_field_data, NamedTuple{name}((Fields.zeros(cs.boundary_space),))) # move this to struct spec
# #         end
# #         # build libary of regrid maps
# #         if !(i.regrid_map in propertynames(cs.regrid_maps))
# #             push!(cs.regrid_maps, get_map(cs, i.regrid_map))
# #         end

# #         # initialize component models exchange fields
# #         field_pull!(cs, i)

# #         field_push!(cs, i)


# #     end
#     # cs.field_data = NamedTuple{coupler_field_names}(ntuple(i -> Fields.zeros(cs.boundary_space), length(coupler_field_names)))

# # end



# function coupler_callbacks()
#     accumulate_flux
# end





# #_________________________________________
# function field_push!(cs, i::_CouplerFieldInfo)
#     writer_path = i.get_reader_path(cs)
#     coupler_path = i.get_writer_path(cs)

#     if isempty(i.regrid_map) && size(reader_path) == 1
#         eval(:($coupler_path = $reader_path))
#     elseif size(reader_path) > 1
#         combine_surfaces!(reader_path...)
#     elseif !isempty(i.regrid_map)
#         regrid!(cs, getproperty(cs.regrid_maps, name), comb_field)
#     end

# end



# function coupler_pull!(cs)
#     combine()

# end


# _________
# """
#     field_get(i)

# Retrieve data array from coupler and save it in targe`target_field`.
# """
# # function field_put(i::_CouplerFieldInfo)
# function field_put(cs, i) # read source
#     i.get_reader_path(cs)
# end

# function field_get!(cs, i::_CouplerFieldInfo)
#     writer_path = i.get_writer_path(cs)
#     :($writer_path .= ) = deepcopy(i.reader_path)
# end

# function coupler_put(cs)
#     for i in cs.exchange_field_info
#         if tag in i.reader_model
#             field_get(i) # this includes regridding of multiple sfc (do after all sfc updated)
#         end
#     end
# end
# function coupler_get(cs, tag)
#     for i in cs.exchange_field_info
#         if tag in i.writer_model
#             field_put!(i)
#         end
#     end
# end


"""
    Clock{T}

Manages a simulation's time information.
"""
mutable struct _Clock{T}
    t_start::T
    t_end::T    # simulation end time
    t::T         # current simulation time
    dt::T           # simulation timestep
end

mutable struct _Calendar{T}
    start_date::T
    current_date::T
    first_day_of_month::T
end

tick!(clock::Clock) = (clock.time += clock.dt)

# stop_time_exceeded(clock::Clock) = (clock.time >= clock.stop_time)





#___












# helpers




# function coupler_add_map!(coupler::CouplerState, map_name::Symbol, map::Operators.LinearRemap)
#     function coupler_add_map!(coupler::CouplerState, map::Operators.LinearRemap)
#         push!(coupler.remap_operators, map_name => map)
#         push!(coupler.remap_operators, (map.target, map.source) => map)
#     end
#     end








#     function get_combined!()
#         function combine_surfaces(cs)
#             surfaces_fields, surfaces_masks = cs.f
#         end
#     end



#     function get_temperature()
#         temperature(cs) = cs.atmos_sim.u.temp
#     end


#     writer_path = cs.model_sims.col.integrator.u.T
#     reader_path = cs.model_sims.slab.integrator.u.T


# test

# function get_reader_path(cs)
#     cs.t1 * 5
# end

# function get_writer_path(cs)
#     cs.t2 .= cs.cf
# end

# cfi = (;path_read = get_reader_path, path_write = get_writer_path)



# function field_put(i)
#     coupler_path = i.coupler_path
#     reader_path = i.reader_path
#     eval(:($coupler_path .= $reader_path))
# end
