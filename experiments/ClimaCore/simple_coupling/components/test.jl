# test.jl
using ClimaCore: Fields

using ClimaCoupler.TestHelper: create_space

using ClimaCoupler.Regridder: combine_surfaces!, nans_to_zero

# coupler regridder
abstract type AbstractRegridMap end
mapname(map::AbstractRegridMap) = Symbol(string(map)[1:end-2])
mapname(::Nothing) = nothing

struct DummyRegridMap <: AbstractRegridMap end

function get_map(::DummyRegridMap)
    function (cs, target, source,)
        parent(target) .= parent(source)
    end
end

function remap!(cs, target, source, map,)
    map(cs, target, source)
end

function apply_mask!(f_t::Fields.Field, f_s::Fields.Field, ::Nothing)
    f_t .= f_s
end
function apply_mask!(f_t_tuple::NamedTuple, f_s::Fields.Field, m_tuple::NamedTuple)
    for (name, f_t) in zip(propertynames(f_t_tuple), f_t_tuple)
        f_t .= f_s .* getproperty(m_tuple, name)
    end
end

# coupler simulation
abstract type AbstractComponentModel end # tag upon which to dispatch

struct ThermalSlab4 <: AbstractComponentModel end
struct DiffusiveColumn4 <: AbstractComponentModel end

abstract type AbstractCouplingType end
struct Sequential <: AbstractCouplingType  end

# utilities
val2sym(v::Val) = Symbol(string(v)[6:end-3])

# exchange field info
struct ExchangeFieldInfo
    # coupler-specific name
    name::Val
    # where the coupler reads/stores info
    coupler_path::Function # TODO change these to printable expressions, use eval
    # source model
    source_models::Tuple
    # source path
    source_path::Function
    # source mask
    source_mask::Function
    # target model
    target_models::Tuple
    # target path
    target_path::Function
    # target mask
    target_mask::Function
    # # if using an exchange grid
    # exchange_grid::Bool
    # map to use for regridding or Nothing (if not needed)
    regrid_source::Union{AbstractRegridMap, Nothing}
    regrid_target::Union{AbstractRegridMap, Nothing}
end

# in experiment adapters
function ExchangeFieldInfo(name::Val{:T_S}) # ; space = false) # option for exchange grid

    name = Val(:T_S) # #todo
    target_models = (DiffusiveColumn4(),)
    source_models = (ThermalSlab4(),) # model that will trigger pull - needs to be the last surface model (TODO: specify `updated`` status to all surfaces instead)

    # set regridding type
    regrid_source = DummyRegridMap() # or nothing
    regrid_target = nothing # or nothing

    # how to obtain source and target
    get_source_path(cs) = (; sfc1 = cs.model_sims.slab.T, sfc2 = cs.model_sims.slab.T)# how coupler views this field
    get_target_path(cs) = cs.model_sims.col.T

    # set source and target masks
    get_source_mask(cs) = (; sfc1 = 0.5 .* ones(cs.boundary_space), sfc2 = 0.5 .* ones(cs.boundary_space))
    get_target_mask(_) = nothing

    get_coupler_path_(cs) = cs.exchange_field_data.T_S # todo: use `name` here # get_target_path # where the coupler stores the variable

    ExchangeFieldInfo(name, get_coupler_path_, source_models, get_source_path, get_source_mask, target_models, get_target_path, get_target_mask, regrid_source, regrid_target)
end

# coupler simulation

struct CouplerSimulation{I,D,R,B,M,L}
    exchange_field_info::I
    exchange_field_data::D
    regrid_maps::R
    boundary_space::B
    model_sims::M
    callbacks::L
end
function CouplerSimulation(exchange_field_info, boundary_space; exchange_field_data = (;), regrid_maps = (;), callbacks = (;), model_sims=(;))

    for i in exchange_field_info
        name = val2sym(i.name)
        # build cached fields (if not pointers)
        if !((i.coupler_path == i.source_path) || (i.coupler_path == i.target_path))
            exchange_field_data = merge(exchange_field_data, NamedTuple{(name,)}((Fields.zeros(boundary_space),))) # move this to struct spec, TODO: consider Fields.FieldVector()
        end
        # build libary of regrid maps
        all_maps = push!([], (i.regrid_source, i.regrid_target)...)
        for map in all_maps
            map_name = mapname(map)
            if !isnothing(map_name) && !(map_name in propertynames(regrid_maps))
                regrid_maps = merge(regrid_maps, NamedTuple{(map_name,)}((get_map(map),))) # TODO push may be cleaner?
            end
        end

    end
    CouplerSimulation(exchange_field_info, exchange_field_data, regrid_maps, boundary_space, model_sims, callbacks)
end

# coupler field exchange
function coupler_reinit(cs)
    # initialize all exchsnge fields
    for i in cs.exchange_field_info
        field_pull!(cs, i)
    end
end

function model_pull(tag::AbstractComponentModel, cs)
    for i in cs.exchange_field_info
        if tag in i.source_models
            field_pull!(cs, i)
        end
    end
end
function field_pull!(cs, i::ExchangeFieldInfo)
    source_path = i.source_path(cs)
    coupler_path = i.coupler_path(cs)

    if isnothing(i.regrid_source) && length(source_path) == 1 # no regridding or combining needed (i.e., same mesh)
        # @assert coupler_space == source_space
        coupler_path .= source_path
    else
        if length(source_path) > 1 # multiple surface models
            combined_field = zeros(cs.boundary_space)
            source_masks = i.source_mask(cs)
            combine_surfaces!(combined_field, source_masks, source_path)
            source_path = combined_field
        end
        if isnothing(i.regrid_source)
            coupler_path .= source_path
        else
            name = mapname(i.regrid_source)
            remap!(cs, coupler_path, source_path, getproperty(cs.regrid_maps, name))
        end
    end
end

function model_push(tag::AbstractComponentModel, cs)
    for i in cs.exchange_field_info
        if tag in (i.target_models)
            field_push!(cs, i)
        end
    end

end
function field_push!(cs, i::ExchangeFieldInfo)
    target_path = i.target_path(cs)
    coupler_path = i.coupler_path(cs)
    target_mask = i.target_mask(cs)

    if isnothing(i.regrid_target) # no regridding or combining needed (i.e., same mesh)
        # @assert coupler_space == source_space
        # eval(:($coupler_path .= getindex($source_path, 1))) # only broadcast
        apply_mask!(target_path, coupler_path, target_mask)
    else
        name = mapname(i.regrid_target)
        combined_field = zeros(cs.boundary_space)
        remap!(cs, combined_field, coupler_path, getproperty(cs.regrid_maps, name))
        apply_mask!(i.target_path(cs), combined_field, target_mask)
    end
end

function coupler_callbacks(cs)
    for cb in cs.callbacks
        cb(cs)
    end
end


# driver
boundary_space = create_space(Float64);

# function accumulate_flux(cs) # this will need to be hard coded
#     cs.accumulated_field_data.F_S += cs.exchange_field_data.F_S
# end

# function calculate_flux(cs) # this will need to be hard coded
#     cs.exchange_field_data.F_S = calcualte_flux(cs)
# end

exchange_field_info = (ExchangeFieldInfo(Val(:T_S)), )
model_sims = (; col = (; tag = DiffusiveColumn4(), T = ones(boundary_space),), slab = (; tag = ThermalSlab4(), T = ones(boundary_space),) );
cs = CouplerSimulation(exchange_field_info, boundary_space, callbacks = (;), model_sims = model_sims);

coupler_reinit(cs)

function step_coupler!(cs, ::Sequential)
    model_tags = values(map(i -> string(getproperty(i, :tag)), cs.model_sims))
    @info "Sequential coupling: \n $model_tags"
    for model in cs.model_sims
        model_push(model.tag, cs)
        # model.step!()
        model_pull(model.tag, cs) #only striger surface pull once
        # coupler_calbacks()
    end
end

# needs to deal with 3 models not 2

step_coupler!(cs, Sequential())