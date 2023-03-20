# test.jl
using ClimaCore: Fields

using ClimaCoupler.TestHelper: create_space

using ClimaCoupler.Regridder: combine_surfaces!, nans_to_zero

# module: coupler simulation
abstract type AbstractComponentModel end # tag upon which to dispatch

abstract type AbstractSimulation end  # Simulation struct containing all sim info

struct ComponentSimulation{T, Y, I, A} <: AbstractSimulation
    tag::T
    Y_init::Y
    integrator::I
    area_fraction::A
end
ComponentSimulation(tag::AbstractComponentModel; Y_init = (;), integrator = (;), area_fraction = nothing) = ComponentSimulation(tag, Y_init, integrator, area_fraction)

struct CouplerSimulation{I,D,R,B,M,L} <: AbstractSimulation
    exchange_field_info::I
    exchange_field_data::D
    regrid_maps::R
    boundary_space::B
    model_sims::M
    callbacks::L
end
function CouplerSimulation(exchange_field_info, boundary_space; exchange_field_data = (;), regrid_maps = (;), callbacks = (;), model_sims=(;))

    for i in exchange_field_info
        name = i.name
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

abstract type AbstractCouplingType end
struct Sequential <: AbstractCouplingType  end

# module: coupler regridder
abstract type AbstractRegridMap end
mapname(map::AbstractRegridMap) = Symbol(string(map)[1:end-2])
mapname(::Nothing) = nothing

struct DummyRegridMap <: AbstractRegridMap end

function get_map(::DummyRegridMap)
    function (cs::CouplerSimulation, target, source,)
        parent(target) .= parent(source)
    end
end

function remap!(cs::CouplerSimulation, target, source, map,)
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

# module: callbacks
function coupler_callbacks(cs::CouplerSimulation)
    for cb in cs.callbacks
        cb(cs)
    end
end

# module: utilities
val2sym(v::Val) = Symbol(string(v)[6:end-3])

# module: model exchange
struct ExchangeFieldInfo
    # coupler-specific name
    name::Symbol
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

function coupler_reinit(cs::CouplerSimulation)
    # initialize all exchsnge fields
    for i in cs.exchange_field_info
        field_pull!(cs, i)
    end
end

function model_pull(tag::AbstractComponentModel, cs::CouplerSimulation)
    for i in cs.exchange_field_info
        if tag in i.source_models
            field_pull!(cs, i)
        end
    end
end
function field_pull!(cs::CouplerSimulation, i::ExchangeFieldInfo)
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

function model_push(tag::AbstractComponentModel, cs::CouplerSimulation)
    for i in cs.exchange_field_info
        if tag in (i.target_models)
            field_push!(cs, i)
        end
    end

end
function field_push!(cs::CouplerSimulation, i::ExchangeFieldInfo)
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

# user spec: experiment/../user_field_exchange.jl
function ExchangeFieldInfo(name::Val{:T_S}) # ; space = false) # option for exchange grid

    name = val2sym(name) # #todo

    # define which models trigger the pull and push exchanges
    target_models = (DiffusiveColumn(),)
    source_models = (ThermalSlab2(),) # only list the model(s) that trigger the exchange (TODO: specify `updated`` status to all surfaces instead)

    # set regridding type
    regrid_source = DummyRegridMap()
    regrid_target = nothing

    # how to obtain source and target
    get_source_path(cs::CouplerSimulation) = (; f1 = cs.model_sims.slab1.integrator.T, f2 = cs.model_sims.slab2.integrator.T) # how coupler views this field; names must correspond to those of masks below TODO: automate
    get_target_path(cs::CouplerSimulation) = cs.model_sims.col.integrator.T

    # set source and target masks
    get_source_mask(cs::CouplerSimulation) = (; f1 = 0.5 .* ones(cs.boundary_space), f2 = 0.5 .* ones(cs.boundary_space))
    get_target_mask(::CouplerSimulation) = nothing

    get_coupler_path_(cs::CouplerSimulation) = cs.exchange_field_data.T_S # todo: use `name` here # get_target_path # where the coupler stores the variable

    ExchangeFieldInfo(name, get_coupler_path_, source_models, get_source_path, get_source_mask, target_models, get_target_path, get_target_mask, regrid_source, regrid_target)
end

# user spec: experiment/../model_init.jl
struct ThermalSlab <: AbstractComponentModel end
struct ThermalSlab2 <: AbstractComponentModel end
struct DiffusiveColumn <: AbstractComponentModel end

# driver
boundary_space = create_space(Float64);

# function accumulate_flux(cs::CouplerSimulation) # this will need to be hard coded
#     cs.accumulated_field_data.F_S += cs.exchange_field_data.F_S
# end

# function calculate_flux(cs::CouplerSimulation) # this will need to be hard coded
#     cs.exchange_field_data.F_S = calcualte_flux(cs)
# end

exchange_field_info = (ExchangeFieldInfo(Val(:T_S)), )

col = ComponentSimulation(DiffusiveColumn(), integrator = (T = ones(boundary_space),), )
slab1 = ComponentSimulation(ThermalSlab(), integrator = (T = zeros(boundary_space),), )
slab2 = ComponentSimulation(ThermalSlab2(), integrator = (T = ones(boundary_space),), )

model_sims = (; col = col, slab1 = slab1, slab2 = slab2);
cs = CouplerSimulation(exchange_field_info, boundary_space, callbacks = (;), model_sims = model_sims);

coupler_reinit(cs)

function step_coupler!(cs::CouplerSimulation, ::Sequential)
    model_tags = values(map(i -> string(getproperty(i, :tag)), cs.model_sims))
    @info "Sequential coupling: \n $model_tags"
    for model in cs.model_sims
        model_push(model.tag, cs)
        # model.step!()
        model_pull(model.tag, cs)
        # coupler_calbacks()
    end
end

step_coupler!(cs, Sequential())

using Test
@test extrema(cs.model_sims.slab1.integrator.T) == (0.0, 0.0)
@test extrema(cs.model_sims.slab2.integrator.T) == (1.0, 1.0)
@test extrema(cs.model_sims.col.integrator.T) == (0.5, 0.5)
