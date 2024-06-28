struct SingleParamChunk <: AbstractParamSpec end

struct MultiParamChunk <: AbstractParamSpec end

struct ParameterHandler{P<:AbstractPointParameters} <: AbstractParameterHandler
    parameters::Vector{P}
    point_mapping::Vector{Int}
end

function ParameterHandler(body::AbstractBody, ch::AbstractChunkHandler)
    parameters = body.point_params
    point_mapping = body.params_map[ch.point_ids]
    return ParameterHandler(parameters, point_mapping)
end

function get_paramsetup(body::AbstractBody, ::AbstractChunkHandler, ::SingleParamChunk)
    return first(body.point_params)
end

function get_paramsetup(body::AbstractBody, ch::AbstractChunkHandler, ::MultiParamChunk)
    return ParameterHandler(body, ch)
end

@inline function get_params(paramhandler::ParameterHandler, point_id::Int)
    return paramhandler.parameters[paramhandler.point_mapping[point_id]]
end

@inline function get_params(params::AbstractPointParameters, ::Int)
    return params
end

@inline function parameter_setup_type(::Body{M,P}, ::SingleParamChunk) where {M,P}
    return P
end

@inline function parameter_setup_type(::Body{M,P}, ::MultiParamChunk) where {M,P}
    return ParameterHandler{P}
end

@inline function get_param_spec(body::AbstractBody)
    return length(body.point_params) == 1 ? SingleParamChunk() : MultiParamChunk()
end
