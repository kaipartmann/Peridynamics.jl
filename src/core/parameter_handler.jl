"""
    SingleParamChunk

$(internal_api_warning())

Type for a body chunk of a body with only one material parameter set.
"""
struct SingleParamChunk <: AbstractParamSpec end

"""
    MultiParamChunk

$(internal_api_warning())

Type for a body chunk of a body with multiple material parameter sets.
"""
struct MultiParamChunk <: AbstractParamSpec end

"""
    ParameterHandler

$(internal_api_warning())

A type used to manage multiple point parameters defined for the same body. It is used to assign different point parameters to the points of a body.

# Type Parameters

- `P<:AbstractPointParameters`: Point parameter type.

# Fields

- `parameters::Vector{P}`: All parameter sets defined in the simulation.
- `point_mapping::Vector{Int}`: Vector assigning the related parameter set to each
    material point.
"""
struct ParameterHandler{P<:AbstractPointParameters} <: AbstractParameterHandler
    parameters::Vector{P}
    point_mapping::Vector{Int}
end

function ParameterHandler(body::AbstractBody, ch::AbstractChunkHandler)
    parameters = body.point_params
    point_mapping = body.params_map[ch.point_ids]
    return ParameterHandler(parameters, point_mapping)
end

"""
    get_paramsetup(body::AbstractBody, ::AbstractChunkHandler, ::SingleParamChunk)
    get_paramsetup(body::AbstractBody, ch::AbstractChunkHandler, ::MultiParamChunk)

$(internal_api_warning())

Return the parameters of a [`BodyChunk`](@ref) if only one parameter set is defined for the
corresponding [`Body`](@ref) or the parameter handler if multiple parameter sets are
defined.
"""
function get_paramsetup(body::AbstractBody, ::AbstractChunkHandler, ::SingleParamChunk)
    return first(body.point_params)
end

function get_paramsetup(body::AbstractBody, ch::AbstractChunkHandler, ::MultiParamChunk)
    return ParameterHandler(body, ch)
end

"""
    get_params(paramhandler::ParameterHandler, point_id::Int)
    get_params(params::AbstractPointParameters, ::Int)
    get_params(chunk::BodyChunk, point_id::Int)

$(extension_api_note())

Return parameters of a specific point with index `point_id` of a `Body` with parameters
`params` or parameter handler `paramhandler` or of the body `chunk`.

# Resolving parameters inside a kernel

[`force_density_point!`](@ref) receives the whole parameter setup of the chunk, which is one
set of parameters or a handler for a body where [`material!`](@ref) was called more than
once. A kernel resolves the set of its point once before the loop, and a material that
averages a parameter across a bond reads the set of the neighbor inside it:

```julia
params_i = get_params(paramsetup, i)
for bond_id in each_bond_idx(system, i)
    j = system.bonds[bond_id].neighbor
    params_j = get_params(paramsetup, j)
    ...
end
```

For a body with a single parameter set the read inside the loop is the identity, so
`params_j` does not depend on the loop and the compiler hoists the averaging out of it. The
averaging therefore costs nothing where there is nothing to average.
"""
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
