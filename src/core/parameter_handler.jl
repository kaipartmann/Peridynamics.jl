struct ParameterHandler{P<:AbstractPointParameters,N} <: AbstractParameterHandler{N}
    parameters::Vector{P}
    point_mapping::Vector{Int}

    function ParameterHandler(body::Body{M,P}, ch::AbstractChunkHandler) where {M,P}
        parameters = body.point_params
        point_mapping = body.params_map[ch.loc_points]
        N = length(body.point_params)
        return new{P,N}(parameters, point_mapping)
    end
end

@inline function get_params(ph::ParameterHandler, point_id::Int)
    return ph.parameters[ph.point_mapping[point_id]]
end

@inline function get_params(ph::ParameterHandler{P,1}, ::Int) where {P}
    return first(ph.parameters)
end

@inline function get_params(ph::ParameterHandler{P,1}) where {P}
    return first(ph.parameters)
end

@inline function parameter_handler_type(body::Body{M,P}) where {M,P}
    N = length(body.point_params)
    return ParameterHandler{P,N}
end
