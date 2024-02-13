
# fallback
function point_param_type(mat::M) where {M<:AbstractMaterial}
    throw(MethodError(point_param_type, mat))
end

# fallback
function allowed_material_kwargs(mat::M) where {M<:AbstractMaterial}
    throw(MethodError(allowed_material_kwargs, mat))
end

# fallback
function get_point_params(mat::M, ::Dict{Symbol,Any}) where {M<:AbstractMaterial}
    throw(MethodError(get_elastic_params, mat))
end
