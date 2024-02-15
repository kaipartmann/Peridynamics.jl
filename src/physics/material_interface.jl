
# fallback
function point_param_type(mat::AbstractMaterial)
    throw(MethodError(point_param_type, mat))
end

# fallback
function allowed_material_kwargs(mat::AbstractMaterial)
    throw(MethodError(allowed_material_kwargs, mat))
end

# fallback
function get_point_params(mat::AbstractMaterial, ::Dict{Symbol,Any})
    throw(MethodError(get_point_params, mat))
end

# fallback
@inline function default_export_fields(::Type{M}) where {M<:AbstractMaterial}
    return (:displacement, :damage)
end
