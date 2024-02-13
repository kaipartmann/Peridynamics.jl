
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
    throw(MethodError(get_elastic_params, mat))
end

# fallback
function default_export_fields(mat::AbstractMaterial)
    throw(MethodError(default_export_fields, mat))
end
