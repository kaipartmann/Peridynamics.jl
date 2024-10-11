function point_param_type(::M) where {M<:AbstractMaterial}
    return point_param_type(M)
end

function point_param_type(mat::Type{T}) where {T}
    return throw(MethodError(point_param_type, mat))
end

function material_type(::P) where {P<:AbstractPointParameters}
    return material_type(P)
end

function material_type(params::Type{T}) where {T}
    return throw(MethodError(material_type, params))
end

function get_point_params(mat::AbstractMaterial, ::Dict{Symbol,Any})
    return throw(MethodError(get_point_params, mat))
end

macro params(material, params)
    macrocheck_input_material(material)
    macrocheck_input_params(params)
    local _checks = quote
        Peridynamics.typecheck_material($(esc(material)))
        Peridynamics.typecheck_params($(esc(params)))
    end
    local _pointparam_type = quote
        Peridynamics.point_param_type(::Base.Type{<:$(esc(material))}) = $(esc(params))
    end
    local _material_type = quote
        Peridynamics.material_type(::Base.Type{<:$(esc(params))}) = $(esc(material))
    end
    local _get_pointparams = quote
        function Peridynamics.get_point_params(m::$(esc(material)), p::Dict{Symbol,Any})
            return $(esc(params))(m, p)
        end
    end
    return Expr(:block, _checks, _pointparam_type, _material_type, _get_pointparams)
end

function macrocheck_input_material(material)
    material isa Symbol && return nothing
    (material isa Expr && material.head === :.) && return nothing
    return throw(ArgumentError("argument `$material` is not a valid material input!\n"))
end

function macrocheck_input_params(params)
    params isa Symbol && return nothing
    (params isa Expr && params.head === :.) && return nothing
    msg = "argument `$params` is not a valid point parameter input!\n"
    return throw(ArgumentError(msg))
end

function typecheck_material(::Type{Material}) where {Material}
    if !(Material <: AbstractMaterial)
        msg = "$Material is not a valid material type!\n"
        throw(ArgumentError(msg))
    end
    return nothing
end

function typecheck_params(::Type{Param}) where {Param}
    if !(Param <: AbstractPointParameters)
        msg = "$Param is not a valid point parameter type!\n"
        throw(ArgumentError(msg))
    end
    material = material_type(Param)
    parameters = fieldnames(Param)
    for req_param in required_point_parameters(material)
        if !in(req_param, parameters)
            msg = "required parameter $req_param not found in $(Param)!\n"
            error(msg)
        end
    end
    return nothing
end
