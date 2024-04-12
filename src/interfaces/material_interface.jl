
function point_param_type(mat::AbstractMaterial)
    throw(MethodError(point_param_type, mat))
end

function get_point_params(mat::AbstractMaterial, ::Dict{Symbol,Any})
    throw(MethodError(get_point_params, mat))
end

function system_type(mat::AbstractMaterial)
    throw(MethodError(system_type, mat))
end

function storage_type(mat::AbstractMaterial, ts::AbstractTimeSolver)
    throw(MethodError(storage_type, mat, ts))
end

function allowed_material_kwargs(::AbstractMaterial)
    return DEFAULT_POINT_KWARGS
end

function get_system(::AbstractBody{M}, args...) where {M}
    msg = "system for material $M not specified!\n"
    return error(msg)
end

macro system(material, system)
    macrocheck_input_material(material)
    macrocheck_input_system(system)
    local _checks = quote
        typecheck_material($(esc(material)))
        typecheck_system($(esc(system)))
    end
    local _system_type = quote
        Peridynamics.system_type(::$(esc(material))) = $(esc(system))
    end
    local _get_system = quote
        function Peridynamics.get_system(body::Peridynamics.AbstractBody{$(esc(material))},
                                         args...)
            return $(esc(system))(body, args...)
        end
    end
    return Expr(:block, _checks, _system_type, _get_system)
end

macro params(material, params)
    macrocheck_input_material(material)
    macrocheck_input_params(params)
    local _checks = quote
        typecheck_material($(esc(material)))
        typecheck_params($(esc(params)))
    end
    local _pointparam_type = quote
        Peridynamics.point_param_type(::$(esc(material))) = $(esc(params))
    end
    local _get_pointparams = quote
        function Peridynamics.get_point_params(m::$(esc(material)), p::Dict{Symbol,Any})
            return $(esc(params))(m, p)
        end
    end
    return Expr(:block, _checks, _pointparam_type, _get_pointparams)
end

function macrocheck_input_material(material)
    material isa Symbol && return nothing
    (material isa Expr && material.head === :.) && return nothing
    return throw(ArgumentError("argument `$material` is not a valid material input!\n"))
end

function macrocheck_input_system(system)
    system isa Symbol && return nothing
    (system isa Expr && system.head === :.) && return nothing
    return throw(ArgumentError("argument `$system` is not a valid system input!\n"))
end

function macrocheck_input_params(params)
    params isa Symbol && return nothing
    (params isa Expr && params.head === :.) && return nothing
    msg = "argument `$params` is not a valid point parameter input!\n"
    return throw(ArgumentError(msg))
end

function typecheck_system(::Type{System}) where {System}
    if !(System <: AbstractSystem)
        msg = "$System is not a valid system type!\n"
        throw(ArgumentError(msg))
    end
    return nothing
end
