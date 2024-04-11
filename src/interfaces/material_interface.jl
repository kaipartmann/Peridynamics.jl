
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


function allowed_material_kwargs(mat::AbstractMaterial)
    return DEFAULT_POINT_KWARGS
end

function get_system(::AbstractBody{M}, args...) where {M}
    msg = "system for material $M not specified!\n"
    return error(msg)
end

macro system(material, system)
    local _system_type = quote
        Peridynamics.system_type(::$(esc(material))) = $(esc(system))
    end
    local _get_system = quote
        function Peridynamics.get_system(body::Peridynamics.AbstractBody{$(esc(material))},
                                         args...)
            return $(esc(system))(body, args...)
        end
    end
    return Expr(:block, _system_type, _get_system)
end


macro pointparams(material, pointparams)
    local _pointparam_type = quote
        Peridynamics.point_param_type(::$(esc(material))) = $(esc(pointparams))
    end
    local _get_pointparams = quote
        function Peridynamics.get_point_params(m::$(esc(material)), p::Dict{Symbol,Any})
            return $(esc(pointparams))(m, p)
        end
    end
    return Expr(:block, _pointparam_type, _get_pointparams)
end
