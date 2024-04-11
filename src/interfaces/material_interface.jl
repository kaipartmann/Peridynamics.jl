
function point_param_type(mat::AbstractMaterial)
    throw(MethodError(point_param_type, mat))
end

function allowed_material_kwargs(mat::AbstractMaterial)
    throw(MethodError(allowed_material_kwargs, mat))
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

function halo_read_fields(s::AbstractStorage)
    throw(MethodError(halo_read_fields, s))
end

function halo_write_fields(s::AbstractStorage)
    throw(MethodError(halo_write_fields, s))
end

function get_system(b::AbstractBody, args...)
    throw(MethodError(get_system, b))
end

function is_halo_field(s::AbstractStorage, v::Val{F}) where {F}
    throw(MethodError(is_halo_field, s, v))
end

##----
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
