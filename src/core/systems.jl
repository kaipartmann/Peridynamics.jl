function system_type(::M) where {M<:AbstractMaterial}
    return system_type(M)
end

function system_type(mat::Type{M}) where {M}
    return throw(MethodError(system_type, mat))
end

function get_system(::AbstractBody{M}, ::PointDecomposition, ::Int) where {M}
    msg = "system for material $M not specified!\n"
    return error(msg)
end

function log_system(options::AbstractJobOptions, dh::AbstractDataHandler)
    log_system(system_type(dh), options, dh)
    return nothing
end

function required_point_parameters(mat::AbstractMaterial)
    return required_point_parameters(system_type(mat))
end

function required_point_parameters(sys::Type{Sys}) where {Sys<:AbstractSystem}
    return throw(MethodError(required_point_parameters, sys))
end
