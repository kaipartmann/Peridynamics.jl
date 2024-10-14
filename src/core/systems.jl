function system_type(mat::AbstractMaterial)
    return system_type(typeof(mat))
end

function system_type(mat::Type{M}) where {M<:AbstractMaterial}
    return throw(InterfaceError(mat, system_type))
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

function required_point_parameters(system::Type{<:AbstractSystem})
    return throw(InterfaceError(system, required_point_parameters))
end
