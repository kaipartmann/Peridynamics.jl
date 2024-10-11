function system_type(mat::AbstractMaterial)
    throw(MethodError(system_type, mat))
end

function get_system(::AbstractBody{M}, ::PointDecomposition, ::Int) where {M}
    msg = "system for material $M not specified!\n"
    return error(msg)
end

function log_system(options::AbstractJobOptions, dh::AbstractDataHandler)
    log_system(system_type(dh), options, dh)
    return nothing
end

function required_point_parameters(sys::AbstractSystem)
    return throw(MethodError(required_point_parameters, sys))
end
