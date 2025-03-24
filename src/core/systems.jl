function system_type(mat::AbstractMaterial)
    return throw(InterfaceError(mat, "system_type"))
end

function get_system(::AbstractBody{M}, ::PointDecomposition, ::Int) where {M}
    msg = "system for material $M not specified!\n"
    return error(msg)
end

function log_system(options::AbstractJobOptions, dh::AbstractDataHandler)
    log_system(system_type(dh), options, dh)
    return nothing
end

@inline function get_n_loc_points(system::AbstractSystem)
    return get_n_loc_points(system.chunk_handler)
end

@inline function get_point_ids(system::AbstractSystem)
    return get_point_ids(system.chunk_handler)
end

@inline function get_loc_points(system::AbstractSystem)
    return get_loc_points(system.chunk_handler)
end

@inline function get_halo_points(system::AbstractSystem)
    return get_halo_points(system.chunk_handler)
end

@inline function get_hidxs_by_src(system::AbstractSystem)
    return get_hidxs_by_src(system.chunk_handler)
end

@inline function get_localizer(system::AbstractSystem)
    return get_localizer(system.chunk_handler)
end

@inline function each_point_idx(system::AbstractSystem)
    return each_point_idx(system.chunk_handler)
end

@inline function each_point_idx_pair(system::AbstractSystem)
    return each_point_idx_pair(system.chunk_handler)
end

@inline function get_loc_view(a::AbstractArray, system::AbstractSystem)
    return get_loc_view(a, system.chunk_handler)
end

@inline function get_n_points(system::AbstractSystem)
    return get_n_points(system.chunk_handler)
end

@inline function init_field_system(system, field)
    return nothing
end

function log_material(mat::M; indentation::Int=2) where {M}
    msg = msg_qty("material type", nameof(M); indentation)
    for prop in fieldnames(M)
        msg *= log_material_property(Val(prop), mat; indentation)
    end
    return msg
end

function log_material_property(::Val{S}, mat; indentation) where {S}
    return ""
end
