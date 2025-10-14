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

@inline function get_n_points(system::AbstractSystem)
    return get_n_points(system.chunk_handler)
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

@inline function each_halo_idx(system::AbstractSystem)
    return each_halo_idx(system.chunk_handler)
end

@inline function each_point_idx_pair(system::AbstractSystem)
    return each_point_idx_pair(system.chunk_handler)
end

@inline function get_loc_view(a::AbstractArray, system::AbstractSystem)
    return get_loc_view(a, system.chunk_handler)
end

@inline function get_n_dim(::AbstractSystem)
    return 3 # 3D only for now
end

@inline function get_n_dof(system::AbstractSystem)
    return get_n_dim(system) * get_n_points(system)
end

@inline function get_n_loc_dof(system::AbstractSystem)
    return get_n_dim(system) * get_n_loc_points(system)
end

@inline function each_dim(system::AbstractSystem)
    return 1:get_n_dim(system)
end

@inline function get_dof(system::AbstractSystem, dim::Int, point_id::Int)
    return get_dof(get_n_dim(system), dim, point_id)
end

@inline function get_dof(n_dim::Int, dim::Int, point_id::Int)
    return (point_id - 1) * n_dim + dim
end

@inline function each_dof_idx(system::AbstractSystem)
    return each_dof_idx(get_n_dim(system), 1:get_n_points(system))
end

@inline function each_loc_dof_idx(system::AbstractSystem)
    return each_dof_idx(get_n_dim(system), 1:get_n_loc_points(system))
end

@inline function each_dof_idx(system::AbstractSystem, idxs::AbstractVector{<:Integer})
    return each_dof_idx(get_n_dim(system), idxs)
end

#=
This function generates a cartesian product of dof indices, dimensions, and point indices.
In Julia it is very convenient, because elements in a 2-dimensional array can be
accessed via:
    A[dof]
or via
    A[dim, i]
=#
@inline function each_dof_idx(n_dim::Int, idxs::AbstractVector{<:Integer})
    return ((get_dof(n_dim, dim, i), dim, i) for i in idxs, dim in 1:n_dim)
end

@inline function each_dof(system::AbstractSystem)
    return each_dof(get_n_dim(system), 1:get_n_points(system))
end

@inline function each_loc_dof(system::AbstractSystem)
    return each_dof(get_n_dim(system), 1:get_n_loc_points(system))
end

@inline function each_dof(system::AbstractSystem, idxs::AbstractVector{<:Integer})
    return each_dof(get_n_dim(system), idxs)
end

@inline function each_dof(n_dim::Int, idxs::AbstractVector{<:Integer})
    return (get_dof(n_dim, dim, i) for i in idxs, dim in 1:n_dim)
end

# Get the point index from a dof index.
get_point(system::AbstractSystem, idx::Int) = get_point(get_n_dim(system), idx)
get_point(n_dim::Int, idx::Int) = div(idx - 1, n_dim) + 1

# Get the dimension index from a dof index.
get_dim(system::AbstractSystem, idx::Int) = get_dim(get_n_dim(system), idx)
get_dim(n_dim::Int, idx::Int) = mod(idx - 1, n_dim) + 1

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
