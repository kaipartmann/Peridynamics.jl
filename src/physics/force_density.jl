function calc_force_density!(b::AbstractBodyChunk)
    _calc_force_density!(b.storage, b.system, b.mat, b.paramsetup, each_point_idx(b))
    return nothing
end

function _calc_force_density!(storage::AbstractStorage, system::AbstractSystem,
                              mat::AbstractMaterial, paramsetup::AbstractParameterHandler,
                              each_point_idx)
    storage.b_int .= 0
    storage.n_active_bonds .= 0
    for point_id in each_point_idx
        params = get_params(paramsetup, point_id)
        force_density_point!(storage, system, mat, params, point_id)
    end
    return nothing
end

function _calc_force_density!(storage::AbstractStorage, system::AbstractSystem,
                              mat::AbstractMaterial,
                              params::AbstractPointParameters, each_point_idx)
    storage.b_int .= 0
    storage.n_active_bonds .= 0
    for point_id in each_point_idx
        force_density_point!(storage, system, mat, params, point_id)
    end
    return nothing
end
