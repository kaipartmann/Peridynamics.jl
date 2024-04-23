function calc_force_density!(b::AbstractBodyChunk)
    _calc_force_density!(b.storage, b.system, b.mat, b.paramhandler, each_point_idx(b))
    return nothing
end

function _calc_force_density!(storage::AbstractStorage, system::AbstractSystem,
                              mat::AbstractMaterial, paramhandler::AbstractParameterHandler,
                              each_point_idx)
    storage.b_int .= 0
    storage.n_active_bonds .= 0
    for point_id in each_point_idx
        param = get_params(paramhandler, point_id)
        force_density_point!(storage, system, mat, param, point_id)
    end
    return nothing
end

function _calc_force_density!(storage::AbstractStorage, system::AbstractSystem,
                              mat::AbstractMaterial,
                              paramhandler::AbstractParameterHandler{1}, each_point_idx)
    storage.b_int .= 0
    storage.n_active_bonds .= 0
    param = get_params(paramhandler)
    for point_id in each_point_idx
        force_density_point!(storage, system, mat, param, point_id)
    end
    return nothing
end
