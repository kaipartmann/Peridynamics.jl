function calc_force_density!(b::AbstractBodyChunk)
    b.store.b_int .= 0
    b.store.n_active_bonds .= 0
    for point_id in each_point_idx(b)
        param = get_param(b, point_id)
        force_density!(b.store, b.dscr, param, point_id)
    end
    return nothing
end
