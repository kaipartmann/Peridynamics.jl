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
    Threads.@threads :static for point_id in each_point_idx
        force_density_point!(storage, system, mat, params, point_id)
    end
    return nothing
end

function calc_force_density!(b::AbstractBodyChunk, nhs)
    _calc_force_density!(b.storage, b.system, b.mat, b.paramsetup, each_point_idx(b), nhs)
    return nothing
end

function calc_force_density!(b::AbstractBodyChunk, nhs::Nothing)
    # Call original function
    _calc_force_density!(b.storage, b.system, b.mat, b.paramsetup, each_point_idx(b))
    return nothing
end

function _calc_force_density!(storage::AbstractStorage, system::AbstractSystem,
                              mat::AbstractMaterial,
                              params::AbstractPointParameters, each_point_idx, nhs)
    storage.b_int .= 0
    storage.n_active_bonds .= 0

    # for_particle_neighbor(storage.position, storage.position, nhs,
    #                       particles=each_point_idx, parallel=false) do i, j, _, L
    Threads.@threads :static for point_id in each_point_idx
    foreach_neighbor(storage.position, storage.position, nhs, point_id) do i, j, _, L
        # current bond length
        Δxijx = storage.position[1, j] - storage.position[1, i]
        Δxijy = storage.position[2, j] - storage.position[2, i]
        Δxijz = storage.position[3, j] - storage.position[3, i]
        l = sqrt(Δxijx * Δxijx + Δxijy * Δxijy + Δxijz * Δxijz)

        # bond strain
        ε = (l - L) / L

        # failure mechanism
        # if ε > params.εc && bond.fail_permit
        #     storage.bond_active[bond_id] = false
        # end
        # storage.n_active_bonds[i] += storage.bond_active[bond_id]

        # update of force density
        scfactor = 1
        bond_fail = true
        temp = bond_fail * scfactor * params.bc * ε / l * system.volume[j]
        storage.b_int[1, i] += temp * Δxijx
        storage.b_int[2, i] += temp * Δxijy
        storage.b_int[3, i] += temp * Δxijz
    end
    end
    return nothing
end
