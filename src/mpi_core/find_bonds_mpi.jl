function find_bonds_mpi(pc::PointCloud, mat::M,
                        loc_points::UnitRange{Int}) where {M <: AbstractMaterial}
    n_bonds_estimate = pc.n_points * 500
    bonds = Vector{Bond}()
    sizehint!(bonds, n_bonds_estimate)
    n_family_members = zeros(Int, length(loc_points))
    li = 0
    for i in loc_points
        li += 1
        n_neighbors = find_bonds!(bonds, pc, mat.δ, i)
        n_family_members[li] = n_neighbors
    end
    return bonds, n_family_members
end
