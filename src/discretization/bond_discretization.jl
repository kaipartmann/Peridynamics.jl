struct Bond
    neighbor::Int
    length::Float64
    fail_permit::Bool
end

struct BondDiscretization <: AbstractDiscretization
    position::Matrix{Float64}
    volume::Vector{Float64}
    bonds::Vector{Bond}
    n_neighbors::Vector{Int}
    bond_ids::Vector{UnitRange{Int}}
end

function init_bond_discretization(body::Body, pd::PointDecomposition, chunk_id::Int)
    check_bond_discretization_compat(body.mat)

    bonds, n_neighbors = find_bonds(body, pd.decomp[chunk_id])

    ch = ChunkHandler(bonds, pd, chunk_id)
    localize!(bonds, ch.localizer)

    bond_ids = find_bond_ids(n_neighbors)
    position, volume = get_pos_and_vol_chunk(body, ch.point_ids)
    bd = BondDiscretization(position, volume, bonds, n_neighbors, bond_ids)

    return bd, ch
end

function check_bond_discretization_compat(mat::M) where {M<:AbstractMaterial}
    if discretization_type(mat) !== BondDiscretization
        msg = "body with $(M) incompatible to BondDiscretization!\n"
        msg *= "Check the method `discretization_type` for $(M)!\n"
        throw(ArgumentError(msg))
    end
    return nothing
end

function find_bonds(body::Body, loc_points::UnitRange{Int})
    balltree = BallTree(body.position)
    bonds = Vector{Bond}()
    sizehint!(bonds, body.n_points * 300)
    n_neighbors = zeros(Int, length(loc_points))
    for (li, i) in enumerate(loc_points)
        n_neighbors[li] = find_bonds!(bonds, balltree, body.position, body.fail_permit,
                                      get_point_param(body, :δ, i), i)
    end
    return bonds, n_neighbors
end

function find_bonds!(bonds::Vector{Bond}, balltree::BallTree, position::Matrix{Float64},
                     fail_permit::Vector{Bool}, δ::Float64, i::Int)
    neigh_idxs_with_i = inrange(balltree, view(position, :, i), δ)
    neigh_idxs = filter(x -> x != i, neigh_idxs_with_i)
    for j in neigh_idxs
        L = sqrt((position[1, j] - position[1, i])^2 +
                 (position[2, j] - position[2, i])^2 +
                 (position[3, j] - position[3, i])^2)
        push!(bonds, Bond(j, L, fail_permit[i] & fail_permit[j]))
    end
    return length(neigh_idxs)
end

function find_bond_ids(n_neighbors::Vector{Int})
    bond_ids = fill(0:0, length(n_neighbors))
    bonds_start, bonds_end = 1, 0
    for i in eachindex(n_neighbors)
        bonds_end = bonds_start + n_neighbors[i] - 1
        bond_ids[i] = bonds_start:bonds_end
        bonds_start += n_neighbors[i]
    end
    return bond_ids
end

function get_pos_and_vol_chunk(body::Body, point_ids::AbstractVector{<:Integer})
    position = @views body.position[:, point_ids]
    volume = @views body.volume[point_ids]
    return position, volume
end

@inline each_bond_idx(bd::BondDiscretization, point_id::Int) = bd.bond_ids[point_id]
