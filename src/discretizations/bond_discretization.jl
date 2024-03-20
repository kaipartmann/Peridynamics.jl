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
    gnhs = GridNeighborhoodSearch{3}(get_point_param(body, :δ, 1), body.n_points)
    TrixiParticles.initialize!(gnhs, body.position)
    bonds = Vector{Bond}()
    sizehint!(bonds, body.n_points * 300)

    n_neighbors = zeros(Int, body.n_points)
    for_particle_neighbor(@views(body.position[:,loc_points]), body.position, gnhs, parallel=false) do i, j, _, L
        i == j && return nothing
        push!(bonds, Bond(j, L, body.fail_permit[i] & body.fail_permit[j]))
        n_neighbors[i] += 1
        return nothing
    end
    return bonds, n_neighbors
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
