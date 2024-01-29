struct Bond
    j::Int
    L::Float64
    failure_allowed::Bool
end

struct PointBondDiscretization <: AbstractDiscretization
    position::Matrix{Float64}
    volume::Vector{Float64}
    bonds::Vector{Bond}
    n_neighbors::Vector{Int}
    bond_range::Vector{UnitRange{Int}}
end

function find_bonds!(bond_buf::Vector{Bond}, pc::PointCloud, δ::Float64, i::Int)
    n_neighbors = 0
    for j in axes(pc.position, 2)
        if i != j
            L = sqrt((pc.position[1, j] - pc.position[1, i])^2 +
                     (pc.position[2, j] - pc.position[2, i])^2 +
                     (pc.position[3, j] - pc.position[3, i])^2)
            if L <= δ
                n_neighbors += 1
                failure_allowed = pc.failure_allowed[i] & pc.failure_allowed[j]
                push!(bond_buf, Bond(j, L, failure_allowed))
            end
        end
    end
    return n_neighbors
end

function find_bonds(pc::PointCloud, mat::AbstractMaterial, loc_points::UnitRange{Int})
    bonds = Vector{Bond}()
    sizehint!(bonds, pc.n_points * 300)
    n_neighbors = zeros(Int, length(loc_points))
    for (li, i) in enumerate(loc_points)
        n_neighbors[li] = find_bonds!(bonds, pc, mat.δ, i)
    end
    return bonds, n_neighbors
end

function localize!(bonds::Vector{Bond}, localizer::Dict{Int,Int})
    for i in eachindex(bonds)
        bond = bonds[i]
        bonds[i] = Bond(localizer[bond.j], bond.L, bond.failure_allowed)
    end
    return nothing
end

function find_bond_range(n_neighbors::Vector{Int})
    bond_range = fill(0:0, length(n_neighbors))
    bonds_start, bonds_end = 1, 0
    for i in eachindex(n_neighbors)
        bonds_end = bonds_start + n_neighbors[i] - 1
        bond_range[i] = bonds_start:bonds_end
        bonds_start += n_neighbors[i]
    end
    return bond_range
end

function find_halo_points(bonds::Vector{Bond}, loc_points::UnitRange{Int})
    _halo_points = Set{Int}()
    for bond in bonds
        neighbor = bond.j
        if !in(neighbor, loc_points)
            push!(_halo_points, neighbor)
        end
    end
    halo_points = [i for i in _halo_points]
    return halo_points
end
