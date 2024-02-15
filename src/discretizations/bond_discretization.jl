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

function find_bonds!(bond_buf::Vector{Bond}, position::Matrix{Float64},
                     possible_neighbors::Vector{Int}, fail_permit::Vector{Bool},
                     δ::Float64, i::Int)
    n_neighbors = 0
    for j in possible_neighbors
        if i != j
            L = sqrt((position[1, j] - position[1, i])^2 +
                     (position[2, j] - position[2, i])^2 +
                     (position[3, j] - position[3, i])^2)
            if L ≤ δ
                n_neighbors += 1
                push!(bond_buf, Bond(j, L, fail_permit[i] & fail_permit[j]))
            end
        end
    end
    return n_neighbors
end
