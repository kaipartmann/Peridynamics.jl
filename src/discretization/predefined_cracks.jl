struct PointSetsPreCrack <: AbstractPredefinedCrack
    set_a::Symbol
    set_b::Symbol
end

function apply_precracks!(b::AbstractBodyChunk, body::AbstractBody)
    for precrack in body.point_sets_precracks
        apply_precrack!(b, body, precrack)
    end
    calc_damage!(b)
    return nothing
end

function apply_precrack!(b::AbstractBodyChunk, body::AbstractBody,
                         precrack::PointSetsPreCrack)
    set_a = filter(x -> in(x, b.ch.point_ids), body.point_sets[precrack.set_a])
    set_b = filter(x -> in(x, b.ch.point_ids), body.point_sets[precrack.set_b])
    localize!(set_a, b.ch.localizer)
    localize!(set_b, b.ch.localizer)
    if isempty(set_a) || isempty(set_b)
        return nothing
    end
    _apply_precrack!(b.storage, b.system, b.ch, set_a, set_b)
    return nothing
end

function _apply_precrack!(s::AbstractStorage, system::AbstractSystem,
                          ch::AbstractChunkHandler, set_a::Vector{Int}, set_b::Vector{Int})
    s.n_active_bonds .= 0
    for point_id in each_point_idx(ch)
        for bond_id in each_bond_idx(system, point_id)
            bond = system.bonds[bond_id]
            neighbor_id = bond.neighbor
            point_in_a = in(point_id, set_a)
            point_in_b = in(point_id, set_b)
            neigh_in_a = in(neighbor_id, set_a)
            neigh_in_b = in(neighbor_id, set_b)
            if (point_in_a && neigh_in_b) || (point_in_b && neigh_in_a)
                s.bond_active[bond_id] = false
            end
            s.n_active_bonds[point_id] += s.bond_active[bond_id]
        end
    end
    return nothing
end

function check_if_sets_intersect(point_sets::Dict{Symbol,Vector{Int}}, key_a::Symbol,
                                 key_b::Symbol)
    set_a, set_b = point_sets[key_a], point_sets[key_b]
    if !isempty(set_a ∩ set_b)
        msg = "set :$key_a and :$key_b intersect!\n"
        msg *= "No point of the first set is allowed in the second!\n"
        throw(ArgumentError(msg))
    end
    return nothing
end

"""
    precrack!(b::AbstractBody, set_a::Symbol, set_b::Symbol)

creates a crack between two point sets by prohibiting interaction between points of
different point sets

# Arguments

- `b::AbstractBody`: peridynamic body
- `set_a::Symbol`: first point set
- `set_b::Symbol`: second point set

# Throws

- error if point set `set_a` or `set_b` does not exist
- error if point sets contain common points

# Example

```julia-repl
julia> precrack!(b, :set_a, :set_b)
julia> b.point_sets_precracks
1-element Vector{Peridynamics.PointSetsPreCrack}:
 Peridynamics.PointSetsPreCrack(:set_a, :set_b)
```
"""
function precrack!(b::AbstractBody, set_a::Symbol, set_b::Symbol)
    check_if_set_is_defined(b.point_sets, set_a)
    check_if_set_is_defined(b.point_sets, set_b)
    check_if_sets_intersect(b.point_sets, set_a, set_b)
    push!(b.point_sets_precracks, PointSetsPreCrack(set_a, set_b))
    return nothing
end
