
struct PointSetsPreCrack <: AbstractPredefinedCrack
    set_a::Symbol
    set_b::Symbol
end

function apply_precracks!(b::AbstractBodyChunk, body::Body)
    for precrack in body.point_sets_precracks
        apply_precrack!(b, body, precrack)
    end
    calc_damage!(b)
    return nothing
end

function apply_precrack!(b::AbstractBodyChunk, body::Body, precrack::PointSetsPreCrack)
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

function _apply_precrack!(s::AbstractStorage, bd::BondSystem, ch::ChunkHandler,
                          set_a::Vector{Int}, set_b::Vector{Int})
    s.n_active_bonds .= 0
    for point_id in each_point_idx(ch)
        for bond_id in each_bond_idx(bd, point_id)
            bond = bd.bonds[bond_id]
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
