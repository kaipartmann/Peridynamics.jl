struct PointSetsPreCrack <: AbstractPredefinedCrack
    set_a::Symbol
    set_b::Symbol
    filter_bonds::Bool
end

@inline filter_bonds(crack::AbstractPredefinedCrack) = crack.filter_bonds

function apply_precracks!(b::AbstractBodyChunk, body::AbstractBody)
    for precrack in body.point_sets_precracks
        apply_precrack!(b, body, precrack)
    end
    calc_damage!(b)
    return nothing
end

function apply_precrack!(b::AbstractBodyChunk, body::AbstractBody,
                         crack::PointSetsPreCrack)
    filter_bonds(crack) && return nothing
    set_a = filter(x -> in(x, b.ch.point_ids), body.point_sets[crack.set_a])
    set_b = filter(x -> in(x, b.ch.point_ids), body.point_sets[crack.set_b])
    localize!(set_a, b.ch.localizer)
    localize!(set_b, b.ch.localizer)
    if isempty(set_a) || isempty(set_b)
        return nothing
    end
    break_bonds!(b.storage, b.system, b.ch, set_a, set_b)
    return nothing
end

"""
    precrack!(b::AbstractBody, set_a::Symbol, set_b::Symbol; update_dmg::Bool=true)

creates a crack between two point sets by prohibiting interaction between points of
different point sets

# Arguments

- `b::AbstractBody`: peridynamic body
- `set_a::Symbol`: first point set
- `set_b::Symbol`: second point set

# Keywords

- `update_dmg::Bool`: if `update_dmg=true`, the material points involved in the predefined
                      are initially damaged. If `update_dmg=false`, the bonds involved are
                      deleted and the material points involved with the predefined crack
                      are not damaged.
                      (default=`true`)

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
function precrack!(b::AbstractBody, set_a::Symbol, set_b::Symbol; update_dmg::Bool=true)
    check_if_set_is_defined(b.point_sets, set_a)
    check_if_set_is_defined(b.point_sets, set_b)
    check_if_sets_intersect(b.point_sets, set_a, set_b)
    push!(b.point_sets_precracks, PointSetsPreCrack(set_a, set_b, !update_dmg))
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
