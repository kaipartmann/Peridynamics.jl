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
    precrack!(body, set_a, set_b; update_dmg=true)

Creates a crack between two point sets by prohibiting interaction between points of
different point sets. The points in `set_a` are not allowed to interact with points in
`set_b`.

# Arguments

- `body::AbstractBody`: [`Body`](@ref).
- `set_a::Symbol`: The name of a point set of this body.
- `set_b::Symbol`: The name of a point set of this body.

# Keywords

- `update_dmg::Bool`: If `true`, the material points involved in the predefined crack are
    initially damaged. If `false`, the bonds involved are deleted and the material points
    involved with the predefined crack are not damaged in the reference results.
    (default: `true`)

# Throws

- Errors if the body does not contain sets with name `set_a` and `set_b`.
- Errors if the point sets intersect and a point is included in both sets.

# Example

```julia-repl
julia> point_set!(body, :a, 1:2)

julia> point_set!(body, :b, 3:4)

julia> precrack!(body, :a, :b)

julia> body
1000-point Body{BBMaterial{NoCorrection}}:
  3 point set(s):
    1000-point set `all_points`
    2-point set `a`
    2-point set `b`
  1 predefined crack(s)
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
    if point_sets_intersect(point_sets, key_a, key_b)
        msg = "point set `:$key_a` and `:$key_b` intersect!\n"
        msg *= "No point of the first set is allowed in the second!\n"
        throw(ArgumentError(msg))
    end
    return nothing
end
