"""
    PointSetsPreCrack

$(internal_api_warning())

Type describing a predefined crack in a peridynamic body.

# Fields

- `set_a::Symbol`: Point set containing points on one side of the crack.
- `set_b::Symbol`: Point set with points on other side of the crack.
- `filter_bonds::Bool`: If true, the involved bonds are filtered out so no damage is
    present at the beginning of the simulation. Else, all involved bonds are marked broken
    from the beginning.
"""
struct PointSetsPreCrack <: AbstractPredefinedCrack
    set_a::Symbol
    set_b::Symbol
    filter_bonds::Bool
end

@inline filter_bonds(crack::AbstractPredefinedCrack) = crack.filter_bonds

"""
    apply_precracks!(chunk, body)

$(internal_api_warning())

Apply all predefined cracks for `chunk` by calling `apply_precrack!` for each crack.
"""
function apply_precracks!(chunk::AbstractBodyChunk, body::AbstractBody)
    for precrack in body.point_sets_precracks
        apply_precrack!(chunk, body, precrack)
    end
    calc_damage!(chunk)
    return nothing
end

"""
    apply_precrack!(chunk, body, crack)

$(internal_api_warning())

Apply the predefined `crack` of the `body` for the considered `chunk` by breaking the
concerned bonds.
"""
function apply_precrack!(chunk::AbstractBodyChunk, body::AbstractBody,
                         crack::PointSetsPreCrack)
    filter_bonds(crack) && return nothing
    point_ids = get_point_ids(chunk.system)
    set_a = filter(x -> in(x, point_ids), body.point_sets[crack.set_a])
    set_b = filter(x -> in(x, point_ids), body.point_sets[crack.set_b])
    localizer = get_localizer(chunk.system)
    localize!(set_a, localizer)
    localize!(set_b, localizer)
    if isempty(set_a) || isempty(set_b)
        return nothing
    end
    break_bonds!(chunk.storage, chunk.system, set_a, set_b)
    return nothing
end

"""
    precrack!(body, set_a, set_b; update_dmg=true)

Create a crack between two point sets by prohibiting interaction between points of
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

- Error if the body does not contain sets with name `set_a` and `set_b`.
- Error if the point sets intersect and a point is included in both sets.

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

"""
    check_if_sets_intersect(point_sets, key_a, key_b)

$(internal_api_warning())

Throw error if two sets chosen for `precrack!` have common points.
"""
function check_if_sets_intersect(point_sets::Dict{Symbol,Vector{Int}}, key_a::Symbol,
                                 key_b::Symbol)
    if point_sets_intersect(point_sets, key_a, key_b)
        msg = "point set `:$key_a` and `:$key_b` intersect!\n"
        msg *= "No point of the first set is allowed in the second!\n"
        throw(ArgumentError(msg))
    end
    return nothing
end
