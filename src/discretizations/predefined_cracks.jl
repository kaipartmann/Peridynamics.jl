"""
    PointSetsPreCrack(set_a::Vector{Int}, set_b::Vector{Int})

Definition of an preexisting crack in the model. Points in `set_a` cannot have
interactions with points in `set_b`.

# Fields
- `set_a::Vector{Int}`: first point-id set
- `set_b::Vector{Int}`: second point-id set
"""
struct PointSetsPreCrack <: AbstractPredefinedCrack
    set_a::Vector{Int}
    set_b::Vector{Int}
end
