# """
#     PointSetsPreCrack(set_a::Vector{Int}, set_b::Vector{Int})

# Definition of an preexisting crack in the model. Points in `set_a` cannot have
# interactions with points in `set_b`.

# # Fields
# - `set_a::Vector{Int}`: first point-id set
# - `set_b::Vector{Int}`: second point-id set
# """
# struct PointSetsPreCrack{U,V<:AbstractVector{<:Integer}} <: AbstractPredefinedCrack
#     set_a::U
#     set_b::V

#     function PointSetsPreCrack(set_a::U, set_b::V) where {U,V<:AbstractVector{<:Integer}}
#         if !isempty(set_a ∩ set_b)
#             msg = "set_a and set_b contain some same indices!\n"
#             msg *= "There should be no intersection between the two sets!\n"
#             throw(ArgumentError(msg))
#         end
#         return new{U,V}(set_a, set_b)
#     end
# end

struct PointSetsPreCrack <: AbstractPredefinedCrack
    set_a::Symbol
    set_b::Symbol
end
