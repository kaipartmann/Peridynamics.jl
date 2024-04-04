const FIND_POINTS_ALLOWED_SYMBOLS = (:x, :y, :z, :p)
const SYMBOL_TO_DIM = Dict(:x => 0x1, :y => 0x2, :z => 0x3)

"""
    point_set!(b::AbstractBody, name::Symbol, points::V) where {V<:AbstractVector}
    point_set!(f::F, b::AbstractBody, name::Symbol) where {F<:Function}

Creates the point set `name` in Body `b` containing all points either defined in `V` or
described by `F`

# Arguments

- `b::AbstractBody`: peridynamic body
- `name::Symbol`: name of the point set
- `points<:AbstractVector`: vector of point indices
- `f<:Function`: function that describes points contained in point set

# Throws

- error if a point set called `name` is already defined
- `BoundsError`: if points in `V` do not exist

# Example

```julia-repl
julia> point_set!(p -> p[1] ≤ -l/2+a && 0 ≤ p[2] ≤ 2δ, b, :set_a)
julia> point_set!(p -> p[1] ≤ -l/2+a && -2δ ≤ p[2] < 0, b, :set_b)
julia> b.point_sets
Dict{Symbol, Vector{Int64}} with 4 entries:
  :set_a      => [1251, 1252, 1253, 1254, 1255, 1256, 1257, 1258, 1259, …
  :set_b      => [951, 952, 953, 954, 955, 956, 957, 958, 959, 960  …  1…
```
"""
function point_set! end

function point_set!(b::AbstractBody, name::Symbol, points::V) where {V<:AbstractVector}
    checkbounds(b.volume, points)
    _point_set!(b.point_sets, name, points)
    return nothing
end

function point_set!(f::F, b::AbstractBody, name::Symbol) where {F<:Function}
    points = find_points(f, b.position)
    point_set!(b, name, points)
    return nothing
end

function _point_set!(point_sets::Dict{Symbol,Vector{Int}}, name::Symbol,
                     points::V) where {V<:AbstractVector{<:Integer}}
    if haskey(point_sets, name)
        error("there is already a point set with name $(name)!")
    end
    point_sets[name] = points
    return nothing
end

function find_points(f::F, position::Matrix{Float64}) where {F<:Function}
    func_method = get_method_of_function(f)
    argnames = get_argument_names_of_function(func_method)
    if isempty(argnames)
        msg =  "no argument in specified function $(nameof(f))!\n\n"
        msg *= "The filter function should contain one argument with name x, y, z or p!\n"
        throw(ArgumentError(msg))
    end
    if length(argnames) > 1
        msg =  "too many arguments in specified function $(nameof(f))!\n\n"
        msg *= "The filter function should contain one argument with name x, y, z or p!\n"
        throw(ArgumentError(msg))
    end
    arg = first(argnames)
    if !in(arg, FIND_POINTS_ALLOWED_SYMBOLS)
        msg =  "unknown argument name in specified function $(nameof(f))!\n\n"
        msg *= "The filter function should contain one argument with name x, y, z or p!\n"
        throw(ArgumentError(msg))
    end
    arg === :p && return findall(f, eachcol(position))
    dim = SYMBOL_TO_DIM[arg]
    return findall(f, @view(position[dim, :]))
end

function check_if_set_is_defined(point_sets::Dict{Symbol,V}, name::Symbol) where {V}
    if !haskey(point_sets, name)
        error("there is no point set with name $(name)!")
    end
    return nothing
end
