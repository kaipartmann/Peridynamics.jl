const FIND_POINTS_ALLOWED_SYMBOLS = (:x, :y, :z, :p)
const SYMBOL_TO_DIM = Dict(:x => 0x1, :y => 0x2, :z => 0x3)

"""
    point_set!(body, set_name, points)
    point_set!(fun, body, set_name)

Add a point set to a [`Body`](@ref). The points of the set can be either specified directly
with the `points::AbstractVector` argument, or as the result of the filter function `fun`.
By default, a body already contains a point set with the name `:all_points`, containg a set
with all points.

# Arguments

- `body::AbstractBody`: [`Body`](@ref) where the set will be added.
- `set_name::Symbol`: Name of the point set.
- `points::AbstractVector`: Some vector containing the point indices of the set.
    The indices have to be in bounds with the `position` and `volume` of `body`.
- `fun::Function`: Function for filtering points. This function accepts only one positional
    argument and will be used in a `findall` call. Depending on the argument name,
    a different input will be processed:
    - `x`: The function will receive the x-coordinate of each point in `position` of `body`:
      ```julia
      points = findall(fun, @view(position[1, :]))
      ```
    - `y`: The function will receive the y-coordinate of each point in `position` of `body`:
      ```julia
      points = findall(fun, @view(position[2, :]))
      ```
    - `z`: The function will receive the z-coordinate of each point in `position` of `body`:
      ```julia
      points = findall(fun, @view(position[3, :]))
      ```
    - `p`: The function will receive the a vector containing each dimension of each point in
      `position` of `body`:
      ```julia
      points = findall(fun, eachcol(position))
      ```

# Throws

- Error if a point set with the same `set_name` already exists.
- Error if `points` are not in bounds with `position` and `volume` of the `body`.

# Examples

Add a point set to `body` with all points that have a x-corrdinate larger than zero:
```julia-repl
julia> point_set!(x -> x > 0, body, :larger_than_zero)

julia> point_sets(body)
Dict{Symbol, Vector{Int64}} with 2 entries:
  :larger_than_zero => [6, 7, 8, 9, 10, 16, 17, 18, 19, 20  …  9…
  :all_points       => [1, 2, 3, 4, 5, 6, 7, 8, 9, 10  …  991, 9…
```

Add a point set to `body` with all points that are positioned inside a sphere with
radius `r` around the center. Note that the `do`-syntax can be used, as `fun` is the first
argument of `point_set!`:
```julia-repl
julia> point_set!(body, :inside_sphere) do p
           sqrt(p[1]^2 + p[2]^2 + p[3]^2) ≤ r
       end

julia> point_sets(body)
Dict{Symbol, Vector{Int64}} with 2 entries:
  :larger_than_zero => [6, 7, 8, 9, 10, 16, 17, 18, 19, 20  …  9…
  :inside_sphere    => [1, 2, 3, 4, 5, 6, 7, 8, 9, 10  …  991, 9…
```
"""
function point_set! end

function point_set!(body::AbstractBody, set_name::Symbol,
                    points::V) where {V<:AbstractVector}
    checkbounds(body.volume, points)
    _point_set!(body.point_sets, set_name, points)
    return nothing
end

function point_set!(f::F, body::AbstractBody, set_name::Symbol) where {F<:Function}
    points = find_points(f, body.position)
    point_set!(body, set_name, points)
    return nothing
end

function _point_set!(point_sets::Dict{Symbol,Vector{Int}}, set_name::Symbol,
                     points::V) where {V<:AbstractVector{<:Integer}}
    if haskey(point_sets, set_name)
        error("there is already a point set with name $(set_name)!")
    end
    point_sets[set_name] = points
    return nothing
end

"""
    find_points(f, position)

$(internal_api_warning())

Find all points whose `position`s meet function `f`.

The function `f` accepts only one positional argument and will be used in a `findall` call.
Depending on the argument name, a different input will be processed.
See [`point_set!`](@ref).
"""
function find_points(f::F, position::AbstractMatrix) where {F<:Function}
    func_method = get_method_of_function(f)
    argnames = get_argument_names_of_function(func_method)
    if isempty(argnames)
        msg = "no argument in specified function $(nameof(f))!\n\n"
        msg *= "The filter function should contain one argument with name x, y, z or p!\n"
        throw(ArgumentError(msg))
    end
    if length(argnames) > 1
        msg = "too many arguments in specified function $(nameof(f))!\n\n"
        msg *= "The filter function should contain one argument with name x, y, z or p!\n"
        throw(ArgumentError(msg))
    end
    arg = first(argnames)
    if !in(arg, FIND_POINTS_ALLOWED_SYMBOLS)
        msg = "unknown argument name in specified function $(nameof(f))!\n\n"
        msg *= "The filter function should contain one argument with name x, y, z or p!\n"
        throw(ArgumentError(msg))
    end
    arg === :p && return findall(f, eachcol(position))
    dim = SYMBOL_TO_DIM[arg]
    return findall(f, @view(position[dim, :]))
end

"""
    check_if_set_is_defined(point_sets, set_name)

$(internal_api_warning())

Throw an error if no point set `set_name` is found in the dictionary `point_sets`.
"""
function check_if_set_is_defined(point_sets::Dict{Symbol,V}, set_name::Symbol) where {V}
    if !haskey(point_sets, set_name)
        error("there is no point set with name $(set_name)!")
    end
    return nothing
end

"""
    point_sets(body)

Return all point sets of `body`.

# Arguments

- `body::AbstractBody`: [`Body`](@ref).

# Example

```julia-repl
julia> body = Body(BBMaterial(), rand(3,100), rand(100))
100-point Body{BBMaterial{NoCorrection}}:
  100-point set `all_points`

julia> point_set!(body, :set_a, 1:10) # first ten points

julia> Peridynamics.point_sets(body)
Dict{Symbol, Vector{Int64}} with 2 entries:
  :all_points => [1, 2, 3, 4, 5, 6, 7, 8, 9, 10  …  91, 92, 93, 94, 95, 96, 97, 98, 9…
  :set_a      => [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
```
"""
function point_sets(body::AbstractBody)
    return body.point_sets
end

@inline function clean_point_set_name(name_str::AbstractString)
    name_str_cleaned = replace(name_str, " " => "_")
    return name_str_cleaned
end

"""
    point_sets_intersect(point_sets, key_a, key_b)

$(internal_api_warning())

Return `true` if point sets `key_a` and `key_b` out of the dictionary `point_sets` have
at least one common point, return false if they do not.
"""
function point_sets_intersect(point_sets::Dict{Symbol,Vector{Int}}, key_a::Symbol,
                              key_b::Symbol)
    isempty(point_sets[key_a] ∩ point_sets[key_b]) && return false
    return true
end
