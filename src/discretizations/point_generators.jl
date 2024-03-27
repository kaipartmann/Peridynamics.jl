"""
    uniform_box(lx::Real, ly::Real, lz::Real, Δx::Real;
                center_x::Real=0, center_y::Real=0, center_z::Real=0)

Creates a grid of points in a cuboid form with lengths `lx`, `ly` and `lz` and the distance
`Δx` between the points.

# Arguments

- `lx::Real`: length in x-dimension
- `ly::Real`: length in y-dimension
- `lz::Real`: length in z-dimension
- `Δx::Real`: distance between neighboring points

# Keywords

- `center_x::Real=0`: center of the cuboid in x-dimension
- `center_y::Real=0`: center of the cuboid in y-dimension
- `center_z::Real=0`: center of the cuboid in z-dimension

# Returns

- `position::Matrix{<:Real}`: 3×n matrix with position of each point
- `volume::Vector{<:Real}`: vector with volume of each point

# Example

```julia-repl
julia> position, volume = uniform_box(10, 10, 10, 2);

julia> position
3×125 Matrix{Float64}:
 -4.0  -2.0   0.0   2.0   4.0  -4.0  -2.0   0.0  …  0.0  2.0  4.0  -4.0  -2.0  0.0  2.0  4.0
 -4.0  -4.0  -4.0  -4.0  -4.0  -2.0  -2.0  -2.0     2.0  2.0  2.0   4.0   4.0  4.0  4.0  4.0
 -4.0  -4.0  -4.0  -4.0  -4.0  -4.0  -4.0  -4.0     4.0  4.0  4.0   4.0   4.0  4.0  4.0  4.0

julia> volume
125-element Vector{Int64}:
 8
 8
 8
 8
 ⋮
 8
 8
 8
 8
```
"""
function uniform_box(lx::Real, ly::Real, lz::Real, Δx::Real;
                     center_x::Real=0, center_y::Real=0, center_z::Real=0)
    _gridx = range((-lx + Δx) / 2, (lx - Δx) / 2; step=Δx)
    gridx = _gridx .- sum(_gridx) / length(_gridx)
    _gridy = range((-ly + Δx) / 2, (ly - Δx) / 2; step=Δx)
    gridy = _gridy .- sum(_gridy) / length(_gridy)
    _gridz = range((-lz + Δx) / 2, (lz - Δx) / 2; step=Δx)
    gridz = _gridz .- sum(_gridz) / length(_gridz)
    _position = vec(collect(Iterators.product(gridx, gridy, gridz)))
    position = copy(reinterpret(reshape, eltype(eltype(_position)), _position))
    volume = fill(Δx^3, size(position, 2))
    center_x != 0 && (position[1, :] .+= center_x)
    center_y != 0 && (position[2, :] .+= center_y)
    center_z != 0 && (position[3, :] .+= center_z)
    return position, volume
end
