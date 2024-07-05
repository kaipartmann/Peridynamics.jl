"""
    uniform_box(lx, ly, lz, ΔX0; kwargs...)

Creates a grid of uniformly distributed points in a cuboid with lengths `lx`, `ly` and `lz`
and point spacing `ΔX0`.

# Arguments
- `lx::Real`: Length in x-dimension.
- `ly::Real`: Length in y-dimension.
- `lz::Real`: Length in z-dimension.
- `ΔX0::Real`: The point spacing of the points.

# Keywords
- `center_x::Real`: Center of the cuboid in x-direction. (default: `0`)
- `center_y::Real`: Center of the cuboid in y-direction. (default: `0`)
- `center_z::Real`: Center of the cuboid in z-direction. (default: `0`)

# Returns
- `position::Matrix{Float64}`: A `3×n_points` matrix with the position of the points.
- `volume::Vector{Float64}`: A vector with the volume of each point.

# Examples

```julia-repl
julia> position, volume = uniform_box(10, 10, 10, 2);

julia> position
3×125 Matrix{Float64}:
 -4.0  -2.0   0.0   2.0   4.0  -4.0  -2.0  …  0.0  2.0  4.0  -4.0  -2.0  0.0  2.0  4.0
 -4.0  -4.0  -4.0  -4.0  -4.0  -2.0  -2.0     2.0  2.0  2.0   4.0   4.0  4.0  4.0  4.0
 -4.0  -4.0  -4.0  -4.0  -4.0  -4.0  -4.0     4.0  4.0  4.0   4.0   4.0  4.0  4.0  4.0

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
function uniform_box(lx::Real, ly::Real, lz::Real, ΔX0::Real;
                     center_x::Real=0, center_y::Real=0, center_z::Real=0)
    _gridx = range((-lx + ΔX0) / 2, (lx - ΔX0) / 2; step=ΔX0)
    gridx = _gridx .- sum(_gridx) / length(_gridx)
    _gridy = range((-ly + ΔX0) / 2, (ly - ΔX0) / 2; step=ΔX0)
    gridy = _gridy .- sum(_gridy) / length(_gridy)
    _gridz = range((-lz + ΔX0) / 2, (lz - ΔX0) / 2; step=ΔX0)
    gridz = _gridz .- sum(_gridz) / length(_gridz)
    _position = vec(collect(Iterators.product(gridx, gridy, gridz)))
    position = copy(reinterpret(reshape, eltype(eltype(_position)), _position))
    volume = fill(ΔX0^3, size(position, 2))
    isapprox(center_x, 0; atol=eps()) || (position[1, :] .+= center_x)
    isapprox(center_y, 0; atol=eps()) || (position[2, :] .+= center_y)
    isapprox(center_z, 0; atol=eps()) || (position[3, :] .+= center_z)
    return position, volume
end

"""
    uniform_sphere(diameter, ΔX0; kwargs...)

Creates a grid of uniformly distributed points in a sphere with a specific `diameter` and
the point spacing `ΔX0`.

# Arguments
- `diameter::Real`: Diameter of the sphere.
- `ΔX0::Real`: The point spacing of the points.

# Keywords
- `center_x::Real`: Center of the cuboid in x-direction. (default: `0`)
- `center_y::Real`: Center of the cuboid in y-direction. (default: `0`)
- `center_z::Real`: Center of the cuboid in z-direction. (default: `0`)

# Returns
- `position::Matrix{Float64}`: A `3×n_points` matrix with the position of the points.
- `volume::Vector{Float64}`: A vector with the volume of each point.

# Examples

```julia-repl
julia> position, volume = uniform_sphere(10, 2);

julia> position
3×81 Matrix{Float64}:
 -2.0   0.0   2.0  -2.0   0.0   2.0  -2.0  …   0.0   2.0  -2.0  0.0  2.0  -2.0  0.0  2.0
 -2.0  -2.0  -2.0   0.0   0.0   0.0   2.0     -2.0  -2.0   0.0  0.0  0.0   2.0  2.0  2.0
 -4.0  -4.0  -4.0  -4.0  -4.0  -4.0  -4.0      4.0   4.0   4.0  4.0  4.0   4.0  4.0  4.0

julia> volume
81-element Vector{Int64}:
 8
 8
 8
 8
 8
 ⋮
 8
 8
 8
 8
 8
```
"""
function uniform_sphere(diameter::Real, ΔX0::Real; center_x::Real=0, center_y::Real=0,
                        center_z::Real=0)
    radius = diameter / 2
    _grid = range(- radius + ΔX0 / 2, radius - ΔX0 / 2; step=ΔX0)
    grid = _grid .- sum(_grid) / length(_grid)
    __position = vec(collect(Iterators.product(grid, grid, grid)))
    _position = reinterpret(reshape, eltype(eltype(__position)), __position)
    sphere_points = find_points(p -> √(p[1]^2 + p[2]^2 + p[3]^2) ≤ radius, _position)
    position = _position[:, sphere_points]
    volume = fill(ΔX0^3, size(position, 2))
    isapprox(center_x, 0; atol=eps()) || (position[1, :] .+= center_x)
    isapprox(center_y, 0; atol=eps()) || (position[2, :] .+= center_y)
    isapprox(center_z, 0; atol=eps()) || (position[3, :] .+= center_z)
    return position, volume
end
