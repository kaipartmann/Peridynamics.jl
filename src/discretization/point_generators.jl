"""
    uniform_box(lx, ly, lz, ΔX0; kwargs...)

Create a grid of uniformly distributed points in a cuboid with lengths `lx`, `ly` and `lz`
and point spacing `ΔX0`.

# Arguments
- `lx::Real`: Length in x-dimension.
- `ly::Real`: Length in y-dimension.
- `lz::Real`: Length in z-dimension.
- `ΔX0::Real`: Spacing of the points.

# Keywords
- `center`: The coordinates of the center of the cuboid. Default: `(0, 0, 0)`

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
function uniform_box(lx::Real, ly::Real, lz::Real, ΔX0::Real; center=(0, 0, 0))
    center_x, center_y, center_z = center
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

Create a grid of uniformly distributed points in a sphere with a specific `diameter` and
the point spacing `ΔX0`. Due to the uniform point spacings, edges on the surface of the
sphere can occur.

# Arguments
- `diameter::Real`: Diameter of the sphere.
- `ΔX0::Real`: Spacing of the points.

# Keywords
- `center`: The coordinates of the center of the sphere. Default: `(0, 0, 0)`

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
function uniform_sphere(diameter::Real, ΔX0::Real; center=(0, 0, 0))
    center_x, center_y, center_z = center
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

"""
    uniform_cylinder(diameter::Real, height::Real, ΔX0::Real; kwargs...)

Create a grid of uniformly distributed points in a cylindrical shape with a specific
`height` in z-dimension, a `diameter` and the point spacing `ΔX0`. Due to the uniform point
spacings, edges on the surface of the cylinder can occur.

# Arguments
- `diameter::Real`: Diameter of the cylinder.
- `height::Real`: Height of the cylinder.
- `ΔX0::Real`: Spacing of the points.

# Keywords
- `center`: The coordinates of the center of the cylinder. Default: `(0, 0, 0)`

# Returns
- `position::Matrix{Float64}`: A `3×n_points` matrix with the position of the points.
- `volume::Vector{Float64}`: A vector with the volume of each point.

# Examples

```julia-repl
julia> position, volume = uniform_cylinder(5, 10, 2);

julia> position
3×20 Matrix{Float64}:
 -1.0   1.0  -1.0   1.0  -1.0   1.0  -1.0  …   1.0  -1.0  1.0  -1.0   1.0  -1.0  1.0
 -1.0  -1.0   1.0   1.0  -1.0  -1.0   1.0     -1.0   1.0  1.0  -1.0  -1.0   1.0  1.0
 -4.0  -4.0  -4.0  -4.0  -2.0  -2.0  -2.0      2.0   2.0  2.0   4.0   4.0   4.0  4.0

julia> volume
20-element Vector{Int64}:
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
function uniform_cylinder(diameter::Real, height::Real, ΔX0::Real; center=(0, 0, 0))
    radius = diameter / 2
    center_x, center_y, center_z = center
    _gridxy = range(-radius + ΔX0 / 2, radius - ΔX0 / 2; step=ΔX0)
    gridxy = _gridxy .- sum(_gridxy) / length(_gridxy)
    _gridz = range((-height + ΔX0) / 2, (height - ΔX0) / 2; step=ΔX0)
    gridz = _gridz .- sum(_gridz) / length(_gridz)
    __position = vec(collect(Iterators.product(gridxy, gridxy, gridz)))
    _position = reinterpret(reshape, eltype(eltype(__position)), __position)
    cylinder_points = find_points(p -> √(p[1]^2 + p[2]^2) ≤ radius, _position)
    position = _position[:, cylinder_points]
    volume = fill(ΔX0^3, size(position, 2))
    isapprox(center_x, 0; atol=eps()) || (position[1, :] .+= center_x)
    isapprox(center_y, 0; atol=eps()) || (position[2, :] .+= center_y)
    isapprox(center_z, 0; atol=eps()) || (position[3, :] .+= center_z)
    return position, volume
end

# ###################################################################


"""
    round_sphere(diameter, ΔX0; kwargs...)

Create a grid of points distributed in a smooth sphere without edges on the surface
with a specific `diameter` and the point spacing `ΔX0`. Internally, some parts of
`TrixiParticles.jl` were copied and adapted for this function.

# Arguments
- `diameter::Real`: Diameter of the sphere.
- `ΔX0::Real`: Spacing of the points.

# Keywords
- `center`: The coordinates of the center of the sphere. Default: `(0, 0, 0)`

# Returns
- `position::Matrix{Float64}`: A `3×n_points` matrix with the position of the points.
- `volume::Vector{Float64}`: A vector with the volume of each point.

# Examples

```julia-repl
julia> position, volume = round_sphere(10, 2);

julia> position
3×63 Matrix{Float64}:
 0.0  1.74586  0.539501  -1.41243  -1.41243  …  -1.95518   -0.48289   1.35303   0.0
 0.0  0.0      1.66041    1.02619  -1.02619     -0.941566  -2.11568  -1.69664   0.0
 2.0  0.97569  0.97569    0.97569   0.97569     -3.36017   -3.36017  -3.36017  -4.0

julia> volume
63-element Vector{Float64}:
 8.311091676163475
 8.311091676163475
 8.311091676163475
 8.311091676163475
 8.311091676163475
 ⋮
 8.311091676163475
 8.311091676163475
 8.311091676163475
 8.311091676163475
```
"""
function round_sphere(diameter::Real, ΔX0::Real; center=(0, 0, 0))
    radius = diameter / 2
    NDIMS = length(center)
    coordinates = sphere_shape_coords(ΔX0, radius, SVector{NDIMS}(center))
    n_points = size(coordinates, 2)
    volumes = zeros(n_points)
    volumes .= 4/3 * π * radius^3 / n_points

    return coordinates, volumes
end

"""
    round_cylinder(diameter::Real, height::Real, ΔX0::Real; kwargs...)

Create a grid of  points distributed in a cylindrical shape with a specific
`height` in z-dimension, a `diameter` and the point spacing `ΔX0`. Due to a concentric
point distribution in the x-y-plane, there are no sharp edges that appear on the surface
of the cylinder.

# Arguments
- `diameter::Real`: Diameter of the cylinder.
- `height::Real`: Height of the cylinder.
- `ΔX0::Real`: Spacing of the points.

# Keywords
- `center`: The coordinates of the center of the cylinder. Default: `(0, 0, 0)`

# Returns
- `position::Matrix{Float64}`: A `3×n_points` matrix with the position of the points.
- `volume::Vector{Float64}`: A vector with the volume of each point.

# Examples

```julia-repl
julia> position, volume = round_cylinder(5, 10, 2);

julia> position
3×30 Matrix{Float64}:
  1.5   0.463525  -1.21353   -1.21353   …  -1.21353   -1.21353    0.463525
  0.0   1.42658    0.881678  -0.881678      0.881678  -0.881678  -1.42658
 -5.0  -5.0       -5.0       -5.0           5.0        5.0        5.0

julia> volume
30-element Vector{Int64}:
 13.089969389957473
 13.089969389957473
 13.089969389957473
 13.089969389957473
 13.089969389957473
  ⋮
 13.089969389957473
 13.089969389957473
 13.089969389957473
 13.089969389957473
```
"""
function round_cylinder(diameter::Real, height::Real, ΔX0::Real; center=(0, 0, 0))
    radius = diameter / 2

    xy = sphere_shape_coords(ΔX0, radius, SVector{2}((center[1], center[2])))
    _z = range(-height/2, height/2, step=ΔX0)
    z = _z .- sum(_z) / length(_z) .+ center[3]
    n_layers = length(z)
    n_points_per_layer = size(xy, 2)
    n_points = n_points_per_layer * n_layers
    coordinates = zeros(3, n_points)
    a, b = 1, 0
    for layer_id in 1:n_layers
        b += n_points_per_layer
        coordinates[1:2, a:b] .= xy
        coordinates[3, a:b] .= z[layer_id]
        a += n_points_per_layer
    end
    volumes = zeros(n_points)
    volumes .= 2 * π * radius^2 * height / n_points
    return coordinates, volumes
end

##########################################

function sphere_shape_coords(particle_spacing, radius, center)
    # Each layer has thickness `particle_spacing`
    n_layers = round(Int, radius / particle_spacing)

    if n_layers < 1
        # Just return one particle at the center
        return collect(reshape(center, (length(center), 1)))
    end

    # Same as above, which puts the inner radius between 0 and `particle_spacing`
    inner_radius = max(0.0, radius - n_layers * particle_spacing + 0.5particle_spacing)

    coords = zeros(length(center), 0)

    for layer in 0:(n_layers - 1)
        sphere_coords = _round_sphere(particle_spacing,
                                     inner_radius + layer * particle_spacing, center)
        coords = hcat(coords, sphere_coords)
    end

    return coords
end

function _round_sphere(particle_spacing, radius, center::SVector{2})

    n_particles = round(Int, 2pi * radius / particle_spacing)

    if n_particles <= 2
        # 2 or less particles produce weird, asymmetric results.
        # Just return one particle at the center.
        return collect(reshape(center, (2, 1)))
    end

    # Remove the last particle at 2pi, which overlaps with the first at 0
    t = LinRange(0, 2pi, n_particles + 1)[1:(end - 1)]

    particle_coords = Array{Float64, 2}(undef, 2, length(t))

    for i in axes(particle_coords, 2)
        particle_coords[:, i] = center + radius * SVector(cos(t[i]), sin(t[i]))
    end

    return particle_coords
end

function _round_sphere(particle_spacing, radius, center::SVector{3})
    # The number of particles can either be calculated in 2D or in 3D.
    # Let δ be the particle spacing and r the sphere radius.
    #
    # The volume of a particle is δ^3 and the volume of the sphere shell with
    # inner radius r - δ/2 and outer radius r + δ/2 is 4pi/3 * ((r + δ/2)^3 - (r - δ/2)^3).
    # The number of particles is then
    # n = 4pi / (3 δ^3) * ((r + δ/2)^3 - (r - δ/2)^3) = 4pi r^2 / δ^2 + pi/3.
    #
    # For small numbers of particles, we get better results without the term pi/3.
    # Omitting the term for the inner layers yields results with only ~5 particles less than
    # the theoretically optimal number of particles for the target density.
    n_particles = round(Int, 4pi * radius^2 / particle_spacing^2 + pi / 3)
    if n_particles < 300
        n_particles = round(Int, 4pi * radius^2 / particle_spacing^2)
    end

    # With fewer than 5 particles, this doesn't work properly
    if n_particles < 5
        if n_particles == 4
            # Return tetrahedron
            return [+1 -1 -1 +1;
                    +1 -1 +1 -1;
                    +1 +1 -1 -1] * radius / sqrt(3) .+ center
        elseif n_particles == 3
            # Return 2D triangle
            y = sin(2pi / 3)
            return [1 -0.5 -0.5;
                    0 y -y;
                    0 0 0] * radius .+ center
        elseif n_particles == 2
            # Return two particles
            return [-1 1;
                    0 0;
                    0 0] * radius .+ center
        else
            return collect(reshape(center, (3, 1)))
        end
    end

    # The following is a slightly adapted version of the "recursive zonal equal area
    # partition" of the sphere as explained by Leopardi (2006).
    #
    # With the equal area partition, the density at the poles is too high.
    # Instead, we slightly increase the area of the poles and modify the algorithm
    # accordingly.
    #
    # References:
    # - Paul Leopardi.
    #   "A partition of the unit sphere into regions of equal area and small diameter".
    #   In: Electronic Transactions on Numerical Analysis 25 (2006), pages 309-327.
    #   [http://eudml.org/doc/129860](http://eudml.org/doc/129860).

    # This is the Θ function, which is defined by Leopardi only as the inverse of V, without
    # giving a closed formula.
    theta(v) = acos(1 - v / 2pi)

    # Ideal area of the equal area partition
    ideal_area = 4pi / n_particles

    # Increase polar area to avoid higher density at the poles
    polar_area = 1.23ideal_area

    polar_radius = theta(polar_area)

    # Divide the remaining surface area equally
    collar_cell_area = (4pi - 2polar_area) / (n_particles - 2)

    # Strictly following Leopardi here. The collars should have equiangular spacing.
    collar_angle = sqrt(collar_cell_area)
    n_collars = max(1, round(Int, (pi - 2polar_radius) / collar_angle))
    fitting_collar_angle = (pi - 2polar_radius) / n_collars

    collar_area = [2pi * (cos(polar_radius + (j - 2) * fitting_collar_angle) -
                    cos(polar_radius + (j - 1) * fitting_collar_angle))
                   for j in 2:(n_collars + 1)]

    # Here, we count the poles as well
    ideal_number_cells = collar_area / collar_cell_area
    pushfirst!(ideal_number_cells, 1)
    push!(ideal_number_cells, 1)

    # Cumulative rounding to maintain the total number of cells
    actual_number_cells = ones(Int, length(ideal_number_cells))
    a = zeros(length(ideal_number_cells))
    for j in 2:(n_collars + 1)
        actual_number_cells[j] = round(Int, ideal_number_cells[j] + a[j - 1])

        a[j] = a[j - 1] + ideal_number_cells[j] - actual_number_cells[j]
    end

    collar_start_latitude = [theta(polar_area +
                                   sum(actual_number_cells[2:(j - 1)]) * collar_cell_area)
                             for j in 2:(n_collars + 2)]

    # Put particles in the center of each collar
    collar_latitude = [0.5 * (collar_start_latitude[i] + collar_start_latitude[i + 1])
                       for i in 1:n_collars]

    # Put the first and last particle on the pole
    pushfirst!(collar_latitude, 0.0)
    push!(collar_latitude, pi)

    # To compute the particle positions in each collar, we use the 2D `round_sphere`
    # function to generate a circle.
    particle_coords = zeros(3, 0)

    for circle in 1:(n_collars + 2)
        z = radius * cos(collar_latitude[circle])
        circle_radius = radius * sin(collar_latitude[circle])

        circle_spacing = 2pi * circle_radius / actual_number_cells[circle]

        # At the poles, `circle_radius` is zero, so we can pass any positive spacing
        if circle_spacing < eps()
            circle_spacing = 1.0
        end

        circle_coords_2d = _round_sphere(circle_spacing, circle_radius,
                                        SVector(center[1], center[2]))
        circle_coords_3d = vcat(circle_coords_2d,
                                center[3] .+ z * ones(1, size(circle_coords_2d, 2)))

        particle_coords = hcat(particle_coords, circle_coords_3d)
    end

    return particle_coords
end

#################################################################

"""
    rotate!(position, dimension, angle)

Rotate a set of 3D points in-place about one of the coordinate axes.

# Arguments
- `position::AbstractMatrix{<:Real}`: A `3×n` matrix of point coordinates (each
    column is a 3-D point). The matrix is modified in-place.
- `dimension::Union{Integer,Symbol}`: The axis to rotate about. Either specified as Symbol or
    integer.
    -  x-direction: `:x` or `1`
    -  y-direction: `:y` or `2`
    -  z-direction: `:z` or `3`
- `angle::Real`: Rotation angle in degrees. A positive angle performs a
    right-handed rotation about the specified axis.

# Examples

```julia-repl
julia> # rotate a point cloud 30° about the z-axis
julia> rotate!(positions, 3, 30.0)
```
"""
function rotate!(position, dimension::Union{Integer,Symbol}, angle)
    # Create rotation matrix
    dim = get_dim(dimension)
    if dim == 1
        R = SMatrix{3,3,Float64,9}(
                    1, 0, 0, 0, cosd(angle), sind(angle), 0, -sind(angle), cosd(angle))
    elseif dim == 2
        R = SMatrix{3,3,Float64,9}(
        cosd(angle), 0, -sind(angle), 0, 1, 0, sind(angle), 0, cosd(angle))
    else
        R = SMatrix{3,3,Float64,9}(
        cosd(angle), sind(angle), 0, -sind(angle), cosd(angle), 0, 0, 0, 1)
    end

    # Apply rotation to all position vectors
    for i in axes(position, 2)
        pos = get_vector(position, i)
        update_vector!(position, i, R * pos)
    end
    return nothing
end

function rotate!(position; angles=(0,0,0))
    for i in 1:3
        rotate!(position, i, angles[i])
    end
    return nothing
end
