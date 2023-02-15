@doc raw"""
    PointCloud

Peridynamic spatial discretization with material points defining a point cloud.

# Fields
- `n_points::Int`: number of material points
- `position::Matrix{Float64}`: coordinates of points in reference configuration
- `volume::Vector{Float64}`: material point volumes
- `failure_flag::BitVector`: if failure of point is possible: element=`true`
- `radius::Vector{Float64}`: radius of the material point sphere

---
```julia
PointCloud(position, volume[, point_sets])
```

Create a `PointCloud` by specifying position and volume of all material points. The radius
of the point sphere is derived from the volume with the function [`sphere_radius`](@ref) and
the `failure_flag` set to `true` for all points.

# Arguments
- `position::Matrix{Float64}`: the position of the material points ($3 \times N$-matrix for
  $N$ material points)
- `volume::Vector{Float64}`: the volume of the material points
- `point_sets::Dict{String,Vector{Int}}=Dict{String,Vector{Int}}()`: optional point sets

# Returns
- `PointCloud`: point cloud with specified material point position and volume

# Examples

Creating a `PointCloud` with 4 manually defined points:
```julia-repl
julia> position = [
           0.0 1.0 0.0 0.0
           0.0 0.0 1.0 0.0
           0.0 0.0 0.0 1.0
       ]
3×4 Matrix{Float64}:
 0.0  1.0  0.0  0.0
 0.0  0.0  1.0  0.0
 0.0  0.0  0.0  1.0

julia> volume = [1.0, 1.0, 1.0, 1.0]
4-element Vector{Float64}:
 1.0
 1.0
 1.0
 1.0

julia> pc = PointCloud(position, volume)
4-points PointCloud

julia> pc.failure_flag
4-element BitVector:
 1
 1
 1
 1

julia> pc.radius
4-element Vector{Float64}:
 0.6203504908994001
 0.6203504908994001
 0.6203504908994001
 0.6203504908994001
```

---
```julia
PointCloud(lx, ly, lz, Δx; center_x=0, center_y=0, center_z=0)
```

Generate a uniformly distributed `PointCloud` with the lengths `lx`, `ly`, `lz`, and
point spacing `Δx`. Optional keyword arguments provide the possibility to set the center.

# Arguments
- `lx::Real`: length of the `PointCloud`-block in x-direction
- `ly::Real`: length of the `PointCloud`-block in y-direction
- `lz::Real`: length of the `PointCloud`-block in z-direction
- `Δx::Real`: point spacing in x-, y- and z- direction
- `center_x::Real=0`: x-coordinate of the `PointCloud`-block center
- `center_y::Real=0`: y-coordinate of the `PointCloud`-block center
- `center_z::Real=0`: z-coordinate of the `PointCloud`-block center

# Returns
- `PointCloud`: point cloud with with `W`, length `L`, height `H` and
  point spacing `Δx`

# Examples

Cube with side length 1 and point spacing $Δx = 0.1$:
```julia-repl
julia> PointCloud(1, 1, 1, 0.1)
1000-points PointCloud
```
"""
struct PointCloud
    n_points::Int
    position::Matrix{Float64}
    volume::Vector{Float64}
    failure_flag::BitVector
    radius::Vector{Float64}
    point_sets::Dict{String,Vector{Int}}
end

function PointCloud(
    position::Matrix{Float64},
    volume::Vector{Float64},
    point_sets::Dict{String,Vector{Int}}=Dict{String,Vector{Int}}(),
)
    if size(position, 1) !== 3 || size(position, 2) !== length(volume)
        throw(DimensionMismatch("size of position: $(size(position)) !== (3, n_points)"))
    end
    n_points = length(volume)
    radius = sphere_radius.(volume)
    failure_flag = BitVector(undef, n_points)
    failure_flag .= true
    return PointCloud(n_points, position, volume, failure_flag, radius, point_sets)
end

function PointCloud(
    lx::Real,
    ly::Real,
    lz::Real,
    Δx::Real;
    center_x::Real=0,
    center_y::Real=0,
    center_z::Real=0,
)
    _gridX = range(; start=(-lx + Δx) / 2, stop=(lx - Δx) / 2, step=Δx)
    gridX = _gridX .- sum(_gridX) / length(_gridX)
    _gridY = range(; start=(-ly + Δx) / 2, stop=(ly - Δx) / 2, step=Δx)
    gridY = _gridY .- sum(_gridY) / length(_gridY)
    _gridZ = range(; start=(-lz + Δx) / 2, stop=(lz - Δx) / 2, step=Δx)
    gridZ = _gridZ .- sum(_gridZ) / length(_gridZ)
    positions = hcat(([x; y; z] for x in gridX for y in gridY for z in gridZ)...)
    if center_x !== 0
        positions[1, :] .+= center_x
    end
    if center_y !== 0
        positions[2, :] .+= center_y
    end
    if center_z !== 0
        positions[3, :] .+= center_z
    end
    n_points = size(positions, 2)
    volumes = fill(Δx^3, n_points)
    radius = sphere_radius.(volumes)
    failure_flag = BitVector(undef, n_points)
    failure_flag .= true
    point_sets = Dict{String,Vector{Int}}()
    return PointCloud(n_points, positions, volumes, failure_flag, radius, point_sets)
end

@doc raw"""
    sphere_radius(vol::T) where {T<:Real}

Calculate the radius $r$ of the sphere by equation
```math
r = \sqrt[3]{\frac{3 \; V}{4 \; \pi}}
```
with specified sphere volume $V$.

# Arguments
- `vol::T where {T<:Real}`: volume $V$ of the sphere

# Returns
- `Float64`: radius $r$ of sphere with volume `vol`
"""
function sphere_radius(vol::T) where {T<:Real}
    r = (3 * vol / (4 * π))^(1 / 3)
    return r
end

function Base.show(io::IO, ::MIME"text/plain", pc::PointCloud)
    print(io, pc.n_points, "-points ", typeof(pc))
    return nothing
end

"""
    pcmerge(v::Vector{PointCloud})

Merge multiple point clouds into one [`PointCloud`](@ref).

# Arguments
- `v::Vector{PointCloud}`: vector of multiple `PointClouds`

# Returns
- `PointCloud`: merged point cloud
"""
function pcmerge(v::Vector{PointCloud})
    n_points = sum([pc.n_points for pc in v])
    position = reduce(hcat, [pc.position for pc in v])
    volume = reduce(vcat, [pc.volume for pc in v])
    failure_flag = reduce(vcat, [pc.failure_flag for pc in v])
    radius = reduce(vcat, [pc.radius for pc in v])
    point_sets = Dict{String,Vector{Int}}()
    n_points_cnt = first(v).n_points
    for pc_id in firstindex(v)+1:lastindex(v)
        for (key, val) in v[pc_id].point_sets
            point_sets[key] = val .+ n_points_cnt
        end
        n_points_cnt += v[pc_id].n_points
    end
    pc = PointCloud(n_points, position, volume, failure_flag, radius, point_sets)
    return pc
end

"""
    PreCrack(point_id_set_a::Vector{Int}, point_id_set_b::Vector{Int})

Definition of an preexisting crack in the model. Points in `point_id_set_a` cannot have
interactions with points in `point_id_set_b`.

# Fields
- `point_id_set_a::Vector{Int}`: first point-id set
- `point_id_set_b::Vector{Int}`: second point-id set
"""
struct PreCrack
    point_id_set_a::Vector{Int}
    point_id_set_b::Vector{Int}
end

function Base.show(io::IO, ::MIME"text/plain", precrack::PreCrack)
    println(io, typeof(precrack), ":")
    println(io, " ", length(precrack.point_id_set_a), " points in set a")
    print(io, " ", length(precrack.point_id_set_b), " points in set b")
    return nothing
end

function define_precrack!(body::AbstractPDBody, precrack::PreCrack)
    @inbounds @threads for tid in 1:body.n_threads
        (@view body.n_active_family_members[:, tid]) .= 0
        for current_one_ni in body.owned_bonds[tid]
            i, j, _, _ = body.bond_data[current_one_ni]
            i_is_in_set_a = in(i, precrack.point_id_set_a)
            i_is_in_set_b = in(i, precrack.point_id_set_b)
            j_is_in_set_a = in(j, precrack.point_id_set_a)
            j_is_in_set_b = in(j, precrack.point_id_set_b)
            if i_is_in_set_a && j_is_in_set_b || i_is_in_set_b && j_is_in_set_a
                body.bond_failure[current_one_ni] = 0
            end
            body.n_active_family_members[i, tid] += body.bond_failure[current_one_ni]
            if body.unique_bonds
                body.n_active_family_members[j, tid] += body.bond_failure[current_one_ni]
            end
        end
    end
    return nothing
end

function calc_damage!(body::AbstractPDBody)
    @inbounds @threads for i in 1:body.n_points
        body.damage[i] = 1 - body.n_active_family_members[i, 1] / body.n_family_members[i]
    end
    return nothing
end

function find_bonds(pc::PointCloud, δ::Float64, owned_points::Vector{UnitRange{Int}})
    n_threads = nthreads()
    _bond_data = fill([(0, 0, 0.0, true)], n_threads)
    n_family_members = zeros(Int, pc.n_points)
    p = Progress(pc.n_points;
        dt=1,
        desc="Search bonds...     ",
        barlen=30,
        color=:normal,
        enabled=!is_logging(stderr),
    )
    @threads for tid in 1:n_threads
        local_bond_data = Vector{Tuple{Int,Int,Float64,Bool}}(undef, 0)
        sizehint!(local_bond_data, pc.n_points * 500)
        idst = 0.0
        num = 0
        fail = false
        for a in owned_points[tid]
            num = 0
            for i in 1:pc.n_points
                if a !== i
                    idst = sqrt(
                        (pc.position[1, i] - pc.position[1, a])^2 +
                        (pc.position[2, i] - pc.position[2, a])^2 +
                        (pc.position[3, i] - pc.position[3, a])^2,
                    )
                    if idst <= δ
                        num += 1
                        fail = pc.failure_flag[a] & pc.failure_flag[i]
                        push!(local_bond_data, (a, i, idst, fail))
                    end
                end
            end
            n_family_members[a] = num
            next!(p)
        end
        _bond_data[tid] = local_bond_data
    end
    finish!(p)
    bond_data = reduce(append!, _bond_data)
    return bond_data, n_family_members
end

function find_unique_bonds(pc::PointCloud, δ::Float64, owned_points::Vector{UnitRange{Int}})
    n_threads = nthreads()
    _bond_data = fill([(0, 0, 0.0, true)], n_threads)
    n_family_members = zeros(Int, pc.n_points)
    p = Progress(pc.n_points;
        dt=1,
        desc="Search bonds...     ",
        barlen=30,
        color=:normal,
        enabled=!is_logging(stderr),
    )
    @threads for tid in 1:n_threads
        local_bonds_data = Vector{Tuple{Int,Int,Float64,Bool}}(undef, 0)
        sizehint!(local_bonds_data, pc.n_points * 500)
        idst = 0.0
        num = 0
        fail = false
        for a in owned_points[tid]
            num = 0
            for i in 1:pc.n_points
                if a !== i
                    idst = sqrt(
                        (pc.position[1, i] - pc.position[1, a])^2 +
                        (pc.position[2, i] - pc.position[2, a])^2 +
                        (pc.position[3, i] - pc.position[3, a])^2,
                    )
                    if idst <= δ
                        num += 1
                        fail = pc.failure_flag[a] & pc.failure_flag[i]
                        push!(
                            local_bonds_data,
                            a < i ? (a, i, idst, fail) : (i, a, idst, fail),
                        )
                    end
                end
            end
            n_family_members[a] = num
            next!(p)
        end
        _bond_data[tid] = local_bonds_data
    end
    finish!(p)
    one_ni_data = unique(reduce(append!, _bond_data))
    return one_ni_data, n_family_members
end
