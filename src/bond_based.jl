"""
    BBMaterial <: AbstractMaterial

Bond based peridynamic material model.

# Fields
- `δ::Float64`: horizon
- `rho::Float64`: density
- `E::Float64`: young's modulus
- `nu::Float64`: poisson ratio
- `G::Float64`: shear modulus
- `K::Float64`: bulk modulus
- `bc::Float64`: bond constant
- `Gc::Float64`: critical energy release rate
- `εc::Float64`: critical bond strain

---
```julia
BBMaterial(; horizon::Real, rho::Real, E::Real, [Gc::Real, epsilon_c::Real])
```

Specify a material with `horizon`, density `rho`, Young's modulus `E` and either the
critical energy release rate `Gc` or the critical bond strain `epsilon_c`.

# Keywords
- `horizon::Real`: horizon
- `rho::Real`:density
- `E::Real`: young's modulus
- `Gc::Real`: critical energy release rate
- `epsilon_c::Real`: critical bond strain
"""
struct BBMaterial <: AbstractMaterial
    δ::Float64
    rho::Float64
    E::Float64
    nu::Float64
    G::Float64
    K::Float64
    bc::Float64
    Gc::Float64
    εc::Float64
end

function BBMaterial(; horizon::Real, rho::Real, E::Real, Gc::Real = -1,
                           epsilon_c::Real = -1)
    nu = 0.25 # limitation of the bond-based formulation
    G = E / (2 * (1 + nu))
    K = E / (3 * (1 - 2 * nu))
    bc = 18 * K / (π * horizon^4)
    if (Gc !== -1) && (epsilon_c == -1)
        epsilon_c = sqrt(5.0 * Gc / (9.0 * K * horizon))
    elseif (Gc == -1) && (epsilon_c !== -1)
        Gc = 9.0 / 5.0 * K * horizon * epsilon_c^2
    elseif (Gc !== -1) && (epsilon_c !== -1)
        msg = "Duplicate definition: define either Gc or epsilon_c, not both!"
        throw(ArgumentError(msg))
    elseif (Gc == -1) && (epsilon_c == -1)
        throw(ArgumentError("Either Gc or epsilon_c have to be defined!"))
    end
    return BBMaterial(horizon, rho, E, nu, G, K, bc, Gc, epsilon_c)
end

function Base.show(io::IO, ::MIME"text/plain", mat::BBMaterial)
    print(io, typeof(mat), ":")
    for field in fieldnames(typeof(mat))
        msg = "\n  " * rpad(string(field) * ":", 5) * string(getfield(mat, field))
        print(io, msg)
    end
    return nothing
end

# """
#     BondBasedBody

# Simulation body that contains all needed informations to run a simulation with bond-based
# peridynamics.

# # Fields
# - `n_points::Int`:
# - `n_bonds::Int`:
# - `n_threads::Int`:
# - `unique_bonds::Bool`:
# - `owned_points::Vector{UnitRange{Int}}`:
# - `owned_bonds::Vector{UnitRange{Int}}`:
# - `single_tids::Vector{Tuple{Int, Int}}`:
# - `multi_tids::Vector{Tuple{Int, Vector{Int}}}`:
# - `bond_data::Vector{Tuple{Int, Int, Float64, Bool}}`:
# - `volume::Vector{Float64}`:
# - `cells::Vector{MeshCell{VTKCellType, Tuple{Int64}}}`:
# - `n_family_members::Vector{Int}`:
# - `n_active_family_members::Matrix{Int}`:
# - `position::Matrix{Float64}`:
# - `displacement::Matrix{Float64}`:
# - `velocity::Matrix{Float64}`:
# - `velocity_half::Matrix{Float64}`:
# - `acceleration::Matrix{Float64}`:
# - `b_int::Array{Float64, 3}`:
# - `b_ext::Matrix{Float64}`:
# - `damage::Vector{Float64}`:
# - `bond_failure::Vector{Int}`:

# ---
# ```julia
# BondBasedBody(mat::PDMaterial{BBMaterial}, pc::PointCloud)
# ```

# # Arguments:
# - `mat::PDMaterial{BBMaterial}`:
# - `pc::PointCloud`:

# """
# struct BondBasedBody <: AbstractPDBody
#     n_points::Int
#     n_bonds::Int
#     n_threads::Int
#     unique_bonds::Bool
#     owned_points::Vector{UnitRange{Int}}
#     owned_bonds::Vector{UnitRange{Int}}
#     single_tids::Vector{Tuple{Int, Int}}
#     multi_tids::Vector{Tuple{Int, Vector{Int}}}
#     bond_data::Vector{Tuple{Int, Int, Float64, Bool}}
#     volume::Vector{Float64}
#     cells::Vector{MeshCell{VTKCellType, Tuple{Int64}}}
#     n_family_members::Vector{Int}
#     n_active_family_members::Matrix{Int}
#     position::Matrix{Float64}
#     displacement::Matrix{Float64}
#     velocity::Matrix{Float64}
#     velocity_half::Matrix{Float64}
#     acceleration::Matrix{Float64}
#     b_int::Array{Float64, 3}
#     b_ext::Matrix{Float64}
#     damage::Vector{Float64}
#     bond_failure::Vector{Int}
# end

# function BondBasedBody(mat::PDMaterial{BBMaterial}, pc::PointCloud)
#     n_threads = nthreads()
#     n_points = pc.n_points
#     @assert n_points>=n_threads "n_points < n_threads"
#     owned_points = defaultdist(n_points, n_threads)
#     volume = pc.volume
#     cells = get_cells(n_points)
#     position = pc.position
#     displacement = zeros(Float64, (3, n_points))
#     velocity = zeros(Float64, (3, n_points))
#     velocity_half = zeros(Float64, (3, n_points))
#     acceleration = zeros(Float64, (3, n_points))
#     forcedensity_ext = zeros(Float64, (3, n_points))
#     damage = zeros(Int, n_points)
#     forcedensity_int = zeros(Float64, (3, n_points, n_threads))
#     unique_bonds = true
#     bond_data, n_family_members = find_unique_bonds(pc, mat, owned_points)
#     n_bonds = length(bond_data)
#     owned_bonds = defaultdist(n_bonds, n_threads)
#     bond_failure = ones(Int, n_bonds)
#     n_active_family_members = zeros(Int, (n_points, n_threads))
#     _sum_tids = zeros(Bool, (n_points, n_threads))
#     _sum_tids .= false
#     @threads for tid in 1:n_threads
#         for current_bond in owned_bonds[tid]
#             (i, j, _, _) = bond_data[current_bond]
#             n_active_family_members[i, tid] += 1
#             n_active_family_members[j, tid] += 1
#             _sum_tids[i, tid] = true
#             _sum_tids[j, tid] = true
#         end
#     end
#     sum_tids = [findall(row) for row in eachrow(_sum_tids)]
#     single_tids, multi_tids = find_tids(sum_tids)
#     return BondBasedBody(n_points, n_bonds, n_threads, unique_bonds, owned_points,
#                          owned_bonds, single_tids, multi_tids, bond_data, volume, cells,
#                          n_family_members, n_active_family_members, position, displacement,
#                          velocity, velocity_half, acceleration, forcedensity_int,
#                          forcedensity_ext, damage, bond_failure)
# end

# function Base.show(io::IO, ::MIME"text/plain", body::BondBasedBody)
#     println(io, body.n_points, "-points ", typeof(body), ":")
#     print(io, "  ", body.n_bonds, " bonds")
#     return nothing
# end

# init_body(mat::PDMaterial{BBMaterial}, pc::PointCloud) = BondBasedBody(mat, pc)

struct BBSimParameters <: AbstractSimParameters
    n_chunks::Int
    mat::PDMaterial{BBMaterial}
    pc::PointCloud
    cells::Vector{MeshCell{VTKCellType, Tuple{Int64}}}
    n_bonds::Int
    bonds::Vector{Int}
    init_dists::Vector{Float64}
    failure_allowed::Vector{Bool}
    n_family_members::Vector{Int}
    hood_range::Vector{UnitRange{Int}}
    bcs::Vector{<:AbstractBC}
    ics::Vector{<:AbstractIC}

    function BBSimParameters(mat::PDMaterial{BBMaterial}, pc::PointCloud,
                             bcs::Vector{<:AbstractBC}, ics::Vector{<:AbstractIC},
                             n_chunks::Int)

        point_dist = defaultdist(pc.n_points, n_chunks)
        @assert length(point_dist) == n_chunks
        bonds, n_family_members, init_dists, failure_allowed = find_bonds(pc, mat.δ, point_dist)
        n_bonds = length(bonds)
        cells = get_cells(pc.n_points)
        hood_range = find_hood_range(n_family_members, pc.n_points)
        return new(n_chunks, mat, pc, cells, n_bonds, bonds, init_dists, failure_allowed,
                   n_family_members, hood_range, bcs, ics)
    end
end

struct BBGlobalStorage <: AbstractGlobalStorage
    position::Matrix{Float64}
    displacement::Matrix{Float64}
    velocity::Matrix{Float64}
    velocity_half::Matrix{Float64}
    acceleration::Matrix{Float64}
    bond_failure::Vector{Int}
    damage::Vector{Float64}
    b_int::Matrix{Float64}
    b_ext::Matrix{Float64}
    n_active_family_members::Vector{Int}

    function BBGlobalStorage(sp::BBSimParameters)
        position = copy(sp.pc.position)
        displacement = zeros(3, sp.pc.n_points)
        velocity = zeros(3, sp.pc.n_points)
        velocity_half = zeros(3, sp.pc.n_points)
        acceleration = zeros(3, sp.pc.n_points)
        b_int = zeros(3, sp.pc.n_points)
        b_ext = zeros(3, sp.pc.n_points)
        bond_failure = ones(Int, sp.n_bonds)
        damage = zeros(sp.pc.n_points)
        n_active_family_members = copy(sp.n_family_members)
        return new(position, displacement, velocity, velocity_half, acceleration,
                   bond_failure, damage, b_int, b_ext, n_active_family_members)
    end
end

struct BBChunkLocalStorage <: AbstractChunkLocalStorage
    cl_points::UnitRange{Int}
    pointmap::Dict{Int,Int}
    b_int::Matrix{Float64}
    n_active_family_members::Vector{Int}
end

function init_chunk_local_storages(sp::BBSimParameters)
    vcls = Vector{BBChunkLocalStorage}(undef, sp.n_chunks)
    @threads :static for cid in 1:n_chunks
        cl_points = sp.point_dist[cid]
        pointmap = Dict{Int,Int}()
        n_local_points = 0
        for i in cl_points
            if i in keys(pointmap)
                error("key $i is already defined in the pointmap dict!\n")
            end
            n_local_points += 1
            pointmap[i] = n_local_points
        end
        @assert length(pointmap) == n_local_points
        b_int_local = zeros(3, n_local_points)
        n_active_family_members_local = zeros(n_local_points)
        for (gi, li) in pointmap
            if li > n_local_points
                error("local_index > n_local_points!\n")
            end
            n_active_family_members_local[li] = n_family_members[gi]
        end
        vcls[cid] = ChunkLocalStorage(cl_points, pointmap, b_int_local,
                                    n_active_family_members_local)
    end
    return vcls
end

function compute_forcedensity!(body::BondBasedBody, mat::PDMaterial{BBMaterial})
    body.b_int .= 0.0
    body.n_active_family_members .= 0
    @inbounds @threads for tid in 1:body.n_threads
        for current_bond in body.owned_bonds[tid]
            (i, j, ξ, failureflag) = body.bond_data[current_bond]
            ϑx = body.position[1, j] - body.position[1, i]
            ϑy = body.position[2, j] - body.position[2, i]
            ϑz = body.position[3, j] - body.position[3, i]
            η = sqrt(ϑx * ϑx + ϑy * ϑy + ϑz * ϑz)
            ε = (η - ξ) / ξ
            temp = body.bond_failure[current_bond] * ε / η
            tempi = temp * mat[i].bc * body.volume[j]
            tempj = temp * mat[j].bc * body.volume[i]
            body.b_int[1, i, tid] += tempi * ϑx
            body.b_int[2, i, tid] += tempi * ϑy
            body.b_int[3, i, tid] += tempi * ϑz
            body.b_int[1, j, tid] -= tempj * ϑx
            body.b_int[2, j, tid] -= tempj * ϑy
            body.b_int[3, j, tid] -= tempj * ϑz
            if ε > mat[i].εc || ε > mat[j].εc
                if failureflag
                    body.bond_failure[current_bond] = 0
                end
            end
            body.n_active_family_members[i, tid] += body.bond_failure[current_bond]
            body.n_active_family_members[j, tid] += body.bond_failure[current_bond]
        end
    end
    return nothing
end
