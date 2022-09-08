"""
    BondBasedMaterial <: AbstractPDMaterial

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
BondBasedMaterial(; horizon::Real, rho::Real, E::Real, Gc::Real)
```

Specify a material only with `horizon`, density `rho`, Young's modulus `E` and critical
energy release rate `Gc`.

# Keywords
- `horizon::Real`: horizon
- `rho::Real`:density
- `E::Real`: young's modulus
- `Gc::Real`: critical energy release rate
"""
struct BondBasedMaterial <: AbstractPDMaterial
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

function BondBasedMaterial(; horizon::Real, rho::Real, E::Real, Gc::Real)
    nu = 0.25 # limitation of the bond-based formulation
    G = E / (2 * (1 + nu))
    K = E / (3 * (1 - 2 * nu))
    bc = 18 * K / (π * horizon^4)
    εc = sqrt(5.0 * Gc / (9.0 * K * horizon))
    return BondBasedMaterial(horizon, rho, E, nu, G, K, bc, Gc, εc)
end

function Base.show(io::IO, ::MIME"text/plain", mat::BondBasedMaterial)
    println(io, typeof(mat), ":")
    for field in fieldnames(typeof(mat))
        @printf(io, "  %s: %g\n", string(field), getfield(mat, field))
    end
    return nothing
end

struct BondBasedBody <: AbstractPDBody
    n_points::Int
    n_bonds::Int
    n_threads::Int
    unique_bonds::Bool
    owned_points::Vector{UnitRange{Int}}
    owned_bonds::Vector{UnitRange{Int}}
    bond_data::Vector{Tuple{Int,Int,Float64,Bool}}
    volumes::Vector{Float64}
    n_family_members::Vector{Int}
    n_active_family_members::Matrix{Int}
    position::Matrix{Float64}
    displacement::Matrix{Float64}
    velocity::Matrix{Float64}
    velocity_half::Matrix{Float64}
    acceleration::Matrix{Float64}
    b_int::Array{Float64,3}
    b_ext::Matrix{Float64}
    damage::Vector{Float64}
    bond_failure::Vector{Int}
end

function BondBasedBody(mat::BondBasedMaterial, pc::PointCloud)
    n_threads = nthreads()
    n_points = pc.n_points
    @assert n_points >= n_threads "n_points < n_threads"
    owned_points = defaultdist(n_points, n_threads)
    volumes = pc.volume
    position = pc.position
    displacement = zeros(Float64, (3, n_points))
    velocity = zeros(Float64, (3, n_points))
    velocity_half = zeros(Float64, (3, n_points))
    acceleration = zeros(Float64, (3, n_points))
    forcedensity_ext = zeros(Float64, (3, n_points))
    damage = zeros(Int, n_points)
    forcedensity_int = zeros(Float64, (3, n_points, n_threads))
    unique_bonds = true
    bond_data, n_family_members = find_unique_bonds(pc, mat.δ, owned_points)
    n_bonds = length(bond_data)
    owned_bonds = defaultdist(n_bonds, n_threads)
    bond_failure = ones(Int, n_bonds)
    n_active_family_members = zeros(Int, (n_points, n_threads))
    @threads :static for tid in 1:n_threads
        for i in owned_points[tid]
            n_active_family_members[i, tid] = n_family_members[i]
        end
    end
    return BondBasedBody(
        n_points,
        n_bonds,
        n_threads,
        unique_bonds,
        owned_points,
        owned_bonds,
        bond_data,
        volumes,
        n_family_members,
        n_active_family_members,
        position,
        displacement,
        velocity,
        velocity_half,
        acceleration,
        forcedensity_int,
        forcedensity_ext,
        damage,
        bond_failure,
    )
end

function Base.show(io::IO, ::MIME"text/plain", body::BondBasedBody)
    println(io, body.n_points, "-points ", typeof(body), ":")
    println(io, "  Number of bonds / 1-NI: ", body.n_bonds)
    return nothing
end

create_simmodel(mat::BondBasedMaterial, pc::PointCloud) = BondBasedBody(mat, pc)

function compute_forcedensity!(body::BondBasedBody, mat::BondBasedMaterial)
    body.b_int .= 0.0
    body.n_active_family_members .= 0
    @threads :static for _ in 1:(body.n_threads)
        tid = threadid()
        for current_bond in body.owned_bonds[tid]
            (i, j, ξ, failureflag) = body.bond_data[current_bond]
            ϑx = body.position[1, j] - body.position[1, i]
            ϑy = body.position[2, j] - body.position[2, i]
            ϑz = body.position[3, j] - body.position[3, i]
            η = sqrt(ϑx * ϑx + ϑy * ϑy + ϑz * ϑz)
            ε = (η - ξ) / ξ
            temp = body.bond_failure[current_bond] * mat.bc * body.volumes[j] * ε / η
            body.b_int[1, i, tid] += temp * ϑx
            body.b_int[2, i, tid] += temp * ϑy
            body.b_int[3, i, tid] += temp * ϑz
            body.b_int[1, j, tid] -= temp * ϑx
            body.b_int[2, j, tid] -= temp * ϑy
            body.b_int[3, j, tid] -= temp * ϑz
            if ε > mat.εc
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
