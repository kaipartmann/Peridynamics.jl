"""
    NOSBMaterial(; maxdmg, maxjacobi, corr)

A material type used to assign the material of a [`Body`](@ref) with the local continuum
consistent (correspondence) formulation of non-ordinary state-based peridynamics.

# Keywords
- `maxdmg::Float64`: Maximum value of damage a point is allowed to obtain. If this value is
    exceeded, all bonds of that point are broken because the deformation gradient would then
    possibly contain `NaN` values.
    (default: `0.95`)
- `maxjacobi::Float64`: Maximum value of the Jacobi determinant. If this value is exceeded,
    all bonds of that point are broken.
    (default: `1.03`)
- `corr::Float64`: Correction factor used for zero-energy mode stabilization. The
    stabilization algorithm of Silling (2017) is used.
    (default: `100.0`)

!!! note "Stability of fracture simulations"
    This formulation is known to be not suitable for fracture simultations without
    stabilization of the zero-energy modes. Therefore be careful when doing fracture
    simulations and try out different paremeters for `maxdmg`, `maxjacobi`, and `corr`.

# Examples

```julia-repl
julia> mat = NOSBMaterial()
NOSBMaterial(maxdmg=0.95, maxjacobi=1.03, corr=100.0)
```

---

```julia
NOSBMaterial
```

Material type for the local continuum consistent (correspondence) formulation of
non-ordinary state-based peridynamics.

# Fields
- `maxdmg::Float64`: Maximum value of damage a point is allowed to obtain. See the
    constructor docs for more informations.
- `maxjacobi::Float64`: Maximum value of the Jacobi determinant. See the constructor docs
    for more informations.
- `corr::Float64`: Correction factor used for zero-energy mode stabilization. See the
    constructor docs for more informations.

# Allowed material parameters
When using [`material!`](@ref) on a [`Body`](@ref) with `NOSBMaterial`, then the following
parameters are allowed:
- `horizon::Float64`: Radius of point interactions
- `rho::Float64`: Density
- `E::Float64`: Young's modulus
- `nu::Float64`: Poisson's ratio
- `Gc::Float64`: Critical energy release rate
- `epsilon_c::Float64`: Critical strain

# Allowed export fields
When specifying the `fields` keyword of [`Job`](@ref) for a [`Body`](@ref) with
`NOSBMaterial`, the following fields are allowed:
- `position::Matrix{Float64}`: Position of each point
- `displacement::Matrix{Float64}`: Displacement of each point
- `velocity::Matrix{Float64}`: Velocity of each point
- `velocity_half::Matrix{Float64}`: Velocity parameter for Verlet time solver
- `acceleration::Matrix{Float64}`: Acceleration of each point
- `b_int::Matrix{Float64}`: Internal force density of each point
- `b_ext::Matrix{Float64}`: External force density of each point
- `damage::Vector{Float64}`: Damage of each point
- `n_active_bonds::Vector{Int}`: Number of intact bonds of each point
"""
Base.@kwdef struct NOSBMaterial <: AbstractBondSystemMaterial{NoCorrection}
    maxdmg::Float64 = 0.95
    maxjacobi::Float64 = 1.03
    corr::Float64 = 100.0
end

function Base.show(io::IO, @nospecialize(mat::NOSBMaterial))
    print(io, typeof(mat))
    print(io, msg_fields_in_brackets(mat))
    return nothing
end

struct NOSBPointParameters <: AbstractPointParameters
    δ::Float64
    rho::Float64
    E::Float64
    nu::Float64
    G::Float64
    K::Float64
    λ::Float64
    μ::Float64
    Gc::Float64
    εc::Float64
    bc::Float64
end

function NOSBPointParameters(::NOSBMaterial, p::Dict{Symbol,Any})
    δ = get_horizon(p)
    rho = get_density(p)
    E, nu, G, K, λ, μ = get_elastic_params(p)
    Gc, εc = get_frac_params(p, δ, K)
    bc = 18 * K / (π * δ^4) # bond constant
    return NOSBPointParameters(δ, rho, E, nu, G, K, λ, μ, Gc, εc, bc)
end

@params NOSBMaterial NOSBPointParameters

struct NOSBVerletStorage <: AbstractStorage
    position::Matrix{Float64}
    displacement::Matrix{Float64}
    velocity::Matrix{Float64}
    velocity_half::Matrix{Float64}
    acceleration::Matrix{Float64}
    b_int::Matrix{Float64}
    b_ext::Matrix{Float64}
    damage::Vector{Float64}
    bond_active::Vector{Bool}
    n_active_bonds::Vector{Int}
    rotation::Matrix{Float64}
    left_stretch::Matrix{Float64}
end

function NOSBVerletStorage(::NOSBMaterial, ::VelocityVerlet, system::BondSystem, ch)
    n_loc_points = length(ch.loc_points)
    position = copy(system.position)
    n_all_points = size(position, 2)
    displacement = zeros(3, n_loc_points)
    velocity = zeros(3, n_all_points)
    velocity_half = zeros(3, n_loc_points)
    acceleration = zeros(3, n_loc_points)
    b_int = zeros(3, length(ch.point_ids))
    b_ext = zeros(3, n_loc_points)
    damage = zeros(n_loc_points)
    bond_active = ones(Bool, length(system.bonds))
    n_active_bonds = copy(system.n_neighbors)
    rotation = zeros(9, n_loc_points)
    rotation[[1, 5, 9], :] .= 1.0
    left_stretch = zeros(9, n_loc_points)
    left_stretch[[1, 5, 9], :] .= 1.0
    s = NOSBVerletStorage(position, displacement, velocity, velocity_half, acceleration,
                          b_int, b_ext, damage, bond_active, n_active_bonds, rotation,
                          left_stretch)
    return s
end

@storage NOSBMaterial VelocityVerlet NOSBVerletStorage

@loc_to_halo_fields NOSBVerletStorage :position :velocity
@halo_to_loc_fields NOSBVerletStorage :b_int

struct NOSBRelaxationStorage <: AbstractStorage
    position::Matrix{Float64}
    displacement::Matrix{Float64}
    velocity::Matrix{Float64}
    velocity_half::Matrix{Float64}
    velocity_half_old::Matrix{Float64}
    b_int::Matrix{Float64}
    b_int_old::Matrix{Float64}
    b_ext::Matrix{Float64}
    density_matrix::Matrix{Float64}
    damage::Vector{Float64}
    bond_active::Vector{Bool}
    n_active_bonds::Vector{Int}
    rotation::Matrix{Float64}
    left_stretch::Matrix{Float64}
end

function NOSBRelaxationStorage(::NOSBMaterial, ::DynamicRelaxation, system::BondSystem, ch)
    n_loc_points = length(ch.loc_points)
    position = copy(system.position)
    n_all_points = size(position, 2)
    displacement = zeros(3, n_loc_points)
    velocity = zeros(3, n_all_points)
    velocity_half = zeros(3, n_loc_points)
    velocity_half_old = zeros(3, n_loc_points)
    b_int = zeros(3, length(ch.point_ids))
    b_int_old = zeros(3, n_loc_points)
    b_ext = zeros(3, n_loc_points)
    density_matrix = zeros(3, n_loc_points)
    damage = zeros(n_loc_points)
    bond_active = ones(Bool, length(system.bonds))
    n_active_bonds = copy(system.n_neighbors)
    rotation = zeros(9, n_loc_points)
    rotation[[1, 5, 9], :] .= 1.0
    left_stretch = zeros(9, n_loc_points)
    left_stretch[[1, 5, 9], :] .= 1.0
    s = NOSBRelaxationStorage(position, displacement, velocity, velocity_half,
                              velocity_half_old, b_int, b_int_old, b_ext, density_matrix,
                              damage, bond_active, n_active_bonds, rotation, left_stretch)
    return s
end

@storage NOSBMaterial DynamicRelaxation NOSBRelaxationStorage

@loc_to_halo_fields NOSBRelaxationStorage :position :velocity
@halo_to_loc_fields NOSBRelaxationStorage :b_int

const NOSBStorage = Union{NOSBVerletStorage,NOSBRelaxationStorage}

function calc_force_density!(chunk::AbstractBodyChunk{S,M},
                             Δt::Float64) where {S<:BondSystem,M<:NOSBMaterial}
    (; system, mat, paramsetup, storage) = chunk
    storage.b_int .= 0
    storage.n_active_bonds .= 0
    for point_id in each_point_idx(chunk)
        force_density_point!(storage, system, mat, paramsetup, Δt, point_id)
    end
    return nothing
end

function force_density_point!(storage::NOSBStorage, system::BondSystem, mat::NOSBMaterial,
                              paramhandler::AbstractParameterHandler, Δt::Float64, i::Int)
    params = get_params(paramhandler, i)
    force_density_point!(storage, system, mat, params, Δt, i)
    return nothing
end

function force_density_point!(storage::NOSBStorage, system::BondSystem, mat::NOSBMaterial,
                              params::NOSBPointParameters, Δt::Float64, i::Int)
    F, Ḟ, Kinv, ω0 = calc_deformation_gradient(storage, system, mat, params, i)
    if storage.damage[i] > mat.maxdmg || containsnan(F)
        kill_point!(storage, system, i)
        return nothing
    end
    σ = calc_cauchy_stress(mat, params, F)
    calc_rod_and_rotation!(storage, F, Ḟ, Δt, i)
    R = get_rotation(storage, i)
    T = R * σ * R'
    P = det(F) * T * inv(F)'
    PKinv = P * Kinv
    for bond_id in each_bond_idx(system, i)
        bond = system.bonds[bond_id]
        j, L = bond.neighbor, bond.length
        ΔXij = get_coordinates_diff(system, i, j)
        Δxij = get_coordinates_diff(storage, i, j)
        l = norm(Δxij)
        ε = (l - L) / L
        stretch_based_failure!(storage, system, bond, params, ε, i, bond_id)

        # stabilization
        ωij = influence_function(mat, params, L) * storage.bond_active[bond_id]
        Tij = mat.corr .* params.bc * ωij / ω0 .* (Δxij .- F * ΔXij)

        # update of force density
        tij = ωij * PKinv * ΔXij + Tij
        update_add_b_int!(storage, i, tij .* system.volume[j])
        update_add_b_int!(storage, j, -tij .* system.volume[i])
    end
    return nothing
end

@inline function influence_function(::NOSBMaterial, params::NOSBPointParameters, L::Float64)
    return params.δ / L
end

function calc_deformation_gradient(storage::NOSBStorage, system::BondSystem,
                                   mat::NOSBMaterial, params::NOSBPointParameters, i::Int)
    K = zeros(SMatrix{3,3})
    _F = zeros(SMatrix{3,3})
    _Fdot = zeros(SMatrix{3,3})
    ω0 = 0.0
    for bond_id in each_bond_idx(system, i)
        bond = system.bonds[bond_id]
        j, L = bond.neighbor, bond.length
        ΔXij = get_coordinates_diff(system, i, j)
        Δxij = get_coordinates_diff(storage, i, j)
        Δvij = get_diff(storage.velocity, i, j)
        ωij = influence_function(mat, params, L) * storage.bond_active[bond_id]
        ω0 += ωij
        temp = ωij * system.volume[j]
        K += temp * ΔXij * ΔXij'
        _F += temp * Δxij * ΔXij'
        _Fdot += temp * Δvij * ΔXij'
    end
    Kinv = inv(K)
    F = _F * Kinv
    Fdot = _Fdot * Kinv
    return F, Fdot, Kinv, ω0
end

# function calc_first_piola_stress(F::SMatrix{3,3}, mat::NOSBMaterial,
#                                  params::NOSBPointParameters)
#     J = det(F)
#     J < eps() && return zero(SMatrix{3,3})
#     J > mat.maxjacobi && return zero(SMatrix{3,3})
#     C = F' * F
#     Cinv = inv(C)
#     S = params.G .* (I - 1 / 3 .* tr(C) .* Cinv) .* J^(-2 / 3) .+
#         params.K / 4 .* (J^2 - J^(-2)) .* Cinv
#     P = F * S
#     return P
# end

function calc_cauchy_stress(mat::NOSBMaterial, params::NOSBPointParameters, F::SMatrix{3,3})
    # ---- Neo-Hookean model
    #      Ψ = G/2 (I₁ - 3) + K/2 (J - 1)^2
    #      P = ∂Ψ/∂F = ∂Ψ/∂I₁ ∂I₁/∂F + ∂Ψ/∂J ∂J/∂F
    #        = G F + K (J - 1) J F^(-T)
    #      σ = 1/J P F^T
    # J = det(F)
    # J < eps() && return zero(SMatrix{3,3})
    # J > mat.maxjacobi && return zero(SMatrix{3,3})
    # B = F * F'
    # σ = params.G / J .* B + params.K * (J - 1) .* I
    # ----

    # ---- Combined model with a Neo-Hookean type deviatoric response and a nonlinear
    #      volumetric response:
    # J = det(F)
    # J < eps() && return zero(SMatrix{3,3})
    # J > mat.maxjacobi && return zero(SMatrix{3,3})
    # C = F' * F
    # Cinv = inv(C)
    # S = params.G .* (I - 1 / 3 .* tr(C) .* Cinv) .* J^(-2 / 3) .+
    #     params.K / 4 .* (J^2 - J^(-2)) .* Cinv
    # P = F * S
    # σ = 1/J .* P * F'
    # ----

    # ---- Saint Venant-Kirchhoff, taken from
    #      https://doi.org/10.1016/j.jfluidstructs.2021.103312
    J = det(F)
    J < eps() && return zero(SMatrix{3,3})
    J > mat.maxjacobi && return zero(SMatrix{3,3})
    E = 0.5 .* (F' * F - I)
    S = params.λ * tr(E) * I + 2 * params.μ * E
    P = F * S
    σ = 1/J .* P * F'
    # ----

    return σ
end

@inline function get_left_stretch(storage::NOSBStorage, i::Int)
    V = storage.left_stretch
    _V = SMatrix{3,3}(V[1,i], V[2,i], V[3,i], V[4,i], V[5,i], V[6,i], V[7,i], V[8,i], V[9,i])
    return _V
end

@inline function get_rotation(storage::NOSBStorage, i::Int)
    R = storage.rotation
    _R = SMatrix{3,3}(R[1,i], R[2,i], R[3,i], R[4,i], R[5,i], R[6,i], R[7,i], R[8,i], R[9,i])
    return _R
end

@inline function update_left_stretch!(storage::NOSBStorage, i::Int, Vₙ₊₁::SMatrix{3,3})
    Vₙ = storage.left_stretch
    Vₙ[1,i] = Vₙ₊₁[1,1]
    Vₙ[2,i] = Vₙ₊₁[1,2]
    Vₙ[3,i] = Vₙ₊₁[1,3]
    Vₙ[4,i] = Vₙ₊₁[2,1]
    Vₙ[5,i] = Vₙ₊₁[2,2]
    Vₙ[6,i] = Vₙ₊₁[2,3]
    Vₙ[7,i] = Vₙ₊₁[3,1]
    Vₙ[8,i] = Vₙ₊₁[3,2]
    Vₙ[9,i] = Vₙ₊₁[3,3]
    return nothing
end

@inline function update_rotation!(storage::NOSBStorage, i::Int, Rₙ₊₁::SMatrix{3,3})
    Rₙ = storage.rotation
    Rₙ[1,i] = Rₙ₊₁[1,1]
    Rₙ[2,i] = Rₙ₊₁[1,2]
    Rₙ[3,i] = Rₙ₊₁[1,3]
    Rₙ[4,i] = Rₙ₊₁[2,1]
    Rₙ[5,i] = Rₙ₊₁[2,2]
    Rₙ[6,i] = Rₙ₊₁[2,3]
    Rₙ[7,i] = Rₙ₊₁[3,1]
    Rₙ[8,i] = Rₙ₊₁[3,2]
    Rₙ[9,i] = Rₙ₊₁[3,3]
    return nothing
end

"""
    calc_rod_and_rotation!()

Calculates the rate of deformation and the rotation tensor needed for the kinematic
computations described
"""
function calc_rod_and_rotation!(storage, F, Ḟ, Δt, i)
    #(storage::NOSBStorage, system::BondSystem, mat::NOSBMaterial, params::NOSBPointParameters, i::Int)

    # inverse of the deformation gradient
    F⁻¹ = inv(F)

    # Eulerian velocity gradient [FT87, eq. (3)]
    L = Ḟ * F⁻¹

    # rate-of-deformation tensor D
    D = 0.5 .* (L + L')

    # spin tensor W
    W = 0.5 .* (L - L')

    # left stretch V
    # 1 -> 1,1
    # 2 -> 1,2
    # 3 -> 1,3
    # 4 -> 2,1
    # 5 -> 2,2
    # 6 -> 2,3
    # 7 -> 3,1
    # 8 -> 3,2
    # 9 -> 3,3
    V = get_left_stretch(storage, i)

    # vector z [FT87, eq. (13)]
    z_x = - V[1,3] * D[2,1] - V[2,3] * D[2,2] -
            V[3,3] * D[2,3] + V[1,2] * D[3,1] +
            V[2,2] * D[3,2] + V[3,2] * D[3,3]
    z_y = V[1,3] * D[1,1] + V[2,3] * D[1,2] +
          V[3,3] * D[1,3] - V[1,1] * D[3,1] -
          V[2,1] * D[3,2] - V[3,1] * D[3,3]
    z_z = - V[1,2] * D[1,1] - V[2,2] * D[1,2] -
            V[3,2] * D[1,3] + V[1,1] * D[2,1] +
            V[2,1] * D[2,2] + V[3,1] * D[2,3]
    z = SVector{3}(z_x, z_y, z_z)

    # w = -1/2 * \epsilon_{ijk} * W_{jk}  [FT87, eq. (11)]
    w = 0.5 .* SVector{3}(W[3,1] - W[2,2], W[1,2] - W[2,3], W[1,3] - W[1,1])

    # ω = w + (I * tr(V) - V)^(-1) * z [FT87, eq. (12)]
    ω = w + inv(I * tr(V) - V) * z

    # Ω [FT87, eq. (10)]
    Ωtens = SMatrix{3,3}(0.0, -ω[3], ω[2], ω[3], 0.0, -ω[1], ω[2], ω[1], 0.0)
    Ωtens² = Ωtens * Ωtens
    Ω² = dot(ω, ω)
    Ω = sqrt(Ω²)

    # compute Q with [FT87, eq. (44)]
    if Ω² > 1e-30 # avoid a potential divide-by-zero
        fac1 = sin(Δt * Ω) / Ω
        fac2 = -(1.0 - cos(Δt * Ω)) / Ω²
        Q = I + fac1 .* Ωtens + fac2 .* Ωtens²
    else
        Q = SMatrix{3,3}(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)
    end

    # compute Rotation of new step [FT87, eq. (36)]
    R = get_rotation(storage, i)
    Rₙ₊₁ = Q * R

    # compute step 4 of [FT87]
    V̇ = L * V - V * Ωtens

    # compute step 5 of [FT87]
    Vₙ₊₁ = V + Δt * V̇

    # update rotation and left stretch
    update_rotation!(storage, i, Rₙ₊₁)
    update_left_stretch!(storage, i, Vₙ₊₁)

    # compute step 6 of [FT87]
    # d = R' * D * R

    return nothing
end

function containsnan(K::T) where {T<:AbstractArray}
    @simd for i in eachindex(K)
        isnan(K[i]) && return true
    end
    return false
end

function kill_point!(s::AbstractStorage, bd::BondSystem, i::Int)
    s.bond_active[each_bond_idx(bd, i)] .= false
    s.n_active_bonds[i] = 0
    return nothing
end
