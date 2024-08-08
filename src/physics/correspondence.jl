"""
    NOSBMaterial(; maxdmg, corr)

A material type used to assign the material of a [`Body`](@ref) with the local continuum
consistent (correspondence) formulation of non-ordinary state-based peridynamics.

# Keywords
- `maxdmg::Float64`: Maximum value of damage a point is allowed to obtain. If this value is
    exceeded, all bonds of that point are broken because the deformation gradient would then
    possibly contain `NaN` values.
    (default: `0.95`)
- `corr::Float64`: Correction factor used for zero-energy mode stabilization. The
    stabilization algorithm of Silling (2017) is used.
    (default: `100.0`)

!!! note "Stability of fracture simulations"
    This formulation is known to be not suitable for fracture simultations without
    stabilization of the zero-energy modes. Therefore be careful when doing fracture
    simulations and try out different paremeters for `maxdmg`, and `corr`.

# Examples

```julia-repl
julia> mat = NOSBMaterial()
NOSBMaterial(maxdmg=0.95, corr=100.0)
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
    rate_of_deformation::Matrix{Float64}
    stress::Matrix{Float64}
    stress_rotated::Matrix{Float64}
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
    rate_of_deformation = zeros(9, n_loc_points)
    rate_of_deformation[[1, 5, 9], :] .= 1.0
    stress = zeros(9, n_loc_points)
    stress_rotated = zeros(9, n_loc_points)
    s = NOSBVerletStorage(position, displacement, velocity, velocity_half, acceleration,
                          b_int, b_ext, damage, bond_active, n_active_bonds, rotation,
                          left_stretch, rate_of_deformation, stress, stress_rotated)
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
    rate_of_deformation::Matrix{Float64}
    stress::Matrix{Float64}
    stress_rotated::Matrix{Float64}
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
    rate_of_deformation = zeros(9, n_loc_points)
    rate_of_deformation[[1, 5, 9], :] .= 1.0
    stress = zeros(9, n_loc_points)
    stress_rotated = zeros(9, n_loc_points)
    s = NOSBRelaxationStorage(position, displacement, velocity, velocity_half,
                              velocity_half_old, b_int, b_int_old, b_ext, density_matrix,
                              damage, bond_active, n_active_bonds, rotation, left_stretch,
                              rate_of_deformation, stress, stress_rotated)
    return s
end

@storage NOSBMaterial DynamicRelaxation NOSBRelaxationStorage

@loc_to_halo_fields NOSBRelaxationStorage :position :velocity
@halo_to_loc_fields NOSBRelaxationStorage :b_int

const NOSBStorage = Union{NOSBVerletStorage,NOSBRelaxationStorage}

function point_data_fields(::Type{S}) where {S<:NOSBStorage}
    fields = (:position, :displacement, :velocity, :velocity_half, :acceleration, :b_int,
              :b_ext, :damage, :n_active_bonds, :left_stretch, :rate_of_deformation,
              :rotation, :stress, :stress_rotated)
    return fields
end

point_data_field(s::NOSBStorage, ::Val{:left_stretch}) = getfield(s, :left_stretch)
point_data_field(s::NOSBStorage, ::Val{:rate_of_deformation}) = getfield(s, :rate_of_deformation)
point_data_field(s::NOSBStorage, ::Val{:rotation}) = getfield(s, :rotation)
point_data_field(s::NOSBStorage, ::Val{:stress}) = getfield(s, :stress)
point_data_field(s::NOSBStorage, ::Val{:stress_rotated}) = getfield(s, :stress_rotated)

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
    calc_rod_and_rotation!(storage, F, Ḟ, Δt, i)
    σ = calc_cauchy_stress(storage, mat, params, F, Δt, i)
    R = get_rotation(storage, i)
    T = R * σ * R'
    # T = σ
    update_tensor!(storage.stress_rotated, i, T)

    # F, Kinv = calc_stress_test!(storage, system, mat, params, Δt, i)
    # T = get_tensor(storage.stress_rotated, i)

    P = det(F) * T * inv(F)'
    PKinv = P * Kinv
    # C_zec = mat.corr * 18 * params.K / (π * params.δ^4)
    for bond_id in each_bond_idx(system, i)
        bond = system.bonds[bond_id]
        j, L = bond.neighbor, bond.length
        ΔXij = get_coordinates_diff(system, i, j)
        Δxij = get_coordinates_diff(storage, i, j)
        l = norm(Δxij)
        ε = (l - L) / L
        stretch_based_failure!(storage, system, bond, params, ε, i, bond_id)

        ωij = influence_function(mat, params, L) * storage.bond_active[bond_id]

        # stabilization Silling
        Tij = mat.corr .* params.bc * ωij / ω0 .* (Δxij .- F * ΔXij)

        # Stabilization peridigm / Littlewood
        # xi = get_coordinates(storage, i)
        # xj = get_coordinates(storage, j)
        # exij = SVector{3}(xi[1] + F[1] * ΔXij[1] + F[2] * ΔXij[2] + F[3] * ΔXij[3],
        #                   xi[2] + F[4] * ΔXij[1] + F[5] * ΔXij[2] + F[6] * ΔXij[3],
        #                   xi[3] + F[7] * ΔXij[1] + F[8] * ΔXij[2] + F[9] * ΔXij[3])
        # h = exij - xj
        # Tij = storage.bond_active[bond_id] * C_zec * (-dot(h, Δxij) / L) * (1.0 / l) * Δxij

        # update of force density
        tij = ωij * PKinv * ΔXij + Tij
        update_add_b_int!(storage, i, tij .* system.volume[j])
        update_add_b_int!(storage, j, -tij .* system.volume[i])
    end

    # for bond_id in each_bond_idx(system, i)
    #     bond = system.bonds[bond_id]
    #     j, L = bond.neighbor, bond.length
    #     Δxij = get_coordinates_diff(storage, i, j)
    #     l = norm(Δxij)
    #     ε = (l - L) / L
    #     stretch_based_failure!(storage, system, bond, params, ε, i, bond_id)
    #     b_int = bond_failure(storage, bond_id) *
    #             params.bc * ε / l * system.volume[j] .* Δxij
    #     update_add_b_int!(storage, i, b_int)
    # end
    return nothing
end

@inline function influence_function(::NOSBMaterial, params::NOSBPointParameters, L::Float64)
    return 1 + params.δ / L
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
#     C = F' * F
#     Cinv = inv(C)
#     S = params.G .* (I - 1 / 3 .* tr(C) .* Cinv) .* J^(-2 / 3) .+
#         params.K / 4 .* (J^2 - J^(-2)) .* Cinv
#     P = F * S
#     return P
# end

function calc_cauchy_stress(storage::NOSBStorage, mat::NOSBMaterial,
                            params::NOSBPointParameters, F::SMatrix{3,3}, Δt::Float64,
                            i::Int)
    # ---- Neo-Hookean model
    #      Ψ = G/2 (I₁ - 3) + K/2 (J - 1)^2
    #      P = ∂Ψ/∂F = ∂Ψ/∂I₁ ∂I₁/∂F + ∂Ψ/∂J ∂J/∂F
    #        = G F + K (J - 1) J F^(-T)
    #      σ = 1/J P F^T
    # J = det(F)
    # J < eps() && return zero(SMatrix{3,3})
    # B = F * F'
    # σ = params.G / J .* B + params.K * (J - 1) .* I
    # ----

    # ---- Combined model with a Neo-Hookean type deviatoric response and a nonlinear
    #      volumetric response:
    J = det(F)
    J < eps() && return zero(SMatrix{3,3})
    C = F' * F
    Cinv = inv(C)
    S = params.G .* (I - 1 / 3 .* tr(C) .* Cinv) .* J^(-2 / 3) .+
        params.K / 4 .* (J^2 - J^(-2)) .* Cinv
    P = F * S
    σ = 1/J .* P * F'
    update_tensor!(storage.stress, i, σ)
    # ----

    # ---- Saint Venant-Kirchhoff, taken from
    #      https://doi.org/10.1016/j.jfluidstructs.2021.103312
    # J = det(F)
    # J < eps() && return zero(SMatrix{3,3})
    # E = 0.5 .* (F' * F - I)
    # S = params.λ * tr(E) * I + 2 * params.μ * E
    # P = F * S
    # σ = 1/J .* P * F'
    # update_tensor!(storage.stress, i, σ)
    # ----

    # ---- taken from peridigm - elastic correspondence
    # if storage.damage[i] > mat.maxdmg
    #     return get_tensor(storage.stress, i)
    # end
    # d = get_tensor(storage.rate_of_deformation, i)
    # strain_inc = d * Δt

    # # v1
    # dilatation_inc = tr(strain_inc)
    # dev_strain_inc = strain_inc - dilatation_inc / 3 * I
    # σ_old = get_tensor(storage.stress, i)
    # σ = σ_old + 2 * params.G * dev_strain_inc + params.K * dilatation_inc * I

    # # v2
    # dil_inc = strain_inc[1,1] + strain_inc[2,2] + strain_inc[3,3]
    # dev_strain_inc = SMatrix{3,3}(strain_inc[1] - dil_inc / 3,
    #                               strain_inc[2],
    #                               strain_inc[3],
    #                               strain_inc[4],
    #                               strain_inc[5] - dil_inc / 3,
    #                               strain_inc[6],
    #                               strain_inc[7],
    #                               strain_inc[8],
    #                               strain_inc[9] - dil_inc / 3)
    # σₙ = get_tensor(storage.stress, i)
    # σₙ₊₁ = SMatrix{3,3}(σₙ[1] + 2 * params.G * dev_strain_inc[1] + params.K * dil_inc,
    #                     σₙ[2] + 2 * params.G * dev_strain_inc[2],
    #                     σₙ[3] + 2 * params.G * dev_strain_inc[3],
    #                     σₙ[4] + 2 * params.G * dev_strain_inc[4],
    #                     σₙ[5] + 2 * params.G * dev_strain_inc[5] + params.K * dil_inc,
    #                     σₙ[6] + 2 * params.G * dev_strain_inc[6],
    #                     σₙ[7] + 2 * params.G * dev_strain_inc[7],
    #                     σₙ[8] + 2 * params.G * dev_strain_inc[8],
    #                     σₙ[9] + 2 * params.G * dev_strain_inc[9] + params.K * dil_inc)
    # # σ = σₙ₊₁
    # if !(σ ≈ σₙ₊₁)
    #     @show d
    #     @show
    #     @show σ
    #     @show σₙ
    #     @show σₙ₊₁
    #     error()
    # end
    # update_tensor!(storage.stress, i, σ)
    # ----

    return σ
end

@inline function get_tensor(T::AbstractMatrix, i::Int)
    tensor = SMatrix{3,3}(T[1,i], T[2,i], T[3,i], T[4,i], T[5,i], T[6,i], T[7,i], T[8,i],
                          T[9,i])
    return tensor
end

@inline function get_left_stretch(storage::NOSBStorage, i::Int)
    return get_tensor(storage.left_stretch, i)
end

@inline function get_rotation(storage::NOSBStorage, i::Int)
    return get_tensor(storage.rotation, i)
end

@inline function update_tensor!(Tₙ::AbstractMatrix, i::Int, Tₙ₊₁::SMatrix{3,3})
    Tₙ[1,i] = Tₙ₊₁[1,1]
    Tₙ[2,i] = Tₙ₊₁[1,2]
    Tₙ[3,i] = Tₙ₊₁[1,3]
    Tₙ[4,i] = Tₙ₊₁[2,1]
    Tₙ[5,i] = Tₙ₊₁[2,2]
    Tₙ[6,i] = Tₙ₊₁[2,3]
    Tₙ[7,i] = Tₙ₊₁[3,1]
    Tₙ[8,i] = Tₙ₊₁[3,2]
    Tₙ[9,i] = Tₙ₊₁[3,3]
    return nothing
end

@inline function update_left_stretch!(storage::NOSBStorage, i::Int, V::SMatrix{3,3})
    update_tensor!(storage.left_stretch, i, V)
    return nothing
end

@inline function update_rotation!(storage::NOSBStorage, i::Int, R::SMatrix{3,3})
    update_tensor!(storage.rotation, i, R)
    return nothing
end

"""
    calc_rod_and_rotation!()

Calculates the rate of deformation and the rotation tensor needed for the kinematic
computations described
"""
function calc_rod_and_rotation!(storage, F, Ḟ, Δt, i)
    # inverse of the deformation gradient
    F⁻¹ = inv(F)

    # Eulerian velocity gradient [FT87, eq. (3)]
    L = Ḟ * F⁻¹

    # rate-of-deformation tensor D
    D = 0.5 .* (L + L')

    # spin tensor W
    W = 0.5 .* (L - L')

    # left stretch V
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
    w = 0.5 .* SVector{3}(W[3,2] - W[2,3], W[1,3] - W[3,1], W[2,1] - W[1,2])

    # ω = w + (I * tr(V) - V)^(-1) * z [FT87, eq. (12)]
    ω = w + inv(I * tr(V) - V) * z

    # Ω [FT87, eq. (10)]
    Ωtens = SMatrix{3,3}(0.0, -ω[3], ω[2], ω[3], 0.0, -ω[1], -ω[2], ω[1], 0.0)
    Ω² = dot(ω, ω)
    Ω = sqrt(Ω²)

    # compute Q with [FT87, eq. (44)]
    if Ω² > 1e-30 # avoid a potential divide-by-zero
        fac1 = sin(Δt * Ω) / Ω
        fac2 = -(1.0 - cos(Δt * Ω)) / Ω²
        Ωtens² = Ωtens * Ωtens
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
    d = Rₙ₊₁' * D * Rₙ₊₁
    update_tensor!(storage.rate_of_deformation, i, d)

    if containsnan(d) && Threads.threadid() == 1
        @show F
        @show Ḟ
        @show F⁻¹
        @show L
        @show V
        @show Rₙ₊₁
        @show d
        error()
    end
    return nothing
end

# function calc_stress_test!(storage::NOSBStorage, system::BondSystem, mat::NOSBMaterial,
#                            params::NOSBPointParameters, Δt::Float64, a::Int)
#     𝐊 = @MMatrix zeros(3,3)
#     inv𝐊 = @MMatrix zeros(3,3)
#     𝐅 = @MMatrix zeros(3,3)
#     𝐅₁ = @MMatrix zeros(3,3)
#     𝐅dot = @MMatrix zeros(3,3)
#     𝐅dot₁ = @MMatrix zeros(3,3)
#     𝐋 = @MMatrix zeros(3,3)
#     RoD = @MMatrix zeros(3,3)
#     Spin = @MMatrix zeros(3,3)
#     left_stretch = @MMatrix zeros(3,3)
#     z = @MVector zeros(3)
#     w = @MVector zeros(3)
#     𝛀Tens = @MMatrix zeros(3,3)
#     Qmatrix = @MMatrix zeros(3,3)
#     rot_tens_ = @MMatrix zeros(3,3)
#     𝛔_unrot = @MMatrix zeros(3,3)
#     tempVec = @MMatrix zeros(3,3)
#     if storage.damage[a] < mat.maxdmg
#         𝐊 .= 0.
#         𝐅₁ .= 0.
#         𝐅dot₁ .= 0.
#         for iInt in system.bond_ids[a]
#             bond = system.bonds[iInt]
#             i, L = bond.neighbor, bond.length
#             Ξaix = system.position[1, i] - system.position[1, a]
#             Ξaiy = system.position[2, i] - system.position[2, a]
#             Ξaiz = system.position[3, i] - system.position[3, a]
#             ξaix = storage.position[1, i] - storage.position[1, a]
#             ξaiy = storage.position[2, i] - storage.position[2, a]
#             ξaiz = storage.position[3, i] - storage.position[3, a]
#             vaix = storage.velocity[1, i] - storage.velocity[1, a]
#             vaiy = storage.velocity[2, i] - storage.velocity[2, a]
#             vaiz = storage.velocity[3, i] - storage.velocity[3, a]
#             temp = storage.bond_active[iInt] * influence_function(mat, params, L) *
#                    system.volume[i]
#             𝐊[1,1] += temp * Ξaix * Ξaix
#             𝐊[1,2] += temp * Ξaix * Ξaiy
#             𝐊[1,3] += temp * Ξaix * Ξaiz
#             𝐊[2,1] += temp * Ξaiy * Ξaix
#             𝐊[2,2] += temp * Ξaiy * Ξaiy
#             𝐊[2,3] += temp * Ξaiy * Ξaiz
#             𝐊[3,1] += temp * Ξaiz * Ξaix
#             𝐊[3,2] += temp * Ξaiz * Ξaiy
#             𝐊[3,3] += temp * Ξaiz * Ξaiz
#             𝐅₁[1,1] += temp * ξaix * Ξaix
#             𝐅₁[1,2] += temp * ξaix * Ξaiy
#             𝐅₁[1,3] += temp * ξaix * Ξaiz
#             𝐅₁[2,1] += temp * ξaiy * Ξaix
#             𝐅₁[2,2] += temp * ξaiy * Ξaiy
#             𝐅₁[2,3] += temp * ξaiy * Ξaiz
#             𝐅₁[3,1] += temp * ξaiz * Ξaix
#             𝐅₁[3,2] += temp * ξaiz * Ξaiy
#             𝐅₁[3,3] += temp * ξaiz * Ξaiz
#             𝐅dot₁[1,1] += temp * vaix * Ξaix
#             𝐅dot₁[1,2] += temp * vaix * Ξaiy
#             𝐅dot₁[1,3] += temp * vaix * Ξaiz
#             𝐅dot₁[2,1] += temp * vaiy * Ξaix
#             𝐅dot₁[2,2] += temp * vaiy * Ξaiy
#             𝐅dot₁[2,3] += temp * vaiy * Ξaiz
#             𝐅dot₁[3,1] += temp * vaiz * Ξaix
#             𝐅dot₁[3,2] += temp * vaiz * Ξaiy
#             𝐅dot₁[3,3] += temp * vaiz * Ξaiz
#         end
#         inv𝐊 = inv(𝐊)
#         𝐅 = 𝐅₁ * inv𝐊
#         𝐅dot = 𝐅dot₁ * inv𝐊
#         𝐋 = 𝐅dot * inv(𝐅)
#         RoD = 1/2 * (𝐋 + transpose(𝐋))
#         Spin = 1/2 * (𝐋 - transpose(𝐋))
#         left_stretch[1,1] = storage.left_stretch[1,a]
#         left_stretch[1,2] = storage.left_stretch[2,a]
#         left_stretch[1,3] = storage.left_stretch[3,a]
#         left_stretch[2,1] = storage.left_stretch[4,a]
#         left_stretch[2,2] = storage.left_stretch[5,a]
#         left_stretch[2,3] = storage.left_stretch[6,a]
#         left_stretch[3,1] = storage.left_stretch[7,a]
#         left_stretch[3,2] = storage.left_stretch[8,a]
#         left_stretch[3,3] = storage.left_stretch[9,a]
#         z[1] = - storage.left_stretch[3, a] * RoD[4] - storage.left_stretch[6, a] * RoD[5] -
#                 storage.left_stretch[9, a] * RoD[6] + storage.left_stretch[2, a] * RoD[7] +
#                 storage.left_stretch[5, a] * RoD[8] + storage.left_stretch[8, a] * RoD[9]
#         z[2] =   storage.left_stretch[3, a] * RoD[1] + storage.left_stretch[6, a] * RoD[2] +
#                 storage.left_stretch[9, a] * RoD[3] - storage.left_stretch[1, a] * RoD[7] -
#                 storage.left_stretch[4, a] * RoD[8] - storage.left_stretch[7, a] * RoD[9]
#         z[3] = - storage.left_stretch[2, a] * RoD[1] - storage.left_stretch[5, a] * RoD[2] -
#                 storage.left_stretch[8, a] * RoD[3] + storage.left_stretch[1, a] * RoD[4] +
#                 storage.left_stretch[4, a] * RoD[5] + storage.left_stretch[7, a] * RoD[6]
#         w[1] = 0.5 * (Spin[3,2] - Spin[2,3])
#         w[2] = 0.5 * (Spin[1,3] - Spin[3,1])
#         w[3] = 0.5 * (Spin[2,1] - Spin[1,2])
#         traceV = storage.left_stretch[1, a] + storage.left_stretch[5, a] + storage.left_stretch[9, a]
#         omega = w + inv(traceV * I - left_stretch) * z
#         𝛀Tens[1,1] = 0.0
#         𝛀Tens[1,2] = -omega[3]
#         𝛀Tens[1,3] = omega[2]
#         𝛀Tens[2,1] = omega[3]
#         𝛀Tens[2,2] = 0.0
#         𝛀Tens[2,3] = -omega[1]
#         𝛀Tens[3,1] = -omega[2]
#         𝛀Tens[3,2] = omega[1]
#         𝛀Tens[3,3] = 0.0
#         ΩSq = omega[1]^2 + omega[2]^2 + omega[3]^2
#         Ω = sqrt(ΩSq)
#         if ΩSq > 1e-16 && Ω !== Inf
#             scfac1 = sin(Δt * Ω) / Ω
#             scfac2 = -(1 - cos(Δt * Ω)) / ΩSq
#             𝛀TensSq = 𝛀Tens * 𝛀Tens
#             Qmatrix = I + scfac1 * 𝛀Tens + scfac2 * 𝛀TensSq
#         else
#             Qmatrix .= 0.
#             Qmatrix[1,1] = 1.0
#             Qmatrix[2,2] = 1.0
#             Qmatrix[3,3] = 1.0
#         end
#         rot_tens_[1,1] = storage.rotation[1,a]
#         rot_tens_[1,2] = storage.rotation[2,a]
#         rot_tens_[1,3] = storage.rotation[3,a]
#         rot_tens_[2,1] = storage.rotation[4,a]
#         rot_tens_[2,2] = storage.rotation[5,a]
#         rot_tens_[2,3] = storage.rotation[6,a]
#         rot_tens_[3,1] = storage.rotation[7,a]
#         rot_tens_[3,2] = storage.rotation[8,a]
#         rot_tens_[3,3] = storage.rotation[9,a]
#         rot_tens = Qmatrix * rot_tens_
#         storage.rotation[1, a] = rot_tens[1,1]
#         storage.rotation[2, a] = rot_tens[1,2]
#         storage.rotation[3, a] = rot_tens[1,3]
#         storage.rotation[4, a] = rot_tens[2,1]
#         storage.rotation[5, a] = rot_tens[2,2]
#         storage.rotation[6, a] = rot_tens[2,3]
#         storage.rotation[7, a] = rot_tens[3,1]
#         storage.rotation[8, a] = rot_tens[3,2]
#         storage.rotation[9, a] = rot_tens[3,3]
#         Vdot = 𝐋 * left_stretch - left_stretch * 𝛀Tens
#         storage.left_stretch[1,a] += Δt * Vdot[1,1]
#         storage.left_stretch[2,a] += Δt * Vdot[1,2]
#         storage.left_stretch[3,a] += Δt * Vdot[1,3]
#         storage.left_stretch[4,a] += Δt * Vdot[2,1]
#         storage.left_stretch[5,a] += Δt * Vdot[2,2]
#         storage.left_stretch[6,a] += Δt * Vdot[2,3]
#         storage.left_stretch[7,a] += Δt * Vdot[3,1]
#         storage.left_stretch[8,a] += Δt * Vdot[3,2]
#         storage.left_stretch[9,a] += Δt * Vdot[3,3]
#         tempVec = RoD * rot_tens
#         UnRotRoD = transpose(rot_tens) * tempVec
#         strainInc = UnRotRoD * Δt
#         deviatoricStrain = copy(strainInc)
#         dilatation = strainInc[1,1] + strainInc[2,2] + strainInc[3,3]
#         deviatoricStrain[1,1] -= dilatation/3
#         deviatoricStrain[2,2] -= dilatation/3
#         deviatoricStrain[3,3] -= dilatation/3
#         storage.stress[1,a] += deviatoricStrain[1,1] * 2 * params.G + params.K * dilatation
#         storage.stress[2,a] += deviatoricStrain[1,2] * 2 * params.G
#         storage.stress[3,a] += deviatoricStrain[1,3] * 2 * params.G
#         storage.stress[4,a] += deviatoricStrain[2,1] * 2 * params.G
#         storage.stress[5,a] += deviatoricStrain[2,2] * 2 * params.G + params.K * dilatation
#         storage.stress[6,a] += deviatoricStrain[2,3] * 2 * params.G
#         storage.stress[7,a] += deviatoricStrain[3,1] * 2 * params.G
#         storage.stress[8,a] += deviatoricStrain[3,2] * 2 * params.G
#         storage.stress[9,a] += deviatoricStrain[3,3] * 2 * params.G + params.K * dilatation
#         𝛔_unrot[1,1] = storage.stress[1,a]
#         𝛔_unrot[1,2] = storage.stress[2,a]
#         𝛔_unrot[1,3] = storage.stress[3,a]
#         𝛔_unrot[2,1] = storage.stress[4,a]
#         𝛔_unrot[2,2] = storage.stress[5,a]
#         𝛔_unrot[2,3] = storage.stress[6,a]
#         𝛔_unrot[3,1] = storage.stress[7,a]
#         𝛔_unrot[3,2] = storage.stress[8,a]
#         𝛔_unrot[3,3] = storage.stress[9,a]
#         tempVec = 𝛔_unrot * transpose(rot_tens)
#         𝛔 = rot_tens * tempVec
#         storage.stress_rotated[1, a] = 𝛔[1, 1]
#         storage.stress_rotated[2, a] = 𝛔[1, 2]
#         storage.stress_rotated[3, a] = 𝛔[1, 3]
#         storage.stress_rotated[4, a] = 𝛔[2, 1]
#         storage.stress_rotated[5, a] = 𝛔[2, 2]
#         storage.stress_rotated[6, a] = 𝛔[2, 3]
#         storage.stress_rotated[7, a] = 𝛔[3, 1]
#         storage.stress_rotated[8, a] = 𝛔[3, 2]
#         storage.stress_rotated[9, a] = 𝛔[3, 3]
#     else
#         for ii = 1:9
#             storage.stress_rotated[ii, a] = 0.
#             storage.left_stretch[ii, a] = 0.
#             storage.rotation[ii, a] = 0.
#         end
#     end
#     return 𝐅, inv𝐊
# end

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
