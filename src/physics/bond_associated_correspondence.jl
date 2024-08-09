Base.@kwdef struct BANOSBMaterial <: AbstractBondAssociatedSystemMaterial
    maxdmg::Float64 = 0.95
end

struct BANOSBPointParameters <: AbstractPointParameters
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

function BANOSBPointParameters(::BANOSBMaterial, p::Dict{Symbol,Any})
    δ = get_horizon(p)
    rho = get_density(p)
    E, nu, G, K, λ, μ = get_elastic_params(p)
    Gc, εc = get_frac_params(p, δ, K)
    bc = 18 * K / (π * δ^4) # bond constant
    return BANOSBPointParameters(δ, rho, E, nu, G, K, λ, μ, Gc, εc, bc)
end

@params BANOSBMaterial BANOSBPointParameters

struct BANOSBVerletStorage <: AbstractStorage
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

function BANOSBVerletStorage(::BANOSBMaterial, ::VelocityVerlet,
                             system::BondAssociatedSystem, ch)
    n_loc_points = length(ch.loc_points)
    position = copy(system.position)
    displacement = zeros(3, n_loc_points)
    velocity = zeros(3, size(position, 2))
    velocity_half = zeros(3, n_loc_points)
    acceleration = zeros(3, n_loc_points)
    b_int = zeros(3, length(ch.point_ids))
    b_ext = zeros(3, n_loc_points)
    damage = zeros(n_loc_points)
    bond_active = ones(Bool, length(system.bonds))
    n_active_bonds = copy(system.n_neighbors)
    rotation = zeros(9, length(system.bonds))
    rotation[[1, 5, 9], :] .= 1.0
    left_stretch = zeros(9, length(system.bonds))
    left_stretch[[1, 5, 9], :] .= 1.0
    s = BANOSBVerletStorage(position, displacement, velocity, velocity_half, acceleration,
                            b_int, b_ext, damage, bond_active, n_active_bonds, rotation,
                            left_stretch)
    return s
end

@storage BANOSBMaterial VelocityVerlet BANOSBVerletStorage

@loc_to_halo_fields BANOSBVerletStorage :position :velocity
@halo_to_loc_fields BANOSBVerletStorage :b_int

struct BANOSBRelaxationStorage <: AbstractStorage
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
end

function BANOSBRelaxationStorage(::BANOSBMaterial, ::DynamicRelaxation,
                                 system::BondAssociatedSystem, ch)
    n_loc_points = length(ch.loc_points)
    position = copy(system.position)
    displacement = zeros(3, n_loc_points)
    velocity = zeros(3, n_loc_points)
    velocity_half = zeros(3, n_loc_points)
    velocity_half_old = zeros(3, n_loc_points)
    b_int = zeros(3, length(ch.point_ids))
    b_int_old = zeros(3, n_loc_points)
    b_ext = zeros(3, n_loc_points)
    density_matrix = zeros(3, n_loc_points)
    damage = zeros(n_loc_points)
    bond_active = ones(Bool, length(system.bonds))
    n_active_bonds = copy(system.n_neighbors)
    s = BANOSBRelaxationStorage(position, displacement, velocity, velocity_half,
                              velocity_half_old, b_int, b_int_old, b_ext, density_matrix,
                              damage, bond_active, n_active_bonds)
    return s
end

@storage BANOSBMaterial DynamicRelaxation BANOSBRelaxationStorage

@loc_to_halo_fields BANOSBRelaxationStorage :position
@halo_to_loc_fields BANOSBRelaxationStorage :b_int

const BANOSBStorage = Union{BANOSBVerletStorage,BANOSBRelaxationStorage}

function calc_force_density!(chunk::AbstractBodyChunk{S,M},
                             Δt::Float64) where {S<:BondAssociatedSystem,M<:BANOSBMaterial}
    (; system, mat, paramsetup, storage) = chunk
    storage.b_int .= 0
    storage.n_active_bonds .= 0
    for point_id in each_point_idx(chunk)
        force_density_point!(storage, system, mat, paramsetup, Δt, point_id)
    end
    return nothing
end

function force_density_point!(storage::BANOSBStorage, system::BondAssociatedSystem,
                              mat::BANOSBMaterial, paramhandler::AbstractParameterHandler,
                              Δt::Float64, i::Int)
    params = get_params(paramhandler, i)
    force_density_point!(storage, system, mat, params, Δt, i)
    return nothing
end

function force_density_point!(storage::BANOSBStorage, system::BondAssociatedSystem,
                              mat::BANOSBMaterial, params::BANOSBPointParameters,
                              Δt::Float64, i::Int)
    for bond_idx in each_bond_idx(system, i)
        force_density_bond!(storage, system, mat, params, Δt, i, bond_idx)
    end
    return nothing
end

function force_density_bond!(storage::BANOSBStorage, system::BondAssociatedSystem,
                             mat::BANOSBMaterial, params::BANOSBPointParameters,
                             Δt::Float64, i::Int, bond_idx::Int)
    F, Ḟ, Kinv = calc_deformation_gradient(storage, system, mat, params, i, bond_idx)
    if containsnan(F) || storage.damage[i] > mat.maxdmg
        storage.bond_active[bond_idx] = false
        return nothing
    end
    σ = calc_cauchy_stress(mat, params, F)
    calc_rod_and_rotation!(storage, F, Ḟ, Δt, bond_idx)
    R = get_rotation(storage, bond_idx)
    T = R * σ * R'
    # T = σ
    P = det(F) * T * inv(F)'
    PKinv = P * Kinv

    bond = system.bonds[bond_idx]
    j, L = bond.neighbor, bond.length
    ΔXij = get_coordinates_diff(system, i, j)
    Δxij = get_coordinates_diff(storage, i, j)
    l = norm(Δxij)
    ε = (l - L) / L
    stretch_based_failure!(storage, system, bond, params, ε, i, bond_idx)

    ωij = influence_function(mat, params, L) * storage.bond_active[bond_idx]
    tij = ωij * PKinv * ΔXij
    # temp_i = volume_fraction_factor(system, i, bond_idx) * system.volume[i]
    temp_i = system.volume[i]
    # temp_j = volume_fraction_factor(system, j, bond_idx) * system.volume[j]
    temp_j = system.volume[j]
    update_add_b_int!(storage, i, tij .* temp_j)
    update_add_b_int!(storage, j, -tij .* temp_i)
    return nothing
end


@inline function influence_function(::BANOSBMaterial, params::BANOSBPointParameters,
                                    L::Float64)
    return params.δ / L
end

function calc_deformation_gradient(storage::BANOSBStorage, system::BondAssociatedSystem,
                                   mat::BANOSBMaterial, params::BANOSBPointParameters,
                                   i::Int, bond_idx::Int)
    (; bonds, volume) = system
    (; bond_active) = storage
    K = zeros(SMatrix{3,3})
    _F = zeros(SMatrix{3,3})
    _Fdot = zeros(SMatrix{3,3})
    for bond_id in each_intersecting_bond_idx(system, i, bond_idx)
        bond = bonds[bond_id]
        j, L = bond.neighbor, bond.length
        ΔXij = get_coordinates_diff(system, i, j)
        Δxij = get_coordinates_diff(storage, i, j)
        Δvij = get_diff(storage.velocity, i, j)
        ωij = influence_function(mat, params, L) * bond_active[bond_id]
        temp = ωij * volume[j]
        K += temp * ΔXij * ΔXij'
        _F += temp * Δxij * ΔXij'
        _Fdot += temp * Δvij * ΔXij'
    end
    Kinv = inv(K)
    F = _F * Kinv
    Fdot = _Fdot * Kinv
    return F, Fdot, Kinv
end

# function calc_first_piola_stress(F::SMatrix{3,3}, ::BANOSBMaterial,
#                                  params::BANOSBPointParameters)
#     J = det(F)
#     J < eps() && return zero(SMatrix{3,3})
#     C = F' * F
#     Cinv = inv(C)
#     S = params.G .* (I - 1 / 3 .* tr(C) .* Cinv) .* J^(-2 / 3) .+
#         params.K / 4 .* (J^2 - J^(-2)) .* Cinv
#     P = F * S
#     return P
# end

function calc_cauchy_stress(mat::BANOSBMaterial, params::BANOSBPointParameters, F::SMatrix{3,3})
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
    J = det(F)
    J < eps() && return zero(SMatrix{3,3})
    C = F' * F
    Cinv = inv(C)
    S = params.G .* (I - 1 / 3 .* tr(C) .* Cinv) .* J^(-2 / 3) .+
        params.K / 4 .* (J^2 - J^(-2)) .* Cinv
    P = F * S
    σ = 1/J .* P * F'
    # ----

    # ---- Saint Venant-Kirchhoff, taken from
    #      https://doi.org/10.1016/j.jfluidstructs.2021.103312
    # J = det(F)
    # J < eps() && return zero(SMatrix{3,3})
    # # J > mat.maxjacobi && return zero(SMatrix{3,3})
    # E = 0.5 .* (F' * F - I)
    # S = params.λ * tr(E) * I + 2 * params.μ * E
    # P = F * S
    # σ = 1/J .* P * F'
    # ----

    return σ
end

@inline function get_tensor(T::AbstractMatrix, i::Int)
    tensor = SMatrix{3,3}(T[1,i], T[2,i], T[3,i], T[4,i], T[5,i], T[6,i], T[7,i], T[8,i],
                          T[9,i])
    return tensor
end

@inline function get_left_stretch(storage::BANOSBStorage, i::Int)
    return get_tensor(storage.left_stretch, i)
end

@inline function get_rotation(storage::BANOSBStorage, i::Int)
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

@inline function update_left_stretch!(storage::BANOSBStorage, i::Int, V::SMatrix{3,3})
    update_tensor!(storage.left_stretch, i, V)
    return nothing
end

@inline function update_rotation!(storage::BANOSBStorage, i::Int, R::SMatrix{3,3})
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
    # d = Rₙ₊₁' * D * Rₙ₊₁
    # update_tensor!(storage.rate_of_deformation, i, d)
    return nothing
end


function kill_point!(s::AbstractStorage, system::BondAssociatedSystem, i::Int)
    s.bond_active[each_bond_idx(system, i)] .= false
    s.n_active_bonds[i] = 0
    return nothing
end
