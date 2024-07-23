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
end

function BANOSBVerletStorage(::BANOSBMaterial, ::VelocityVerlet,
                             system::BondAssociatedSystem, ch)
    n_loc_points = length(ch.loc_points)
    position = copy(system.position)
    displacement = zeros(3, n_loc_points)
    velocity = zeros(3, n_loc_points)
    velocity_half = zeros(3, n_loc_points)
    acceleration = zeros(3, n_loc_points)
    b_int = zeros(3, length(ch.point_ids))
    b_ext = zeros(3, n_loc_points)
    damage = zeros(n_loc_points)
    bond_active = ones(Bool, length(system.bonds))
    n_active_bonds = copy(system.n_neighbors)
    s = BANOSBVerletStorage(position, displacement, velocity, velocity_half, acceleration,
                          b_int, b_ext, damage, bond_active, n_active_bonds)
    return s
end

@storage BANOSBMaterial VelocityVerlet BANOSBVerletStorage

@loc_to_halo_fields BANOSBVerletStorage :position
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

function force_density_point!(storage::BANOSBStorage, system::BondAssociatedSystem,
                              mat::BANOSBMaterial, paramhandler::AbstractParameterHandler,
                              i::Int)
    params = get_params(paramhandler, i)
    force_density_point!(storage, system, mat, params, i)
    return nothing
end

function force_density_point!(storage::BANOSBStorage, system::BondAssociatedSystem,
                              mat::BANOSBMaterial, params::BANOSBPointParameters, i::Int)
    for bond_idx in each_bond_idx(system, i)
        force_density_bond!(storage, system, mat, params, i, bond_idx)
    end
    return nothing
end

function force_density_bond!(storage::BANOSBStorage, system::BondAssociatedSystem,
                             mat::BANOSBMaterial, params::BANOSBPointParameters, i::Int,
                             bond_idx::Int)
    F, Kinv = calc_deformation_gradient(storage, system, mat, params, i, bond_idx)
    if containsnan(F) || storage.damage[i] > mat.maxdmg
        storage.bond_active[bond_idx] = false
        return nothing
    end
    if containsnan(Kinv)
        # @warn "Kinv contains NaN!"
        storage.bond_active[bond_idx] = false
        return nothing
    end
    P = calc_first_piola_stress(F, mat, params)
    if containsnan(P)
        # @warn "P contains NaN!"
        storage.bond_active[bond_idx] = false
        return nothing
    end
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
    temp_i = volume_fraction_factor(system, i, bond_idx) * system.volume[i]
    temp_j = volume_fraction_factor(system, j, bond_idx) * system.volume[j]
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
    for bond_id in each_intersecting_bond_idx(system, i, bond_idx)
        bond = bonds[bond_id]
        j, L = bond.neighbor, bond.length
        ΔXij = get_coordinates_diff(system, i, j)
        Δxij = get_coordinates_diff(storage, i, j)
        ωij = influence_function(mat, params, L) * bond_active[bond_id]
        temp = ωij * volume[j]
        K += temp * ΔXij * ΔXij'
        _F += temp * Δxij * ΔXij'
    end
    Kinv = inv(K)
    F = _F * Kinv
    return F, Kinv
end

function calc_first_piola_stress(F::SMatrix{3,3}, ::BANOSBMaterial,
                                 params::BANOSBPointParameters)
    J = det(F)
    J < eps() && return zero(SMatrix{3,3})
    C = F' * F
    Cinv = inv(C)
    S = params.G .* (I - 1 / 3 .* tr(C) .* Cinv) .* J^(-2 / 3) .+
        params.K / 4 .* (J^2 - J^(-2)) .* Cinv
    P = F * S
    return P
end

function kill_point!(s::AbstractStorage, system::BondAssociatedSystem, i::Int)
    s.bond_active[each_bond_idx(system, i)] .= false
    s.n_active_bonds[i] = 0
    return nothing
end
