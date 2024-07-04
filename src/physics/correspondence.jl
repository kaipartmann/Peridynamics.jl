"""
    NOSBMaterial(; maxdmg, maxjacobi, corr)

Construct the type `NOSBMaterial` used to specify the material of a [`Body`](@ref).

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

function Base.show(io::IO, mat::NOSBMaterial)
    print(io, typeof(mat))
    print(io, msg_fields_in_brackets(mat))
    return nothing
end

"""
TODO
"""
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
end

function NOSBVerletStorage(::NOSBMaterial, ::VelocityVerlet, system::BondSystem, ch)
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
    s = NOSBVerletStorage(position, displacement, velocity, velocity_half, acceleration,
                          b_int, b_ext, damage, bond_active, n_active_bonds)
    return s
end

@storage NOSBMaterial VelocityVerlet NOSBVerletStorage

@loc_to_halo_fields NOSBVerletStorage :position
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
end

function NOSBRelaxationStorage(::NOSBMaterial, ::DynamicRelaxation, system::BondSystem, ch)
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
    s = NOSBRelaxationStorage(position, displacement, velocity, velocity_half,
                              velocity_half_old, b_int, b_int_old, b_ext, density_matrix,
                              damage, bond_active, n_active_bonds)
    return s
end

@storage NOSBMaterial DynamicRelaxation NOSBRelaxationStorage

@loc_to_halo_fields NOSBRelaxationStorage :position
@halo_to_loc_fields NOSBRelaxationStorage :b_int

const NOSBStorage = Union{NOSBVerletStorage,NOSBRelaxationStorage}

function force_density_point!(storage::NOSBStorage, system::BondSystem, mat::NOSBMaterial,
                              paramhandler::AbstractParameterHandler, i::Int)
    params = get_params(paramhandler, i)
    force_density_point!(storage, system, mat, params, i)
    return nothing
end

function force_density_point!(storage::NOSBStorage, system::BondSystem, mat::NOSBMaterial,
                              params::NOSBPointParameters, i::Int)
    F, Kinv, ω0 = calc_deformation_gradient(storage, system, mat, params, i)
    if storage.damage[i] > mat.maxdmg || containsnan(F)
        kill_point!(storage, system, i)
        return nothing
    end
    P = calc_first_piola_stress(F, mat, params)
    if iszero(P) || containsnan(P)
        kill_point!(storage, system, i)
        return nothing
    end
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
        if containsnan(tij)
            tij = zero(SMatrix{3,3})
        end
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
    ω0 = 0.0
    for bond_id in each_bond_idx(system, i)
        bond = system.bonds[bond_id]
        j, L = bond.neighbor, bond.length
        ΔXij = get_coordinates_diff(system, i, j)
        Δxij = get_coordinates_diff(storage, i, j)
        Vj = system.volume[j]
        ωij = influence_function(mat, params, L) * storage.bond_active[bond_id]
        ω0 += ωij
        temp = ωij * Vj
        K += temp * ΔXij * ΔXij'
        _F += temp * Δxij * ΔXij'
    end
    Kinv = inv(K)
    F = _F * Kinv
    return F, Kinv, ω0
end

function calc_first_piola_stress(F::SMatrix{3,3}, mat::NOSBMaterial,
                                 params::NOSBPointParameters)
    J = det(F)
    J < eps() && return zero(SMatrix{3,3})
    J > mat.maxjacobi && return zero(SMatrix{3,3})
    C = F' * F
    Cinv = inv(C)
    S = params.G .* (I - 1 / 3 .* tr(C) .* Cinv) .* J^(-2 / 3) .+
        params.K / 4 .* (J^2 - J^(-2)) .* Cinv
    P = F * S
    return P
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
