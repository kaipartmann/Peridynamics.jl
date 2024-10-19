"""
    CCMaterial(; maxdmg, maxjacobi, corr)

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
julia> mat = CCMaterial()
CCMaterial(maxdmg=0.95, maxjacobi=1.03, corr=100.0)
```

---

```julia
CCMaterial
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
When using [`material!`](@ref) on a [`Body`](@ref) with `CCMaterial`, then the following
parameters are allowed:
- `horizon::Float64`: Radius of point interactions
- `rho::Float64`: Density
- `E::Float64`: Young's modulus
- `nu::Float64`: Poisson's ratio
- `Gc::Float64`: Critical energy release rate
- `epsilon_c::Float64`: Critical strain

# Allowed export fields
When specifying the `fields` keyword of [`Job`](@ref) for a [`Body`](@ref) with
`CCMaterial`, the following fields are allowed:
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
Base.@kwdef struct CCMaterial <: AbstractBondSystemMaterial{NoCorrection}
    maxdmg::Float64 = 0.95
    maxjacobi::Float64 = 1.03
    corr::Float64 = 100.0
end

function Base.show(io::IO, @nospecialize(mat::CCMaterial))
    print(io, typeof(mat))
    print(io, msg_fields_in_brackets(mat))
    return nothing
end

struct CCPointParameters <: AbstractPointParameters
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

function CCPointParameters(mat::CCMaterial, p::Dict{Symbol,Any})
    (; δ, rho, E, nu, G, K, λ, μ) = get_required_point_parameters(mat, p)
    (; Gc, εc) = get_frac_params(p, δ, K)
    bc = 18 * K / (π * δ^4) # bond constant
    return CCPointParameters(δ, rho, E, nu, G, K, λ, μ, Gc, εc, bc)
end

@params CCMaterial CCPointParameters

@storage CCMaterial struct CCStorage <: AbstractStorage
    @lthfield position::Matrix{Float64}
    @pointfield displacement::Matrix{Float64}
    @pointfield velocity::Matrix{Float64}
    @pointfield velocity_half::Matrix{Float64}
    @pointfield velocity_half_old::Matrix{Float64}
    @pointfield acceleration::Matrix{Float64}
    @htlfield b_int::Matrix{Float64}
    @pointfield b_int_old::Matrix{Float64}
    @pointfield b_ext::Matrix{Float64}
    @pointfield density_matrix::Matrix{Float64}
    @pointfield damage::Vector{Float64}
    bond_active::Vector{Bool}
    @pointfield n_active_bonds::Vector{Int}
end

function init_field(::CCMaterial, ::AbstractTimeSolver, system::BondSystem, ::Val{:b_int})
    return zeros(3, get_n_points(system))
end

function force_density_point!(storage::CCStorage, system::BondSystem, mat::CCMaterial,
                              paramhandler::AbstractParameterHandler, t, Δt, i)
    params = get_params(paramhandler, i)
    force_density_point!(storage, system, mat, params, t, Δt, i)
    return nothing
end

function force_density_point!(storage::CCStorage, system::BondSystem, mat::CCMaterial,
                              params::CCPointParameters, t, Δt, i)
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

@inline function influence_function(::CCMaterial, params::CCPointParameters, L)
    return params.δ / L
end

function calc_deformation_gradient(storage::CCStorage, system::BondSystem, mat::CCMaterial,
                                   params::CCPointParameters, i)
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

function calc_first_piola_stress(F::SMatrix{3,3}, mat::CCMaterial,
                                 params::CCPointParameters)
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

function kill_point!(s::AbstractStorage, bd::BondSystem, i)
    s.bond_active[each_bond_idx(bd, i)] .= false
    s.n_active_bonds[i] = 0
    return nothing
end
