"""
    CCMaterial(; maxdmg, zem)

A material type used to assign the material of a [`Body`](@ref) with the local continuum
consistent (correspondence) formulation of non-ordinary state-based peridynamics.

# Keywords
- `maxdmg::Float64`: Maximum value of damage a point is allowed to obtain. If this value is
    exceeded, all bonds of that point are broken because the deformation gradient would then
    possibly contain `NaN` values.
    (default: `0.95`)
- `zem::AbstractZEMStabilization`: zero-energy mode stabilization. The
    stabilization algorithm of Silling (2017) is used as default.
    (default: `ZEMSilling()`)

!!! note "Stability of fracture simulations"
    This formulation is known to be not suitable for fracture simultations without
    stabilization of the zero-energy modes. Therefore be careful when doing fracture
    simulations and try out different parameters for `maxdmg` and `zem`.

# Examples

```julia-repl
julia> mat = CCMaterial()
CCMaterial(maxdmg=0.95, zem_fac=ZEMSilling())
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
- `zem_fac::Float64`: Correction factor used for zero-energy mode stabilization. See the
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
struct BACCMaterial{CM} <: AbstractBondAssociatedSystemMaterial
    constitutive_model::CM
    maxdmg::Float64
    function BACCMaterial(cm::CM, maxdmg::Real) where {CM}
        return new{CM}(cm, maxdmg)
    end
end

function Base.show(io::IO, @nospecialize(mat::BACCMaterial))
    print(io, typeof(mat))
    print(io, msg_fields_in_brackets(mat, (:maxdmg,)))
    return nothing
end

function BACCMaterial(; model::AbstractConstitutiveModel=NeoHookeNonlinear(),
                      maxdmg::Real=0.85)
    return BACCMaterial(model, maxdmg)
end

struct BACCPointParameters <: AbstractPointParameters
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

function BACCPointParameters(mat::BACCMaterial, p::Dict{Symbol,Any})
    (; δ, rho, E, nu, G, K, λ, μ) = get_required_point_parameters(mat, p)
    (; Gc, εc) = get_frac_params(p, δ, K)
    bc = 18 * K / (π * δ^4) # bond constant
    return BACCPointParameters(δ, rho, E, nu, G, K, λ, μ, Gc, εc, bc)
end

@params BACCMaterial BACCPointParameters

@storage BACCMaterial struct BACCStorage
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
    @pointfield stress::Matrix{Float64}
    @pointfield von_mises_stress::Vector{Float64}
end

function init_field(::BACCMaterial, ::AbstractTimeSolver, system::BondAssociatedSystem,
                    ::Val{:b_int})
    return zeros(3, get_n_points(system))
end

function init_field(::BACCMaterial, ::AbstractTimeSolver, system::BondAssociatedSystem,
                    ::Val{:stress})
    return zeros(9, get_n_loc_points(system))
end

function init_field(::BACCMaterial, ::AbstractTimeSolver, system::BondAssociatedSystem,
                    ::Val{:von_mises_stress})
    return zeros(get_n_loc_points(system))
end

# function calc_force_density!(chunk::AbstractBodyChunk{S,M},
#                              Δt::Float64) where {S<:BondAssociatedSystem,M<:BANOSBMaterial}
#     (; system, mat, paramsetup, storage) = chunk
#     storage.b_int .= 0
#     storage.n_active_bonds .= 0
#     for point_id in each_point_idx(chunk)
#         force_density_point!(storage, system, mat, paramsetup, Δt, point_id)
#     end
#     return nothing
# end

function force_density_point!(storage::BACCStorage, system::BondAssociatedSystem,
                              mat::BACCMaterial, paramhandler::AbstractParameterHandler,
                              t, Δt, i)
    params = get_params(paramhandler, i)
    force_density_point!(storage, system, mat, params, t, Δt, i)
    return nothing
end

function force_density_point!(storage::BACCStorage, system::BondAssociatedSystem,
                              mat::BACCMaterial, params::BACCPointParameters, t, Δt, i)
    for bond_idx in each_bond_idx(system, i)
        force_density_bond!(storage, system, mat, params, t, Δt, i, bond_idx)
    end
    return nothing
end

function force_density_bond!(storage::BACCStorage, system::BondAssociatedSystem,
                             mat::BACCMaterial, params::BACCPointParameters, t, Δt, i,
                             bond_idx)
    defgrad_res = calc_deformation_gradient(storage, system, mat, params, i, bond_idx)
    (; F) = defgrad_res
    if containsnan(F) || storage.damage[i] > mat.maxdmg
        storage.bond_active[bond_idx] = false
        return nothing
    end
    PKinv = calc_first_piola_kirchhoff!(storage, mat, params, defgrad_res, Δt, i, bond_idx)

    bond = system.bonds[bond_idx]
    j, L = bond.neighbor, bond.length
    ΔXij = get_coordinates_diff(system, i, j)
    Δxij = get_coordinates_diff(storage, i, j)
    l = norm(Δxij)
    ε = (l - L) / L
    stretch_based_failure!(storage, system, bond, params, ε, i, bond_idx)

    ωij = influence_function(mat, params, L) * storage.bond_active[bond_idx]
    ϕi = volume_fraction_factor(system, i, bond_idx)
    tij = ϕi * ωij * PKinv * ΔXij
    # tij = ωij * PKinv * ΔXij
    update_add_b_int!(storage, i, tij .* system.volume[j])
    update_add_b_int!(storage, j, -tij .* system.volume[i])
    return nothing
end


@inline function influence_function(::BACCMaterial, params::BACCPointParameters,
                                    L::Float64)
    return params.δ / L
end

function calc_deformation_gradient(storage::BACCStorage, system::BondAssociatedSystem,
                                   mat::BACCMaterial, params::BACCPointParameters, i,
                                   bond_idx)
    (; bonds, volume) = system
    (; bond_active) = storage
    K = zero(SMatrix{3,3,Float64,9})
    _F = zero(SMatrix{3,3,Float64,9})
    for bond_id in each_intersecting_bond_idx(system, i, bond_idx)
        bond = bonds[bond_id]
        j, L = bond.neighbor, bond.length
        ΔXij = get_diff(system.position, i, j)
        Δxij = get_diff(storage.position, i, j)

        temp = influence_function(mat, params, L) * volume[j]
        K += temp * (ΔXij * ΔXij')

        temp = influence_function(mat, params, L) * bond_active[bond_id] * volume[j]
        _F += temp * (Δxij * ΔXij')
    end
    Kinv = inv(K)
    F = _F * Kinv
    return (; F, Kinv)
end

function calc_first_piola_kirchhoff!(storage::BACCStorage, mat::BACCMaterial,
                                     params::BACCPointParameters, defgrad_res, Δt, i,
                                     bond_idx)
    (; F, Kinv) = defgrad_res
    σ = cauchy_stress(mat.constitutive_model, storage, params, F)
    # update_tensor!(storage.stress, i, σ)
    # storage.von_mises_stress[i] = von_mises_stress(σ)
    P = det(F) * σ * inv(F)'
    PKinv = P * Kinv
    return PKinv
end
