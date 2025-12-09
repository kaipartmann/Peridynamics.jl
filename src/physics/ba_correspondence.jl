"""
    BACMaterial(; kernel, model, dmgmodel, maxdmg)

A material type used to assign the material of a [`Body`](@ref) with the bond-associated
correspondence formulation of Chen and Spencer (2019).

# Keywords
- `kernel::Function`: Kernel function used for weighting the interactions between points.
    (default: `linear_kernel`)
- `model::AbstractConstitutiveModel`: Constitutive model defining the material behavior.
    (default: `SaintVenantKirchhoff()`)
- `dmgmodel::AbstractDamageModel`: Damage model defining the fracture behavior.
    (default: `CriticalStretch()`)
- `maxdmg::Float64`: Maximum value of damage a point is allowed to obtain. If this value is
    exceeded, all bonds of that point are broken because the deformation gradient would then
    possibly contain `NaN` values.
    (default: `0.85`)

# Examples

```julia-repl
julia> mat = BACMaterial()
BACMaterial{SaintVenantKirchhoff, typeof(linear_kernel), CriticalStretch}()
```
---

```julia
BACMaterial{CM,K,DM}
```

Material type for the bond-associated correspondence formulation of Chen and Spencer (2019).

# Type Parameters
- `CM`: A constitutive model type. See the constructor docs for more informations.
- `K`: A kernel function type. See the constructor docs for more informations.
- `DM`: A damage model type.

# Fields
- `kernel::Function`: Kernel function used for weighting the interactions between points.
- `model::AbstractConstitutiveModel`: Constitutive model defining the material behavior.
- `dmgmodel::AbstractDamageModel`: Damage model defining the fracture behavior.
- `maxdmg::Float64`: Maximum value of damage a point is allowed to obtain. See the
    constructor docs for more informations.

# Allowed material parameters
When using [`material!`](@ref) on a [`Body`](@ref) with `BACMaterial`, then the following
parameters are allowed:
- `horizon::Float64`: Radius of point interactions.
- `rho::Float64`: Density.
- `E::Float64`: Young's modulus.
- `nu::Float64`: Poisson's ratio.
- `Gc::Float64`: Critical energy release rate.
- `epsilon_c::Float64`: Critical strain.

# Allowed export fields
When specifying the `fields` keyword of [`Job`](@ref) for a [`Body`](@ref) with
`BACMaterial`, the following fields are allowed:
- `position::Matrix{Float64}`: Position of each point.
- `displacement::Matrix{Float64}`: Displacement of each point.
- `velocity::Matrix{Float64}`: Velocity of each point.
- `velocity_half::Matrix{Float64}`: Velocity parameter for Verlet time solver.
- `acceleration::Matrix{Float64}`: Acceleration of each point.
- `b_int::Matrix{Float64}`: Internal force density of each point.
- `b_ext::Matrix{Float64}`: External force density of each point.
- `damage::Vector{Float64}`: Damage of each point.
- `n_active_bonds::Vector{Int}`: Number of intact bonds of each point.
"""
struct BACMaterial{CM,K,DM} <: AbstractBondAssociatedSystemMaterial
    kernel::K
    constitutive_model::CM
    dmgmodel::DM
    maxdmg::Float64
    function BACMaterial(kernel::K, cm::CM, dmgmodel::DM, maxdmg::Real) where {K,CM,DM}
        return new{CM,K,DM}(kernel, cm, dmgmodel, maxdmg)
    end
end

function BACMaterial(; kernel::Function=linear_kernel,
                     model::AbstractConstitutiveModel=SaintVenantKirchhoff(),
                     dmgmodel::AbstractDamageModel=CriticalStretch(),
                     maxdmg::Real=0.85)
    return BACMaterial(kernel, model, dmgmodel, maxdmg)
end

function Base.show(io::IO, @nospecialize(mat::BACMaterial))
    print(io, typeof(mat))
    print(io, msg_fields_in_brackets(mat, (:maxdmg,)))
    return nothing
end

function log_material_property(::Val{:kernel}, mat::BACMaterial; indentation::Int)
    msg = msg_qty("kernel function", mat.kernel; indentation)
    return msg
end

function log_material_property(::Val{:dmgmodel}, mat::BACMaterial; indentation::Int)
    msg = msg_qty("damage model type", typeof(mat.dmgmodel); indentation)
    return msg
end

function log_material_property(::Val{:constitutive_model}, mat::BACMaterial;
                               indentation::Int)
    msg = msg_qty("constitutive model", mat.constitutive_model; indentation)
    return msg
end

function log_material_property(::Val{:maxdmg}, mat::BACMaterial; indentation::Int)
    msg = msg_qty("maximum damage", mat.maxdmg; indentation)
    return msg
end

"""
    BACPointParameters

$(internal_api_warning())

Type containing the material parameters for a peridynamics model using the bond-associated
correspondence formulation of Chen and Spencer.

# Fields

- `δ::Float64`: Horizon.
- `δb::Float64`: Bond-associated horizon.
- `rho::Float64`: Density.
- `E::Float64`: Young's modulus.
- `nu::Float64`: Poisson's ratio.
- `G::Float64`: Shear modulus.
- `K::Float64`: Bulk modulus.
- `λ::Float64`: 1st Lamé parameter.
- `μ::Float64`: 2nd Lamé parameter.
- `Gc::Float64`: Critical energy release rate.
- `εc::Float64`: Critical strain.
- `bc::Float64`: Bond constant.
"""
struct BACPointParameters <: AbstractPointParameters
    δ::Float64
    δb::Float64
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

function BACPointParameters(mat::BACMaterial, p::Dict{Symbol,Any})
    (; δ, δb, rho, E, nu, G, K, λ, μ) = get_required_point_parameters(mat, p)
    (; Gc, εc) = get_frac_params(mat.dmgmodel, p, δ, K)
    # for a cubic neighborhood, the weighted volume is the integral
    # 8 δ^4 ∫_0^1 ∫_0^1 ∫_0^1 √(x² + y² + z²) dx dy dz
    # which results in 8 * δ^4 * 0.9605919564548167 (solved numerically)
    bc = 18 * K / (8 * δ^4 * 0.9605919564548167) # bond constant
    return BACPointParameters(δ, δb, rho, E, nu, G, K, λ, μ, Gc, εc, bc)
end

@params BACMaterial BACPointParameters

@storage BACMaterial struct BACStorage
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
    @pointfield n_active_bonds::Vector{Int}
    @pointfield stress::Matrix{Float64}
    @pointfield von_mises_stress::Vector{Float64}
    bond_active::Vector{Bool}
    residual::Vector{Float64}
    displacement_copy::Matrix{Float64}
    b_int_copy::Matrix{Float64}
    temp_force::Vector{Float64}
    Δu::Vector{Float64}
    v_temp::Vector{Float64}
    Jv_temp::Vector{Float64}
    precond_diag::Vector{Float64}
end

function init_field(::BACMaterial, ::AbstractTimeSolver, system::BondAssociatedSystem,
                    ::Val{:b_int})
    return zeros(3, get_n_points(system))
end

function init_field(::BACMaterial, ::AbstractTimeSolver, system::BondAssociatedSystem,
                    ::Val{:stress})
    return zeros(9, get_n_loc_points(system))
end

function init_field(::BACMaterial, ::AbstractTimeSolver, system::BondAssociatedSystem,
                    ::Val{:von_mises_stress})
    return zeros(get_n_loc_points(system))
end

function force_density_point!(storage::BACStorage, system::BondAssociatedSystem,
                              mat::BACMaterial, params::BACPointParameters, t, Δt, i)
    for bond_idx in each_bond_idx(system, i)
        force_density_bond!(storage, system, mat, params, t, Δt, i, bond_idx)
    end
    return nothing
end

function force_density_bond!(storage::BACStorage, system::BondAssociatedSystem,
                             mat::BACMaterial, params::BACPointParameters, t, Δt, i,
                             bond_idx)
    defgrad_res = calc_deformation_gradient(storage, system, mat, params, i, bond_idx)
    (; F) = defgrad_res
    if containsnan(F) || storage.damage[i] > mat.maxdmg
        storage.bond_active[bond_idx] = false
        return nothing
    end
    PKinv = calc_first_piola_kirchhoff!(storage, mat, params, defgrad_res, Δt, i, bond_idx)

    bond = system.bonds[bond_idx]
    j = bond.neighbor
    ΔXij = get_vector_diff(system.position, i, j)

    ωij = kernel(system, bond_idx) * storage.bond_active[bond_idx]
    ϕi = volume_fraction_factor(system, i, bond_idx)
    tij = ϕi * ωij * PKinv * ΔXij
    update_add_vector!(storage.b_int, i, tij .* system.volume[j])
    update_add_vector!(storage.b_int, j, -tij .* system.volume[i])
    return nothing
end

function calc_deformation_gradient(storage::BACStorage, system::BondAssociatedSystem,
                                   mat::BACMaterial, params::BACPointParameters, i,
                                   bond_idx)
    (; bonds, volume) = system
    (; bond_active) = storage
    K = zero(SMatrix{3,3,Float64,9})
    _F = zero(SMatrix{3,3,Float64,9})
    for bond_id in each_intersecting_bond_idx(system, i, bond_idx)
        bond = bonds[bond_id]
        j = bond.neighbor
        ΔXij = get_vector_diff(system.position, i, j)
        Δxij = get_vector_diff(storage.position, i, j)
        ωijV = kernel(system, bond_id) * volume[j]
        ωijωDV = ωijV * bond_active[bond_id]
        K += ωijV * (ΔXij * ΔXij')
        _F += ωijωDV * (Δxij * ΔXij')
    end
    Kinv = inv(K)
    F = _F * Kinv
    return (; F, Kinv)
end

function calc_first_piola_kirchhoff!(storage::BACStorage, mat::BACMaterial,
                                     params::BACPointParameters, defgrad_res, Δt, i,
                                     bond_idx)
    (; F, Kinv) = defgrad_res
    P = first_piola_kirchhoff(mat.constitutive_model, storage, params, F)
    PKinv = P * Kinv
    return PKinv
end
