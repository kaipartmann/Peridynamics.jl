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
- `maxdmg::Float64`: Maximum value of damage a point or a bond-associated family is allowed
    to obtain. If this value is exceeded, all bonds of that point are broken or the family
    stops carrying stress, because the deformation gradient would then possibly contain
    `NaN` values.
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
- `maxdmg::Float64`: Maximum value of damage a point or a bond-associated family is allowed
    to obtain. See the constructor docs for more informations.

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

@inline get_constitutive_model(mat::BACMaterial) = mat.constitutive_model

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
    bc = 18 * K / (π * δ^4) # bond constant
    return BACPointParameters(δ, δb, rho, E, nu, G, K, λ, μ, Gc, εc, bc)
end

@params BACMaterial BACPointParameters

@storage BACMaterial struct BACStorage
    @inherit VelocityVerletFields DynamicRelaxationFields NewtonKrylovFields
    @inherit BondFracFields
    @htl b_int::PointVector
    stress::PointTensor
    von_mises_stress::PointScalar
    cm_state::ConstitutiveState
    # scratch space for a single point, see `init_field` below
    bond_stress::Matrix{Float64}
end

# scratch space for a single point, filled and used up inside `force_density_point!`
function init_field(::BACMaterial, ::AbstractTimeSolver, system::BondAssociatedSystem,
                    ::Val{:bond_stress})
    return zeros(9, max(1, maximum(system.n_neighbors; init=0)))
end

# First collect the stress of all bond-associated families in `bond_stress`, then turn it
# into forces. Same sum as writing the force of every family directly to all of its bonds,
# but with a factor of the family size fewer scattered writes to `b_int`.
function force_density_point!(storage::BACStorage, system::BondAssociatedSystem,
                              mat::BACMaterial, params::BACPointParameters, t, Δt, i)
    bond_ids_of_i = each_bond_idx(system, i)
    for k in eachindex(bond_ids_of_i)
        zero_tensor!(storage.bond_stress, k)
    end
    for bond_idx in bond_ids_of_i
        collect_bond_stress!(storage, system, mat, params, t, Δt, i, bond_idx)
    end
    # `intersection_bond_ids` indexes the bonds of a point from one, so this is the shift
    # between that numbering and the bond indices of the chunk
    offset = first(bond_ids_of_i) - 1
    for k in eachindex(bond_ids_of_i)
        bond_idx = offset + k
        storage.bond_active[bond_idx] || continue
        j = system.bonds[bond_idx].neighbor
        ΔXij = get_vector_diff(system.position, i, j)
        tij = kernel(system, bond_idx) * (get_tensor(storage.bond_stress, k) * ΔXij)
        update_add_vector!(storage.b_int, i, tij .* system.volume[j])
        update_add_vector!(storage.b_int, j, -tij .* system.volume[i])
    end
    return nothing
end

# The stress of a bond-associated family belongs to all bonds of that family and not only to
# the bond it was calculated for, because `F` is built from the deformation of the whole
# family. Summed over the bonds of a point this recovers `Σ_b w_b P_b`, which is what makes
# the formulation conserve energy and momentum.
function collect_bond_stress!(storage::BACStorage, system::BondAssociatedSystem,
                              mat::BACMaterial, params::BACPointParameters, t, Δt, i,
                              bond_idx)
    if storage.damage[i] > mat.maxdmg
        storage.bond_active[bond_idx] = false
        return nothing
    end
    # a broken bond has no share of the energy, so its family is never needed
    storage.bond_active[bond_idx] || return nothing
    defgrad_res = calc_deformation_gradient(storage, system, mat, params, i, bond_idx)
    (; F, too_damaged) = defgrad_res
    # without a usable deformation gradient there is no stress, but breaking the bond is
    # still up to the damage model
    (too_damaged || containsnan(F)) && return nothing
    PKinv = calc_first_piola_kirchhoff!(storage, mat, params, defgrad_res, Δt, i, bond_idx)

    wPKinv = volume_fraction_factor(system, i, bond_idx) * PKinv
    offset = first(each_bond_idx(system, i)) - 1
    for k in system.intersection_bond_ids[bond_idx]
        storage.bond_active[offset + k] || continue
        update_add_tensor!(storage.bond_stress, k, wPKinv)
    end
    return nothing
end

# Lower limit for the shape quality `det(K) / (tr(K)/3)^3` of a bond-associated family,
# which is one for bonds spread evenly over all directions and goes to zero as they collapse
# onto a plane or a line. Without this limit a nearly singular `K` produces force densities
# of `1e19` that break everything around them. Intact families are at `0.289` and above,
# degenerate ones many orders of magnitude below, so anything from here down to `1e-14`
# gives the same results. `containsnan(F)` is no substitute, it only fires if `K` is exactly
# singular.
const BA_MIN_SHAPE_QUALITY = 1e-3

# Broken bonds are left out of both sums, so that `F` stays the identity for an undeformed
# body. A bond-associated family is much smaller than the family of a point and degenerates
# long before the point damage reaches `maxdmg`, so `too_damaged` applies the same criterion
# per family: too much lost volume or a collapsed shape and the family carries no stress.
function calc_deformation_gradient(storage::BACStorage, system::BondAssociatedSystem,
                                   mat::BACMaterial, params::BACPointParameters, i,
                                   bond_idx)
    (; bonds, volume, ba_hood_volume) = system
    (; bond_active) = storage
    K = zero(SMatrix{3,3,Float64,9})
    _F = zero(SMatrix{3,3,Float64,9})
    intact_volume = 0.0
    for bond_id in each_intersecting_bond_idx(system, i, bond_idx)
        bond_active[bond_id] || continue
        bond = bonds[bond_id]
        j = bond.neighbor
        ΔXij = get_vector_diff(system.position, i, j)
        Δxij = get_vector_diff(storage.position, i, j)
        ωijV = kernel(system, bond_id) * volume[j]
        K += ωijV * (ΔXij * ΔXij')
        _F += ωijV * (Δxij * ΔXij')
        intact_volume += volume[j]
    end
    family_damage = 1 - intact_volume / ba_hood_volume[bond_idx]
    quality = det(K) / (tr(K) / 3)^3
    # `!(quality > …)` and not `quality ≤ …`, so that a `NaN` counts as degenerate
    too_damaged = family_damage > mat.maxdmg || !(quality > BA_MIN_SHAPE_QUALITY)
    Kinv = inv(K)
    F = _F * Kinv
    return (; F, Kinv, too_damaged)
end

function calc_first_piola_kirchhoff!(storage::BACStorage, mat::BACMaterial,
                                     params::BACPointParameters, defgrad_res, Δt, i,
                                     bond_idx)
    (; F, Kinv) = defgrad_res
    P = first_piola_kirchhoff(mat.constitutive_model, storage, params, F, bond_idx, Δt)
    PKinv = P * Kinv
    return PKinv
end
