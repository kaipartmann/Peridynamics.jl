"""
    GBBMaterial()
    GBBMaterial{Correction}()

A material type used to assign the material of a [`Body`](@ref) with the generalized
bond-based formulation of peridynamics.

# Keywords
- `dmgmodel::AbstractDamageModel`: Damage model defining the fracture behavior.
    (default: `CriticalStretch()`)

Possible correction methods are:
- [`NoCorrection`](@ref): No correction is applied. (default)
- [`EnergySurfaceCorrection`](@ref): The energy based surface correction method of
    Le and Bobaru (2018) is applied.

# Examples

```julia-repl
julia> mat = GBBMaterial()
GBBMaterial{NoCorrection}()

julia> mat = GBBMaterial{EnergySurfaceCorrection}()
GBBMaterial{EnergySurfaceCorrection}()
```
---

```julia
GBBMaterial{Correction}
```

Material type for the dual-horizon bond-based peridynamics formulation.

# Type Parameters
- `Correction`: A correction algorithm type. See the constructor docs for more informations.
- `DM`: A damage model type.

# Allowed material parameters
When using [`material!`](@ref) on a [`Body`](@ref) with `GBBMaterial`, then the following
parameters are allowed:
Material parameters:
- `horizon::Float64`: Radius of point interactions.
- `rho::Float64`: Density.
Elastic parameters:
- `E::Float64`: Young's modulus.
- `G::Float64`: Shear modulus.
- `K::Float64`: Bulk modulus.
- `lambda::Float64`: 1st Lamé parameter.
- `mu::Float64`: 2nd Lamé parameter.
Fracture parameters:
- `Gc::Float64`: Critical energy release rate.
- `epsilon_c::Float64`: Critical strain.

!!! note "Poisson's ratio and bond-based peridynamics"
    In bond-based peridynamics, the Poisson's ratio is limited to 1/4 for 3D simulations.
    Therefore, only one additional elastic parameter is required.
    Optionally, the specification of a second keyword is allowed, if the parameter
    combination results in `nu = 1/4`.

# Allowed export fields
When specifying the `fields` keyword of [`Job`](@ref) for a [`Body`](@ref) with
`GBBMaterial`, the following fields are allowed:
- `position::Matrix{Float64}`: Position of each point.
- `displacement::Matrix{Float64}`: Displacement of each point.
- `velocity::Matrix{Float64}`: Velocity of each point.
- `velocity_half::Matrix{Float64}`: Velocity parameter for Verlet time solver.
- `acceleration::Matrix{Float64}`: Acceleration of each point.
- `b_int::Matrix{Float64}`: Internal force density of each point.
- `b_ext::Matrix{Float64}`: External force density of each point.
- `damage::Vector{Float64}`: Damage of each point.
- `n_active_bonds::Vector{Int}`: Number of intact bonds of each point.
- `strain_energy_density::Vector{Float64}`: Strain energy density of each point.
"""
struct GBBMaterial{Correction,DM} <: AbstractBondBasedMaterial{Correction}
    dmgmodel::DM
    function GBBMaterial{C}(dmgmodel::DM) where {C,DM}
        new{C,DM}(dmgmodel)
    end
end

function GBBMaterial{C}(; dmgmodel::AbstractDamageModel=CriticalStretch()) where {C}
    return GBBMaterial{C}(dmgmodel)
end
GBBMaterial(; kwargs...) = GBBMaterial{NoCorrection}(; kwargs...)

function StandardPointParameters(mat::GBBMaterial{C,D}, p::Dict{Symbol,Any}) where {C,D}
    (; δ, rho, E, nu, G, K, λ, μ) = get_required_point_parameters_bb(mat, p)
    (; Gc, εc) = get_frac_params(mat.dmgmodel, p, δ, K)
    # for a cubic neighborhood, the weighted volume is the integral
    # 8 δ^4 ∫_0^1 ∫_0^1 ∫_0^1 √(x² + y² + z²) dx dy dz
    # which results in 8 * δ^4 * 0.9605919564548167 (solved numerically)
    bc = 18 * K / (8 * δ^4 * 0.9605919564548167) # bond constant
    return StandardPointParameters(δ, rho, E, nu, G, K, λ, μ, Gc, εc, bc)
end

@params GBBMaterial StandardPointParameters

@storage GBBMaterial struct GBBStorage <: AbstractStorage
    @lthfield position::Matrix{Float64}
    @pointfield displacement::Matrix{Float64}
    @pointfield velocity::Matrix{Float64}
    @pointfield velocity_half::Matrix{Float64}
    @pointfield velocity_half_old::Matrix{Float64}
    @pointfield acceleration::Matrix{Float64}
    @pointfield b_int::Matrix{Float64}
    @pointfield b_int_old::Matrix{Float64}
    @pointfield b_ext::Matrix{Float64}
    @pointfield density_matrix::Matrix{Float64}
    @pointfield damage::Vector{Float64}
    @pointfield n_active_bonds::Vector{Int}
    @pointfield strain_energy_density::Vector{Float64}
    @pointfield weighted_volume::Vector{Float64}
    bond_length::Vector{Float64}
    bond_active::Vector{Bool}
    residual::Vector{Float64}
    jacobian::Matrix{Float64}
    displacement_copy::Matrix{Float64}
    b_int_copy::Matrix{Float64}
    temp_force_a::Vector{Float64}
    temp_force_b::Vector{Float64}
    Δu::Vector{Float64}
    affected_points::Vector{Vector{Int}}
end

function init_field(::GBBMaterial, ::AbstractTimeSolver, system::BondSystem,
                    ::Val{:weighted_volume})
    return zeros(get_n_loc_points(system))
end

function calc_weighted_volume!(storage::GBBStorage, system::BondSystem,
                               ::GBBMaterial, ::AbstractParameterSetup, i)
    (; bonds, correction, volume) = system
    wvol = 0.0
    for bond_id in each_bond_idx(system, i)
        bond = bonds[bond_id]
        j, L = bond.neighbor, bond.length
        ω = surface_correction_factor(correction, bond_id)
        wvol += ω * L * volume[j]
    end
    storage.weighted_volume[i] = wvol
    return wvol
end

function force_density_point!(storage::GBBStorage, system::BondSystem, mat::GBBMaterial,
                              params::StandardPointParameters, t, Δt, i)
    (; position, bond_length, bond_active, b_int) = storage
    (; bonds, correction, volume) = system
    wvol = calc_weighted_volume!(storage, system, mat, params, i)
    iszero(wvol) && return nothing
    bond_constant = 18 * params.K / wvol
    for bond_id in each_bond_idx(system, i)
        bond = bonds[bond_id]
        j, L = bond.neighbor, bond.length
        Δxij = get_vector_diff(position, i, j)
        l = bond_length[bond_id]
        ε = (l - L) / L
        ω = bond_active[bond_id] * surface_correction_factor(correction, bond_id)
        b = ω * bond_constant * ε * volume[j] .* Δxij / l
        update_add_vector!(b_int, i, b)
    end
    return nothing
end

function force_density_point!(storage::GBBStorage, system::BondSystem, mat::GBBMaterial,
                              paramhandler::ParameterHandler, t, Δt, i)
    (; position, bond_length, bond_active, b_int) = storage
    (; bonds, correction, volume) = system
    wvol = calc_weighted_volume!(storage, system, mat, paramhandler, i)
    iszero(wvol) && return nothing
    params_i = get_params(paramhandler, i)
    for bond_id in each_bond_idx(system, i)
        bond = bonds[bond_id]
        j, L = bond.neighbor, bond.length
        Δxij = get_vector_diff(position, i, j)
        l = bond_length[bond_id]
        ε = (l - L) / L
        params_j = get_params(paramhandler, j)
        ω = bond_active[bond_id] * surface_correction_factor(correction, bond_id)
        bond_constant = 9 * (params_i.K + params_j.K) / wvol
        b = ω * bond_constant * ε * volume[j] .* Δxij / l
        update_add_vector!(b_int, i, b)
    end
    return nothing
end

# Do not rely on any custom pre-stored properties here!
function strain_energy_density_point!(storage::AbstractStorage, system::BondSystem,
                                      mat::GBBMaterial, paramsetup::AbstractParameterSetup,
                                      i)
    (; position, bond_active, strain_energy_density) = storage
    (; bonds, correction, volume) = system
    params_i = get_params(paramsetup, i)
    wvol = calc_weighted_volume!(storage, system, mat, paramsetup, i)
    iszero(wvol) && return nothing
    Ψ = 0.0
    for bond_id in each_bond_idx(system, i)
        bond = bonds[bond_id]
        j, L = bond.neighbor, bond.length
        Δxij = get_vector_diff(position, i, j)
        l = norm(Δxij) # do not rely on the stored bond length here!
        ε = (l - L) / L
        params_j = get_params(paramsetup, j)
        ωij = bond_active[bond_id] * surface_correction_factor(correction, bond_id)
        bond_constant = 9 * (params_i.K + params_j.K) / wvol
        Ψ += 0.25 * ωij * bond_constant * ε * ε * L * volume[j]
    end
    strain_energy_density[i] = Ψ
    return nothing
end
