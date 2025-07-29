"""
    DHBBMaterial()
    DHBBMaterial{Correction}()

A material type used to assign the material of a [`Body`](@ref) with the dual-horizon
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
julia> mat = DHBBMaterial()
DHBBMaterial{NoCorrection}()

julia> mat = DHBBMaterial{EnergySurfaceCorrection}()
DHBBMaterial{EnergySurfaceCorrection}()
```
---

```julia
DHBBMaterial{Correction}
```

Material type for the dual-horizon bond-based peridynamics formulation.

# Type Parameters
- `Correction`: A correction algorithm type. See the constructor docs for more informations.
- `DM`: A damage model type.

# Allowed material parameters
When using [`material!`](@ref) on a [`Body`](@ref) with `DHBBMaterial`, then the following
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
`DHBBMaterial`, the following fields are allowed:
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
struct DHBBMaterial{Correction,DM} <: AbstractBondBasedMaterial{Correction}
    dmgmodel::DM
    function DHBBMaterial{C}(dmgmodel::DM) where {C,DM}
        new{C,DM}(dmgmodel)
    end
end

function DHBBMaterial{C}(; dmgmodel::AbstractDamageModel=CriticalStretch()) where {C}
    return DHBBMaterial{C}(dmgmodel)
end
DHBBMaterial(; kwargs...) = DHBBMaterial{NoCorrection}(; kwargs...)

function StandardPointParameters(mat::DHBBMaterial{C,D}, p::Dict{Symbol,Any}) where {C,D}
    (; δ, rho, E, nu, G, K, λ, μ) = get_required_point_parameters_bb(mat, p)
    (; Gc, εc) = get_frac_params(mat.dmgmodel, p, δ, K)
    # for a cubic neighborhood, the weighted volume is the integral
    # 8 δ^4 ∫_0^1 ∫_0^1 ∫_0^1 √(x² + y² + z²) dx dy dz
    # which results in 8 * δ^4 * 0.9605919564548167 (solved numerically)
    bc = 0.5 * 18 * K / (8 * δ^4 * 0.9605919564548167) # half of the normal bond constant
    return StandardPointParameters(δ, rho, E, nu, G, K, λ, μ, Gc, εc, bc)
end

@params DHBBMaterial StandardPointParameters

@storage DHBBMaterial struct DHBBStorage <: AbstractStorage
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
    bond_stretch::Vector{Float64}
    bond_active::Vector{Bool}
    @pointfield n_active_bonds::Vector{Int}
end

function init_field(::DHBBMaterial, ::AbstractTimeSolver, system::BondSystem, ::Val{:b_int})
    return zeros(3, get_n_points(system))
end

function init_field(::DHBBMaterial, ::AbstractTimeSolver, system::BondSystem,
                    ::Val{:bond_stretch})
    return zeros(get_n_bonds(system))
end

function force_density_point!(storage::DHBBStorage, system::BondSystem, ::DHBBMaterial,
                              params::StandardPointParameters, t, Δt, i)
    (; position, bond_stretch, bond_active, b_int) = storage
    (; bonds, correction, volume) = system
    for bond_id in each_bond_idx(system, i)
        bond = bonds[bond_id]
        j = bond.neighbor
        Δxij = get_vector_diff(position, i, j)
        ε = bond_stretch[bond_id]
        ω = bond_active[bond_id] * surface_correction_factor(correction, bond_id)
        b = ω * params.bc * ε .* Δxij
        update_add_vector!(b_int, i, b * volume[j])
        update_add_vector!(b_int, j, -b * volume[i])
    end
    return nothing
end

function force_density_point!(storage::DHBBStorage, system::BondSystem, ::DHBBMaterial,
                              paramhandler::ParameterHandler, t, Δt, i)
    (; position, bond_stretch, bond_active, b_int) = storage
    (; bonds, correction, volume) = system
    params_i = get_params(paramhandler, i)
    for bond_id in each_bond_idx(system, i)
        bond = bonds[bond_id]
        j = bond.neighbor
        Δxij = get_vector_diff(position, i, j)
        ε = bond_stretch[bond_id]
        params_j = get_params(paramhandler, j)
        ω = bond_active[bond_id] * surface_correction_factor(correction, bond_id)
        b = ω * (params_i.bc + params_j.bc) / 2 * ε .* Δxij
        update_add_vector!(b_int, i, b * volume[j])
        update_add_vector!(b_int, j, -b * volume[i])
    end
    return nothing
end
