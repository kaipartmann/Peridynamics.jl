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
- `strain_energy_density::Vector{Float64}`: Strain energy density of each point.
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

"""
    DHBBPointParameters

$(extension_api_note())

Point parameters of the dual-horizon bond-based material: the [`BBPointParameters`](@ref)
with half of the bond constant, because every bond is visited from both of its points.

$(block_table(DHBBPointParameters))
"""
@params DHBBMaterial struct DHBBPointParameters
    @inherit BBPointParameters
    @derived bc = 0.5 * 18 * K / (π * δ^4) # half of the normal bond constant
end

"""
    DHBBStorage

$(extension_api_note())

Storage of [`DHBBMaterial`](@ref): the fields of [`BBStorage`](@ref), with `b_int` exchanged
halo to local, because the dual-horizon formulation accumulates force density into the
neighbors of a point.

$(block_table(DHBBStorage))
"""
@storage DHBBMaterial struct DHBBStorage <: AbstractStorage
    @inherit VelocityVerletFields DynamicRelaxationFields NewtonKrylovFields
    @inherit BondLengthCache BondFracFields
    @htl b_int::PointVector
    strain_energy_density::PointScalar
    dmg_state::DamageState
end

function force_density_point!(storage::DHBBStorage, system::BondSystem, ::DHBBMaterial,
                              paramsetup::AbstractParameterSetup, t, Δt, i)
    (; position, bond_active, b_int) = storage
    (; bonds, correction, volume) = system
    params_i = get_params(paramsetup, i)
    for bond_id in each_bond_idx(system, i)
        bond = bonds[bond_id]
        j, L = bond.neighbor, bond.length
        Δxij = get_vector_diff(position, i, j)
        l = current_bond_length(storage, system, i, bond_id)
        ε = (l - L) / L
        params_j = get_params(paramsetup, j)
        ω = bond_active[bond_id] * surface_correction_factor(correction, bond_id)
        b = ω * (params_i.bc + params_j.bc) / 2 * ε .* Δxij / l
        update_add_vector!(b_int, i, b * volume[j])
        update_add_vector!(b_int, j, -b * volume[i])
    end
    return nothing
end

function strain_energy_density_point!(storage::AbstractStorage, system::BondSystem,
                                      ::DHBBMaterial, paramsetup::AbstractParameterSetup, i)
    (; bond_active, strain_energy_density) = storage
    (; bonds, correction, volume) = system
    update_bond_lengths!(storage, system, i)
    params_i = get_params(paramsetup, i)
    Ψ = 0.0
    for bond_id in each_bond_idx(system, i)
        bond = bonds[bond_id]
        j, L = bond.neighbor, bond.length
        ε = bond_stretch(storage, system, i, bond_id)
        params_j = get_params(paramsetup, j)
        ωij = bond_active[bond_id] * surface_correction_factor(correction, bond_id)
        bc = (params_i.bc + params_j.bc) / 2
        Ψ += 0.5 * ωij * bc * ε * ε * L * volume[j] # added factor 2 here!
    end
    strain_energy_density[i] = Ψ
    return nothing
end
