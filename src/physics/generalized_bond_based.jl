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

# the point parameters of the bond-based material; `bc` is its bond constant
@params GBBMaterial BBPointParameters

"""
    GBBStorage

$(extension_api_note())

Storage of [`GBBMaterial`](@ref): the fields of [`BBStorage`](@ref) and the weighted volume
of every point.

$(block_table(GBBStorage))
"""
@storage GBBMaterial struct GBBStorage <: AbstractStorage
    @inherit VelocityVerletFields DynamicRelaxationFields NewtonKrylovFields
    @inherit BondLengthCache
    strain_energy_density::PointScalar
    weighted_volume::PointScalar
    dmg_state::DamageState
end

function calc_weighted_volume!(storage::GBBStorage, system::BondSystem,
                               ::GBBMaterial, ::AbstractParameterSetup, i)
    (; volume) = system
    wvol = 0.0
    for bond_id in each_bond_idx(system, i)
        j = get_neighbor(system, bond_id)
        L = reference_bond_length(system, bond_id)
        ω = surface_correction_factor(system, bond_id)
        wvol += ω * L * volume[j]
    end
    storage.weighted_volume[i] = wvol
    return wvol
end

function force_density_point!(storage::GBBStorage, system::BondSystem, mat::GBBMaterial,
                              paramsetup::AbstractParameterSetup, t, Δt, i)
    (; position, b_int) = storage
    (; volume) = system
    wvol = calc_weighted_volume!(storage, system, mat, paramsetup, i)
    iszero(wvol) && return nothing
    params_i = get_params(paramsetup, i)
    for bond_id in each_bond_idx(system, i)
        j = get_neighbor(system, bond_id)
        L = reference_bond_length(system, bond_id)
        Δxij = get_vector_diff(position, i, j, dims(system))
        l = current_bond_length(storage, system, i, bond_id)
        ε = (l - L) / L
        params_j = get_params(paramsetup, j)
        ω = bond_is_active(storage, system, bond_id) *
            surface_correction_factor(system, bond_id)
        bond_constant = 9 * (params_i.K + params_j.K) / wvol
        b = ω * bond_constant * ε * volume[j] .* Δxij / l
        update_add_vector!(b_int, i, b, dims(system))
    end
    return nothing
end

function strain_energy_density_point!(storage::AbstractStorage, system::BondSystem,
                                      mat::GBBMaterial, paramsetup::AbstractParameterSetup,
                                      i)
    (; strain_energy_density) = storage
    (; volume) = system
    update_bond_lengths!(storage, system, i)
    params_i = get_params(paramsetup, i)
    wvol = calc_weighted_volume!(storage, system, mat, paramsetup, i)
    iszero(wvol) && return nothing
    Ψ = 0.0
    for bond_id in each_bond_idx(system, i)
        j = get_neighbor(system, bond_id)
        L = reference_bond_length(system, bond_id)
        ε = bond_stretch(storage, system, i, bond_id)
        params_j = get_params(paramsetup, j)
        ωij = bond_is_active(storage, system, bond_id) *
              surface_correction_factor(system, bond_id)
        bond_constant = 9 * (params_i.K + params_j.K) / wvol
        Ψ += 0.25 * ωij * bond_constant * ε * ε * L * volume[j]
    end
    strain_energy_density[i] = Ψ
    return nothing
end
