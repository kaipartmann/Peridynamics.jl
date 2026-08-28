"""
    BBMaterial()
    BBMaterial{Correction}()

A material type used to assign the material of a [`Body`](@ref) with the standard bond-based
formulation of peridynamics.

# Keywords
- `dmgmodel::AbstractDamageModel`: Damage model defining the fracture behavior.
    (default: `CriticalStretch()`)

Possible correction methods are:
- [`NoCorrection`](@ref): No correction is applied. (default)
- [`EnergySurfaceCorrection`](@ref): The energy based surface correction method of
    Le and Bobaru (2018) is applied.

# Examples

```julia-repl
julia> mat = BBMaterial()
BBMaterial{NoCorrection}()

julia> mat = BBMaterial{EnergySurfaceCorrection}()
BBMaterial{EnergySurfaceCorrection}()
```
---

```julia
BBMaterial{Correction}
```

Material type for the bond-based peridynamics formulation.

# Type Parameters
- `Correction`: A correction algorithm type. See the constructor docs for more informations.
- `DM`: A damage model type.

# Allowed material parameters
When using [`material!`](@ref) on a [`Body`](@ref) with `BBMaterial`, then the following
parameters are allowed:
Material parameters:
- `horizon::Float64`: Radius of point interactions.
- `rho::Float64`: Density.
Elastic parameters
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
`BBMaterial`, the following fields are allowed:
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
struct BBMaterial{Correction,DM} <: AbstractBondBasedMaterial{Correction}
    dmgmodel::DM
    function BBMaterial{C}(dmgmodel::DM) where {C,DM}
        new{C,DM}(dmgmodel)
    end
end

function BBMaterial{C}(; dmgmodel::AbstractDamageModel=CriticalStretch()) where {C}
    return BBMaterial{C}(dmgmodel)
end
BBMaterial(; kwargs...) = BBMaterial{NoCorrection}(; kwargs...)

"""
    BBElasticParameters

$(extension_api_note())

Parameter block of the six elastic parameters of a bond-based material. Bond-based
peridynamics fixes the Poisson's ratio at `nu = 0.25`, so a single elastic keyword is enough
and any combination that results in another value is rejected. See [`@params_fields`](@ref).

$(block_table(BBElasticParameters))
"""
@params_fields BBElasticParameters begin
    @derived (; E, nu, G, K, λ, μ) = get_elastic_params_bb(; E, nu, G, K, lambda, mu)
    @log "Young's modulus" E
    @log "Poisson's ratio" nu
    @log "shear modulus" G
    @log "bulk modulus" K
end

function get_elastic_params_bb(; E=nothing, nu=nothing, G=nothing, K=nothing,
                               lambda=nothing, mu=nothing)
    given = get_given_elastic_params(; E, nu, G, K, lambda, mu)
    if isfinite(given.nu) && !isapprox(given.nu, 0.25)
        throw(ArgumentError(bb_poissons_ratio_msg()))
    elseif !isfinite(given.nu) && length(findall(isfinite, given)) == 1
        given = merge(given, (; nu=0.25))
    end
    elastic_params = resolve_elastic_params(given)
    if !isapprox(elastic_params.nu, 0.25)
        msg = bb_poissons_ratio_msg()
        msg *= "The submitted parameter combination results in an illegal value for nu!\n"
        msg *= "Please define either only one or two fitting elastic parameters!\n"
        throw(ArgumentError(msg))
    end
    return elastic_params
end

function bb_poissons_ratio_msg()
    msg = "Bond-based peridynamics has a limitation on the Poisson's ratio!\n"
    msg *= "With BBMaterial, no other values than nu=0.25 are allowed!\n"
    return msg
end

"""
    BBPointParameters

$(extension_api_note())

Point parameters of the bond-based material: the discretization parameters, the elastic
parameters with the Poisson's ratio of bond-based peridynamics, the bond constant `bc` and
the parameters of the damage model. [`GBBMaterial`](@ref) uses them as they are,
[`DHBBMaterial`](@ref) inherits them and halves the bond constant.

$(block_table(BBPointParameters))
"""
@params BBMaterial struct BBPointParameters
    @inherit DiscretizationParameters BBElasticParameters
    @derived bc = 18 * K / (π * δ^4)
    dmg_params::DamageParameters
end

"""
    BBStorage

$(extension_api_note())

Storage of [`BBMaterial`](@ref): the fields of the three time solvers and of the fracture
bookkeeping, the strain energy density of every point, the current length of every bond
and the state of the damage model.

$(block_table(BBStorage))
"""
@storage BBMaterial struct BBStorage <: AbstractStorage
    @inherit VelocityVerletFields DynamicRelaxationFields NewtonKrylovFields
    @inherit BondFracFields
    strain_energy_density::PointScalar
    bond_length::BondScalar
    dmg_state::DamageState
end

# Customized calc_failure to save the bond stretch ε for force density calculation
function calc_failure!(storage::AbstractStorage, system::BondSystem,
                       ::AbstractBondBasedMaterial, ::CriticalStretch,
                       paramsetup::AbstractParameterSetup, t, Δt, i)
    (; εc) = get_params(paramsetup, i)
    (; position, n_active_bonds, bond_active, bond_length) = storage
    (; bonds) = system
    for bond_id in each_bond_idx(system, i)
        bond = bonds[bond_id]
        j, L = bond.neighbor, bond.length
        Δxij = get_vector_diff(position, i, j)
        l = norm(Δxij)
        ε = (l - L) / L
        bond_length[bond_id] = l # store current bond length
        if ε > εc && bond.fail_permit
            bond_active[bond_id] = false
        end
        n_active_bonds[i] += bond_active[bond_id]
    end
    return nothing
end

function force_density_point!(storage::BBStorage, system::BondSystem, ::BBMaterial,
                              params::BBPointParameters, t, Δt, i)
    (; position, bond_length, bond_active, b_int) = storage
    (; bonds, correction, volume) = system
    for bond_id in each_bond_idx(system, i)
        bond = bonds[bond_id]
        j, L = bond.neighbor, bond.length
        Δxij = get_vector_diff(position, i, j)
        l = bond_length[bond_id]
        ε = (l - L) / L
        ω = bond_active[bond_id] * surface_correction_factor(correction, bond_id)
        b = ω * params.bc * ε * volume[j] .* Δxij / l
        update_add_vector!(b_int, i, b)
    end
    return nothing
end

function force_density_point!(storage::BBStorage, system::BondSystem, ::BBMaterial,
                              paramhandler::ParameterHandler, t, Δt, i)
    (; position, bond_length, bond_active, b_int) = storage
    (; bonds, correction, volume) = system
    params_i = get_params(paramhandler, i)
    for bond_id in each_bond_idx(system, i)
        bond = bonds[bond_id]
        j, L = bond.neighbor, bond.length
        Δxij = get_vector_diff(position, i, j)
        l = bond_length[bond_id]
        ε = (l - L) / L
        params_j = get_params(paramhandler, j)
        ω = bond_active[bond_id] * surface_correction_factor(correction, bond_id)
        b = ω * (params_i.bc + params_j.bc) / 2 * ε * volume[j] .* Δxij / l
        update_add_vector!(b_int, i, b)
    end
    return nothing
end

# Do not rely on any custom pre-stored properties here!
function strain_energy_density_point!(storage::AbstractStorage, system::BondSystem,
                                      ::BBMaterial, paramsetup::AbstractParameterSetup, i)
    (; position, bond_active, strain_energy_density) = storage
    (; bonds, correction, volume) = system
    params_i = get_params(paramsetup, i)
    Ψ = 0.0
    for bond_id in each_bond_idx(system, i)
        bond = bonds[bond_id]
        j, L = bond.neighbor, bond.length
        Δxij = get_vector_diff(position, i, j)
        l = norm(Δxij) # do not rely on the stored bond length here!
        ε = (l - L) / L
        params_j = get_params(paramsetup, j)
        ωij = bond_active[bond_id] * surface_correction_factor(correction, bond_id)
        bc = (params_i.bc + params_j.bc) / 2
        Ψ += 0.25 * ωij * bc * ε * ε * L * volume[j]
    end
    strain_energy_density[i] = Ψ
    return nothing
end

function export_field(::Val{:strain_energy_density}, mat::AbstractBondBasedMaterial,
                      system::BondSystem, storage::AbstractStorage,
                      paramsetup::AbstractParameterSetup, t)
    for i in each_point_idx(system)
        strain_energy_density_point!(storage, system, mat, paramsetup, i)
    end
    return storage.strain_energy_density
end
