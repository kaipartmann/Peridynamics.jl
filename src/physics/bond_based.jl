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

function StandardPointParameters(mat::BBMaterial, p::Dict{Symbol,Any})
    (; δ, rho, E, nu, G, K, λ, μ) = get_required_point_parameters_bb(mat, p)
    (; Gc, εc) = get_frac_params(mat.dmgmodel, p, δ, K)
    bc = 18 * K / (π * δ^4) # bond constant
    return StandardPointParameters(δ, rho, E, nu, G, K, λ, μ, Gc, εc, bc)
end

function get_required_point_parameters_bb(mat::AbstractBondBasedMaterial,
                                          p::Dict{Symbol,Any})
    par = get_given_elastic_params(p)
    (; E, nu, G, K, λ, μ) = par
    if isfinite(nu) && !isapprox(nu, 0.25)
        msg = "Bond-based peridynamics has a limitation on the Poisson's ratio!\n"
        msg *= "With BBMaterial, no other values than nu=0.25 are allowed!\n"
        throw(ArgumentError(msg))
    elseif !isfinite(nu) && length(findall(isfinite, par)) == 1
        p[:nu] = 0.25
    end
    (; δ, rho, E, nu, G, K, λ, μ) = get_required_point_parameters(mat, p)
    if !isapprox(nu, 0.25)
        msg = "Bond-based peridynamics has a limitation on the Poisson's ratio!\n"
        msg *= "With BBMaterial, no other values than nu=0.25 are allowed!\n"
        msg *= "The submitted parameter combination results in an illegal value for nu!\n"
        msg *= "Please define either only one or two fitting elastic parameters!\n"
        throw(ArgumentError(msg))
    end
    return (; δ, rho, E, nu, G, K, λ, μ)
end

@params BBMaterial StandardPointParameters

@storage BBMaterial struct BBStorage <: AbstractStorage
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

function init_field(::AbstractBondBasedMaterial, ::AbstractTimeSolver, system::BondSystem,
                    ::Val{:bond_length})
    return zeros(get_n_bonds(system))
end

function init_field(::AbstractBondBasedMaterial, ::AbstractTimeSolver, system::BondSystem,
                    ::Val{:strain_energy_density})
    return zeros(get_n_loc_points(system))
end

# Customized calc_failure to save the bond stretch ε for force density calculation
function calc_failure!(storage::AbstractStorage, system::BondSystem,
                       ::AbstractBondBasedMaterial, ::CriticalStretch,
                       paramsetup::AbstractParameterSetup, i)
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
                              params::StandardPointParameters, t, Δt, i)
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

function strain_energy_density_point!(storage::AbstractStorage, system::BondSystem,
                                      ::BBMaterial, paramsetup::AbstractParameterSetup, i)
    (; bond_active, bond_length, strain_energy_density) = storage
    (; bonds, correction, volume) = system
    params_i = get_params(paramsetup, i)
    Ψ = 0.0
    for bond_id in each_bond_idx(system, i)
        bond = bonds[bond_id]
        j, L = bond.neighbor, bond.length
        l = bond_length[bond_id]
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
