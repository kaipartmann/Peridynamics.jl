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
"""
struct BBMaterial{Correction,DM} <: AbstractBondSystemMaterial{Correction}
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
    (; Gc, εc) = get_frac_params(mat.dmgmodel, p, δ, K)
    bc = 18 * K / (π * δ^4) # bond constant
    return StandardPointParameters(δ, rho, E, nu, G, K, λ, μ, Gc, εc, bc)
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
    bond_stretch::Vector{Float64}
    bond_active::Vector{Bool}
    @pointfield n_active_bonds::Vector{Int}
end

function init_field(::BBMaterial, ::AbstractTimeSolver, system::BondSystem,
                    ::Val{:bond_stretch})
    return zeros(get_n_bonds(system))
end

# Customized calc_failure to save the bond stretch ε for force density calculation
function calc_failure!(storage::BBStorage, system::BondSystem,
                       ::BBMaterial, ::CriticalStretch,
                       paramsetup::AbstractParameterSetup, i)
    (; εc) = get_params(paramsetup, i)
    (; position, n_active_bonds, bond_active, bond_stretch) = storage
    (; bonds) = system
    for bond_id in each_bond_idx(system, i)
        bond = bonds[bond_id]
        j, L = bond.neighbor, bond.length
        Δxij = get_vector_diff(position, i, j)
        l = norm(Δxij)
        ε = (l - L) / L
        bond_stretch[bond_id] = ε / l # note that this is  ε / l!
        if ε > εc && bond.fail_permit
            bond_active[bond_id] = false
        end
        n_active_bonds[i] += bond_active[bond_id]
    end
    return nothing
end

function force_density_point!(storage::BBStorage, system::BondSystem, ::BBMaterial,
                              params::StandardPointParameters, t, Δt, i)
    (; position, bond_stretch, bond_active, b_int) = storage
    (; bonds, correction, volume) = system
    for bond_id in each_bond_idx(system, i)
        bond = bonds[bond_id]
        j = bond.neighbor
        Δxij = get_vector_diff(position, i, j)
        ε = bond_stretch[bond_id]
        ω = bond_active[bond_id] * surface_correction_factor(correction, bond_id)
        b = ω * params.bc * ε * volume[j] .* Δxij
        update_add_vector!(b_int, i, b)
    end
    return nothing
end

function force_density_point!(storage::BBStorage, system::BondSystem, ::BBMaterial,
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
        b = ω * (params_i.bc + params_j.bc) / 2 * ε * volume[j] .* Δxij
        update_add_vector!(b_int, i, b)
    end
    return nothing
end
