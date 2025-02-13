
"""
    BBMaterial()
    BBMaterial{Correction}()

A material type used to assign the material of a [`Body`](@ref) with the standard bond-based
formulation of peridynamics.

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

# Allowed material parameters
When using [`material!`](@ref) on a [`Body`](@ref) with `BBMaterial`, then the following
parameters are allowed:
Material parameters:
- `horizon::Float64`: Radius of point interactions
- `rho::Float64`: Density
Elastic parameters
- `E::Float64`: Young's modulus
- `G::Float64`: Shear modulus
- `K::Float64`: Bulk modulus
- `lambda::Float64`: 1st Lamé parameter
- `mu::Float64`: 2nd Lamé parameter
Fracture parameters:
- `Gc::Float64`: Critical energy release rate
- `epsilon_c::Float64`: Critical strain

!!! note "Poisson's ratio and bond-based peridynamics"
    In bond-based peridynamics, the Poisson's ratio is limited to 1/4 for 3D simulations.
    Therefore the specification of this keyword is not allowed when using `material!`, as it
    is hardcoded to `nu = 1/4`.
    Therefore, only one additional elastic parameter is required.

# Allowed export fields
When specifying the `fields` keyword of [`Job`](@ref) for a [`Body`](@ref) with
`BBMaterial`, the following fields are allowed:
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
struct BBMaterial{Correction} <: AbstractBondSystemMaterial{Correction} end

BBMaterial() = BBMaterial{NoCorrection}()

struct BBPointParameters <: AbstractPointParameters
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

function BBPointParameters(mat::BBMaterial, p::Dict{Symbol,Any})
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
        msg *= "Please define only one or two fitting elastic parameters!\n"
        throw(ArgumentError(msg))
    end
    (; Gc, εc) = get_frac_params(p, δ, K)
    bc = 18 * K / (π * δ^4) # bond constant
    return BBPointParameters(δ, rho, E, nu, G, K, λ, μ, Gc, εc, bc)
end

@params BBMaterial BBPointParameters

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
    bond_active::Vector{Bool}
    @pointfield n_active_bonds::Vector{Int}
end

function force_density_point!(storage::BBStorage, system::BondSystem, ::BBMaterial,
                              params::BBPointParameters, i::Int)
    for bond_id in each_bond_idx(system, i)
        bond = system.bonds[bond_id]
        j, L = bond.neighbor, bond.length
        Δxij = get_coordinates_diff(storage, i, j)
        l = norm(Δxij)
        ε = (l - L) / L
        stretch_based_failure!(storage, system, bond, params, ε, i, bond_id)
        b_int = bond_failure(storage, bond_id) *
                surface_correction_factor(system.correction, bond_id) *
                params.bc * ε / l * system.volume[j] .* Δxij
        update_add_b_int!(storage, i, b_int)
    end
    return nothing
end

function force_density_point!(storage::BBStorage, system::BondSystem, ::BBMaterial,
                              paramhandler::ParameterHandler, i::Int)
    params_i = get_params(paramhandler, i)
    for bond_id in each_bond_idx(system, i)
        bond = system.bonds[bond_id]
        j, L = bond.neighbor, bond.length
        Δxij = get_coordinates_diff(storage, i, j)
        l = norm(Δxij)
        ε = (l - L) / L
        stretch_based_failure!(storage, system, bond, params_i, ε, i, bond_id)
        params_j = get_params(paramhandler, j)
        b_int = bond_failure(storage, bond_id) *
                surface_correction_factor(system.correction, bond_id) *
                (params_i.bc + params_j.bc) / 2 * ε / l * system.volume[j] .* Δxij
        update_add_b_int!(storage, i, b_int)
    end
    return nothing
end
