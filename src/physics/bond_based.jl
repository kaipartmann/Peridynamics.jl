
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
- `horizon::Float64`: Radius of point interactions
- `rho::Float64`: Density
- `E::Float64`: Young's modulus
- `Gc::Float64`: Critical energy release rate
- `epsilon_c::Float64`: Critical strain

!!! note "Poisson's ratio and bond-based peridynamics"
    In bond-based peridynamics, the Poisson's ratio is limited to 1/4 for 3D simulations.
    Therefore the specification of this keyword is not allowed when using `material!`, as it
    is hardcoded to `nu = 1/4`.

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

function BBPointParameters(::BBMaterial, p::Dict{Symbol,Any})
    δ = get_horizon(p)
    rho = get_density(p)
    if haskey(p, :nu)
        msg = "keyword `nu` is not allowed for BBMaterial!\n"
        msg *= "Bond-based peridynamics has a limitation on the possion ratio.\n"
        msg *= "Therefore, when using BBMaterial, `nu` is hardcoded as 1/4.\n"
        throw(ArgumentError(msg))
    else
        p[:nu] = 0.25
    end
    E, nu, G, K, λ, μ = get_elastic_params(p)
    Gc, εc = get_frac_params(p, δ, K)
    bc = 18 * K / (π * δ^4) # bond constant
    return BBPointParameters(δ, rho, E, nu, G, K, λ, μ, Gc, εc, bc)
end

@params BBMaterial BBPointParameters

struct BBVerletStorage <: AbstractStorage
    position::Matrix{Float64}
    displacement::Matrix{Float64}
    velocity::Matrix{Float64}
    velocity_half::Matrix{Float64}
    acceleration::Matrix{Float64}
    b_int::Matrix{Float64}
    b_ext::Matrix{Float64}
    damage::Vector{Float64}
    bond_active::Vector{Bool}
    n_active_bonds::Vector{Int}
end

function BBVerletStorage(::BBMaterial, ::VelocityVerlet, system::BondSystem, ch)
    n_loc_points = length(ch.loc_points)
    position = copy(system.position)
    displacement = zeros(3, n_loc_points)
    velocity = zeros(3, n_loc_points)
    velocity_half = zeros(3, n_loc_points)
    acceleration = zeros(3, n_loc_points)
    b_int = zeros(3, n_loc_points)
    b_ext = zeros(3, n_loc_points)
    damage = zeros(n_loc_points)
    bond_active = ones(Bool, length(system.bonds))
    n_active_bonds = copy(system.n_neighbors)
    s = BBVerletStorage(position, displacement, velocity, velocity_half, acceleration,
                        b_int, b_ext, damage, bond_active, n_active_bonds)
    return s
end

@storage BBMaterial VelocityVerlet BBVerletStorage

@loc_to_halo_fields BBVerletStorage :position

struct BBRelaxationStorage <: AbstractStorage
    position::Matrix{Float64}
    displacement::Matrix{Float64}
    velocity::Matrix{Float64}
    velocity_half::Matrix{Float64}
    velocity_half_old::Matrix{Float64}
    b_int::Matrix{Float64}
    b_int_old::Matrix{Float64}
    b_ext::Matrix{Float64}
    density_matrix::Matrix{Float64}
    damage::Vector{Float64}
    bond_active::Vector{Bool}
    n_active_bonds::Vector{Int}
end

function BBRelaxationStorage(::BBMaterial, ::DynamicRelaxation, system::BondSystem, ch)
    n_loc_points = length(ch.loc_points)
    position = copy(system.position)
    displacement = zeros(3, n_loc_points)
    velocity = zeros(3, n_loc_points)
    velocity_half = zeros(3, n_loc_points)
    velocity_half_old = zeros(3, n_loc_points)
    b_int = zeros(3, n_loc_points)
    b_int_old = zeros(3, n_loc_points)
    b_ext = zeros(3, n_loc_points)
    density_matrix = zeros(3, n_loc_points)
    damage = zeros(n_loc_points)
    bond_active = ones(Bool, length(system.bonds))
    n_active_bonds = copy(system.n_neighbors)
    s = BBRelaxationStorage(position, displacement, velocity, velocity_half,
                            velocity_half_old, b_int, b_int_old, b_ext, density_matrix,
                            damage, bond_active, n_active_bonds)
    return s
end

@storage BBMaterial DynamicRelaxation BBRelaxationStorage

@loc_to_halo_fields BBRelaxationStorage :position

const BBStorage = Union{BBVerletStorage,BBRelaxationStorage}

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
