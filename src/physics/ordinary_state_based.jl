"""
    OSBMaterial{Correction} <: AbstractMaterial

Material type for ordinary state-based peridynamic simulations

# Type Parameters

- `Correction`: Applied surface correction method
    Available methods: - `NoCorrection`: No surface corrections are applied (default)
                       - `EnergySurfaceCorrection`: Surface correction factors are applied
                           based on the strain energy density

# Allowed material parameters

- `horizon::Float64`: Radius of point interactions
- `rho::Float64`: Density
- `E::Float64`: Young's modulus
- `nu::Float64`: Poisson's ratio
- `Gc::Float64`: Critical energy release rate
- `epsilon_c::Float64`: Critical strain

# Allowed export fields

- `position::Matrix{Float64}`: Position of each point
- `displacement::Matrix{Float64}`: Displacement of each point
- `velocity::Matrix{Float64}`: Velocity of each point
- `velocity_half::Matrix{Float64}`: Velocity parameter for Verlet time solver
- `acceleration::Matrix{Float64}`: Acceleration of each point
- `b_int::Matrix{Float64}`: Internal force density of each point
- `b_ext::Matrix{Float64}`: External force density of each point
- `damage::Vector{Float64}`: Damage of each point
- `n_active_bonds::Vector{Int}`: Number of intact bonds for each point

# Examples

```julia-repl
julia> mat = OSBMaterial()
OSBMaterial{NoCorrection}()

julia> mat = OSBMaterial{EnergySurfaceCorrection}()
OSBMaterial{EnergySurfaceCorrection}()
```

"""
struct OSBMaterial{Correction} <: AbstractBondSystemMaterial{Correction} end

OSBMaterial() = OSBMaterial{NoCorrection}()

"""
TODO
"""
struct OSBPointParameters <: AbstractPointParameters
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

function OSBPointParameters(::OSBMaterial, p::Dict{Symbol,Any})
    δ = get_horizon(p)
    rho = get_density(p)
    E, nu, G, K, λ, μ = get_elastic_params(p)
    Gc, εc = get_frac_params(p, δ, K)
    bc = 18 * K / (π * δ^4)
    return OSBPointParameters(δ, rho, E, nu, G, K, λ, μ, Gc, εc, bc)
end

@params OSBMaterial OSBPointParameters

struct OSBVerletStorage <: AbstractStorage
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

function OSBVerletStorage(::OSBMaterial, ::VelocityVerlet, system::BondSystem, ch)
    n_loc_points = length(ch.loc_points)
    position = copy(system.position)
    displacement = zeros(3, n_loc_points)
    velocity = zeros(3, n_loc_points)
    velocity_half = zeros(3, n_loc_points)
    acceleration = zeros(3, n_loc_points)
    b_int = zeros(3, length(ch.point_ids))
    b_ext = zeros(3, n_loc_points)
    damage = zeros(n_loc_points)
    bond_active = ones(Bool, length(system.bonds))
    n_active_bonds = copy(system.n_neighbors)
    s = OSBVerletStorage(position, displacement, velocity, velocity_half, acceleration,
                         b_int, b_ext, damage, bond_active, n_active_bonds)
    return s
end

@storage OSBMaterial VelocityVerlet OSBVerletStorage

@loc_to_halo_fields OSBVerletStorage :position
@halo_to_loc_fields OSBVerletStorage :b_int

struct OSBRelaxationStorage <: AbstractStorage
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

function OSBRelaxationStorage(::OSBMaterial, ::DynamicRelaxation, system::BondSystem, ch)
    n_loc_points = length(ch.loc_points)
    position = copy(system.position)
    displacement = zeros(3, n_loc_points)
    velocity = zeros(3, n_loc_points)
    velocity_half = zeros(3, n_loc_points)
    velocity_half_old = zeros(3, n_loc_points)
    b_int = zeros(3, length(ch.point_ids))
    b_int_old = zeros(3, n_loc_points)
    b_ext = zeros(3, n_loc_points)
    density_matrix = zeros(3, n_loc_points)
    damage = zeros(n_loc_points)
    bond_active = ones(Bool, length(system.bonds))
    n_active_bonds = copy(system.n_neighbors)
    s = OSBRelaxationStorage(position, displacement, velocity, velocity_half,
                             velocity_half_old, b_int, b_int_old, b_ext, density_matrix,
                             damage, bond_active, n_active_bonds)
    return s
end

@storage OSBMaterial DynamicRelaxation OSBRelaxationStorage

@loc_to_halo_fields OSBRelaxationStorage :position
@halo_to_loc_fields OSBRelaxationStorage :b_int

const OSBStorage = Union{OSBVerletStorage,OSBRelaxationStorage}

function force_density_point!(storage::OSBStorage, system::BondSystem, ::OSBMaterial,
                              params::OSBPointParameters, i::Int)
    # weighted volume
    wvol = calc_weighted_volume(storage, system, params, i)
    iszero(wvol) && return nothing
    # dilatation
    dil = calc_dilatation(storage, system, params, wvol, i)
    # force density
    c1 = 15.0 * params.G / wvol
    c2 = dil * (3.0 * params.K / wvol - c1 / 3.0)
    for bond_id in each_bond_idx(system, i)
        bond = system.bonds[bond_id]
        j, L = bond.neighbor, bond.length
        Δxijx = storage.position[1, j] - storage.position[1, i]
        Δxijy = storage.position[2, j] - storage.position[2, i]
        Δxijz = storage.position[3, j] - storage.position[3, i]
        l = sqrt(Δxijx * Δxijx + Δxijy * Δxijy + Δxijz * Δxijz)
        ε = (l - L) / L

        # failure mechanism
        if ε > params.εc && bond.fail_permit
            storage.bond_active[bond_id] = false
        end
        storage.n_active_bonds[i] += storage.bond_active[bond_id]

        # update of force density
        scfactor = surface_correction_factor(system.correction, bond_id)
        ωij = (1 + params.δ / L) * storage.bond_active[bond_id] * scfactor
        temp = ωij * (c2 * L + c1 * (l - L)) / l
        storage.b_int[1, i] += temp * Δxijx * system.volume[j]
        storage.b_int[2, i] += temp * Δxijy * system.volume[j]
        storage.b_int[3, i] += temp * Δxijz * system.volume[j]
        storage.b_int[1, j] -= temp * Δxijx * system.volume[i]
        storage.b_int[2, j] -= temp * Δxijy * system.volume[i]
        storage.b_int[3, j] -= temp * Δxijz * system.volume[i]
    end
    return nothing
end

function calc_weighted_volume(storage::OSBStorage, system::BondSystem,
                              params::OSBPointParameters, i::Int)
    wvol = 0.0
    for bond_id in each_bond_idx(system, i)
        bond = system.bonds[bond_id]
        j, L = bond.neighbor, bond.length
        ΔXijx = system.position[1, j] - system.position[1, i]
        ΔXijy = system.position[2, j] - system.position[2, i]
        ΔXijz = system.position[3, j] - system.position[3, i]
        ΔXij_sq = ΔXijx * ΔXijx + ΔXijy * ΔXijy + ΔXijz * ΔXijz
        scfactor = surface_correction_factor(system.correction, bond_id)
        ωij = (1 + params.δ / L) * storage.bond_active[bond_id] * scfactor
        wvol += ωij * ΔXij_sq * system.volume[j]
    end
    return wvol
end

function calc_dilatation(storage::OSBStorage, system::BondSystem,
                         params::OSBPointParameters, wvol::Float64, i::Int)
    dil = 0.0
    c1 = 3.0 / wvol
    for bond_id in each_bond_idx(system, i)
        bond = system.bonds[bond_id]
        j, L = bond.neighbor, bond.length
        Δxijx = storage.position[1, j] - storage.position[1, i]
        Δxijy = storage.position[2, j] - storage.position[2, i]
        Δxijz = storage.position[3, j] - storage.position[3, i]
        l = sqrt(Δxijx * Δxijx + Δxijy * Δxijy + Δxijz * Δxijz)
        scfactor = surface_correction_factor(system.correction, bond_id)
        ωij = (1 + params.δ / L) * storage.bond_active[bond_id] * scfactor
        dil += ωij * c1 * L * (l - L) * system.volume[j]
    end
    return dil
end
