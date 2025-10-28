"""
    OSBMaterial(; kernel, dmgmodel)
    OSBMaterial{Correction}(; kernel, dmgmodel)

A material type used to assign the material of a [`Body`](@ref) with the ordinary
state-based formulation of peridynamics.

Possible correction methods are:
- [`NoCorrection`](@ref): No correction is applied. (default)
- [`EnergySurfaceCorrection`](@ref): The energy based surface correction method of
    Le and Bobaru (2018) is applied.

# Keywords
- `kernel::Function`: Kernel function used for weighting the interactions between points. \\
    (default: `linear_kernel`)
- `dmgmodel::AbstractDamageModel`: Damage model defining the damage behavior. \\
    (default: `CriticalStretch()`)

# Examples

```julia-repl
julia> mat = OSBMaterial()
OSBMaterial{NoCorrection}(dmgmodel=CriticalStretch())

julia> mat = OSBMaterial{EnergySurfaceCorrection}()
OSBMaterial{EnergySurfaceCorrection}(dmgmodel=CriticalStretch())
```

---

```julia
OSBMaterial{Correction,K,DM}
```

Material type for the ordinary state-based peridynamics formulation.

# Type Parameters
- `Correction`: A correction algorithm type. See the constructor docs for more informations.
- `K`: A kernel function type. See the constructor docs for more informations.
- `DM`: A damage model type. See the constructor docs for more informations.

# Fields
- `kernel::Function`: Kernel function used for weighting the interactions between points.
- `dmgmodel::AbstractDamageModel`: Damage model defining the damage behavior. See the
    constructor docs for more informations.

# Allowed material parameters
When using [`material!`](@ref) on a [`Body`](@ref) with `OSBMaterial`, then the following
parameters are allowed:
Material parameters:
- `horizon::Float64`: Radius of point interactions.
- `rho::Float64`: Density.
Elastic parameters:
- `E::Float64`: Young's modulus.
- `nu::Float64`: Poisson's ratio.
- `G::Float64`: Shear modulus.
- `K::Float64`: Bulk modulus.
- `lambda::Float64`: 1st Lamé parameter.
- `mu::Float64`: 2nd Lamé parameter.
Fracture parameters:
- `Gc::Float64`: Critical energy release rate.
- `epsilon_c::Float64`: Critical strain.

!!! note "Elastic parameters"
    Note that exactly two elastic parameters are required to specify a material.
    Please choose two out of the six allowed elastic parameters.

# Allowed export fields
When specifying the `fields` keyword of [`Job`](@ref) for a [`Body`](@ref) with
`OSBMaterial`, the following fields are allowed:
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
struct OSBMaterial{Correction,K,DM} <: AbstractBondSystemMaterial{Correction}
    kernel::K
    dmgmodel::DM
    function OSBMaterial{C}(kernel::K, dmgmodel::DM) where {C,K,DM}
        return new{C,K,DM}(kernel, dmgmodel)
    end
end

function OSBMaterial{C}(; kernel::F=linear_kernel,
                        dmgmodel::AbstractDamageModel=CriticalStretch()) where{C,F}
    return OSBMaterial{C}(kernel, dmgmodel)
end
OSBMaterial(; kwargs...) = OSBMaterial{NoCorrection}(; kwargs...)

@params OSBMaterial StandardPointParameters

@storage OSBMaterial struct OSBStorage <: AbstractStorage
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

function init_field(::OSBMaterial, ::AbstractTimeSolver, system::BondSystem, ::Val{:b_int})
    return zeros(3, get_n_points(system))
end

function init_field(::OSBMaterial, ::AbstractTimeSolver, system::BondSystem,
                    ::Val{:bond_length})
    return zeros(get_n_bonds(system))
end

function init_field(::OSBMaterial, ::AbstractTimeSolver, system::BondSystem,
                    ::Val{:strain_energy_density})
    return zeros(get_n_loc_points(system))
end

# Customized calc_failure to save the bond length for force density calculation
function calc_failure!(storage::OSBStorage, system::BondSystem,
                       ::OSBMaterial, ::CriticalStretch,
                       paramsetup::AbstractParameterSetup, i)
    (; εc) = get_params(paramsetup, i)
    (; position, n_active_bonds, bond_active, bond_length) = storage
    (; bonds) = system
    for bond_id in each_bond_idx(system, i)
        bond = bonds[bond_id]
        j, L = bond.neighbor, bond.length
        Δxij = get_vector_diff(position, i, j)
        l = norm(Δxij)
        bond_length[bond_id] = l # this is customized!
        ε = (l - L) / L
        if ε > εc && bond.fail_permit
            bond_active[bond_id] = false
        end
        n_active_bonds[i] += bond_active[bond_id]
    end
    return nothing
end

function force_density_point!(storage::OSBStorage, system::BondSystem, mat::OSBMaterial,
                              params::StandardPointParameters, t, Δt, i)
    wvol = calc_weighted_volume(storage, system, mat, params, i)
    iszero(wvol) && return nothing
    dil = calc_dilatation(storage, system, mat, params, wvol, i)
    (; position, bond_active, b_int, bond_length) = storage
    (; bonds, correction, volume) = system
    c1 = 15.0 * params.G / wvol
    c2 = dil * (3.0 * params.K / wvol - c1 / 3.0)
    for bond_id in each_bond_idx(system, i)
        bond = bonds[bond_id]
        j, L = bond.neighbor, bond.length
        Δxij = get_vector_diff(position, i, j)
        l = bond_length[bond_id]
        ωij = kernel(system, bond_id) * bond_active[bond_id]
        β = surface_correction_factor(correction, bond_id)
        p = ωij * β * (c2 * L + c1 * (l - L)) / l .* Δxij
        update_add_vector!(b_int, i, p .* volume[j])
        update_add_vector!(b_int, j, -p .* volume[i])
    end
    return nothing
end

function force_density_point!(storage::OSBStorage, system::BondSystem, mat::OSBMaterial,
                              paramhandler::ParameterHandler, t, Δt, i)
    params_i = get_params(paramhandler, i)
    wvol = calc_weighted_volume(storage, system, mat, params_i, i)
    iszero(wvol) && return nothing
    dil = calc_dilatation(storage, system, mat, params_i, wvol, i)
    (; position, bond_active, b_int, bond_length) = storage
    (; bonds, correction, volume) = system
    for bond_id in each_bond_idx(system, i)
        bond = bonds[bond_id]
        j, L = bond.neighbor, bond.length
        Δxij = get_vector_diff(position, i, j)
        l = bond_length[bond_id]
        params_j = get_params(paramhandler, j)
        c1 = 15.0 * (params_i.G + params_j.G) / (2 * wvol)
        c2 = dil * (3.0 * (params_i.K + params_j.K) / (2 * wvol) - c1 / 3.0)
        ωij = kernel(system, bond_id) * bond_active[bond_id]
        β = surface_correction_factor(correction, bond_id)
        p = ωij * β * (c2 * L + c1 * (l - L)) / l .* Δxij
        update_add_vector!(b_int, i, p .* volume[j])
        update_add_vector!(b_int, j, -p .* volume[i])
    end
    return nothing
end

function calc_weighted_volume(storage::OSBStorage, system::BondSystem, mat::OSBMaterial,
                              params::StandardPointParameters, i)
    wvol = 0.0
    for bond_id in each_bond_idx(system, i)
        bond = system.bonds[bond_id]
        j = bond.neighbor
        ΔXij = get_vector_diff(system.position, i, j)
        ΔXij_sq = dot(ΔXij, ΔXij)
        ωij = kernel(system, bond_id) * storage.bond_active[bond_id]
        β = surface_correction_factor(system.correction, bond_id)
        wvol += ωij * β * ΔXij_sq * system.volume[j]
    end
    return wvol
end

function calc_dilatation(storage::OSBStorage, system::BondSystem, mat::OSBMaterial,
                         params::StandardPointParameters, wvol, i)
    dil = 0.0
    c1 = 3.0 / wvol
    for bond_id in each_bond_idx(system, i)
        bond = system.bonds[bond_id]
        j, L = bond.neighbor, bond.length
        l = storage.bond_length[bond_id]
        ωij = kernel(system, bond_id) * storage.bond_active[bond_id]
        β = surface_correction_factor(system.correction, bond_id)
        dil += ωij * β * c1 * L * (l - L) * system.volume[j]
    end
    return dil
end

function strain_energy_density_point!(storage::AbstractStorage, system::BondSystem,
                                      mat::OSBMaterial, paramsetup::AbstractParameterSetup,
                                      i)
    (; bond_active, bond_length, strain_energy_density) = storage
    (; bonds, correction, volume) = system
    params_i = get_params(paramsetup, i)
    wvol = calc_weighted_volume(storage, system, mat, params_i, i)
    iszero(wvol) && return nothing
    dil = calc_dilatation(storage, system, mat, params_i, wvol, i)
    Ψvol = 0.5 * params_i.K * dil^2
    Ψdev = 0.0
    for bond_id in each_bond_idx(system, i)
        bond = bonds[bond_id]
        j, L = bond.neighbor, bond.length
        l = bond_length[bond_id]
        e = l - L
        edev = e - 1/3 * dil * L
        ωij = kernel(system, bond_id) * bond_active[bond_id]
        β = surface_correction_factor(correction, bond_id)
        params_j = get_params(paramsetup, j)
        G = (params_i.G + params_j.G) / 2
        cdev = 15.0 * G / (2 * wvol)
        Ψdev += cdev * ωij * β * edev * edev * volume[j]
    end
    Ψ = Ψvol + Ψdev
    strain_energy_density[i] = Ψ
    return nothing
end

function export_field(::Val{:strain_energy_density}, mat::OSBMaterial, system::BondSystem,
                      storage::AbstractStorage, paramsetup::AbstractParameterSetup, t)
    for i in each_point_idx(system)
        strain_energy_density_point!(storage, system, mat, paramsetup, i)
    end
    return storage.strain_energy_density
end
