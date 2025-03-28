"""
    NSCMaterial(; kernel, model, dmgmodel, maxdmg)

A material type used to assign the material of a [`Body`](@ref) with the local continuum
consistent (correspondence) formulation of non-ordinary state-based peridynamics.

# Keywords
- `kernel::Function`: Kernel function used for weighting the interactions between points. \\
    (default: `cubic_b_spline_kernel`) \\
    The following kernels can be used:
    - [`cubic_b_spline_kernel`](@ref)
- `model::AbstractConstitutiveModel`: Constitutive model defining the material behavior. \\
    (default: `LinearElastic()`) \\
    The following models can be used:
    - [`LinearElastic`](@ref)
    - [`NeoHooke`](@ref)
    - [`MooneyRivlin`](@ref)
    - [`SaintVenantKirchhoff`](@ref)
- `dmgmodel::AbstractDamageModel`: Damage model defining the damage behavior. \\
    (default: [`CriticalStretch`](@ref))
- `maxdmg::Float64`: Maximum value of damage a point is allowed to obtain. If this value is
    exceeded, all bonds of that point are broken because the deformation gradient would then
    possibly contain `NaN` values. \\
    (default: `0.85`)

!!! note "Stability of fracture simulations"
    This formulation is known to be not suitable for fracture simulations without
    stabilization of the zero-energy modes. Therefore be careful when doing fracture
    simulations and try out different parameters for `maxdmg` and `zem`.

# Examples

```julia-repl
julia> mat = NSCMaterial()
NSCMaterial{LinearElastic, typeof(linear_kernel), CriticalStretch}(maxdmg=0.85)
```

---

```julia
NSCMaterial{CM,K,DM}
```

Material type for the naturally stabilized local continuum consistent (correspondence) formulation of
non-ordinary state-based peridynamics.

# Type Parameters
- `CM`: A constitutive model type. See the constructor docs for more informations.
- `K`: A kernel function type. See the constructor docs for more informations.
- `DM`: A damage model type. See the constructor docs for more informations.

# Fields
- `kernel::Function`: Kernel function used for weighting the interactions between points.
    See the constructor docs for more informations.
- `model::AbstractConstitutiveModel`: Constitutive model defining the material behavior. See
    the constructor docs for more informations.
- `dmgmodel::AbstractDamageModel`: Damage model defining the damage behavior. See the
    constructor docs for more informations.
- `maxdmg::Float64`: Maximum value of damage a point is allowed to obtain. See the
    constructor docs for more informations.

# Allowed material parameters
When using [`material!`](@ref) on a [`Body`](@ref) with `NSCMaterial`, then the following
parameters are allowed:
Material parameters:
- `horizon::Float64`: Radius of point interactions
- `rho::Float64`: Density
Elastic parameters:
- `E::Float64`: Young's modulus
- `nu::Float64`: Poisson's ratio
- `G::Float64`: Shear modulus
- `K::Float64`: Bulk modulus
- `lambda::Float64`: 1st Lamé parameter
- `mu::Float64`: 2nd Lamé parameter
Fracture parameters:
- `Gc::Float64`: Critical energy release rate
- `epsilon_c::Float64`: Critical strain

!!! note "Elastic parameters"
    Note that exactly two elastic parameters are required to specify a material.
    Please choose two out of the six allowed elastic parameters.

# Allowed export fields
When specifying the `fields` keyword of [`Job`](@ref) for a [`Body`](@ref) with
`NSCMaterial`, the following fields are allowed:
- `position::Matrix{Float64}`: Position of each point
- `displacement::Matrix{Float64}`: Displacement of each point
- `velocity::Matrix{Float64}`: Velocity of each point
- `velocity_half::Matrix{Float64}`: Velocity parameter for Verlet time solver
- `acceleration::Matrix{Float64}`: Acceleration of each point
- `b_int::Matrix{Float64}`: Internal force density of each point
- `b_ext::Matrix{Float64}`: External force density of each point
- `damage::Vector{Float64}`: Damage of each point
- `n_active_bonds::Vector{Int}`: Number of intact bonds of each point
- `stress::Matrix{Float64}`: Stress tensor of each point
- `von_mises_stress::Vector{Float64}`: Von Mises stress of each point
"""
struct NSCMaterial{CM,K,DM} <: AbstractNSCMaterial{CM,NoCorrection}
    kernel::K
    constitutive_model::CM
    dmgmodel::DM
    maxdmg::Float64
    function NSCMaterial(kernel::K, cm::CM, dmgmodel::DM, maxdmg::Real) where {CM,K,DM}
        return new{CM,K,DM}(kernel, cm, dmgmodel, maxdmg)
    end
end

function NSCMaterial(; kernel::Function=cubic_b_spline_kernel,
                     model::AbstractConstitutiveModel=LinearElastic(),
                     dmgmodel::AbstractDamageModel=CriticalStretch(), maxdmg::Real=0.85)
    if kernel !== cubic_b_spline_kernel
        msg = "The kernel $kernel is not recommended for use with NSCMaterial."
        msg *= "Use the `cubic_b_spline_kernel` function instead.\n"
        throw(ArgumentError(msg))
    end
    return NSCMaterial(kernel, model, dmgmodel, maxdmg)
end

function Base.show(io::IO, @nospecialize(mat::NSCMaterial))
    print(io, typeof(mat))
    print(io, msg_fields_in_brackets(mat, (:maxdmg,)))
    return nothing
end

@params NSCMaterial StandardPointParameters

@storage NSCMaterial struct NSCStorage
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
    bond_active::Vector{Bool}
    @pointfield n_active_bonds::Vector{Int}
    @pointfield damage_changed::Vector{Bool}
    @pointfield stress::Matrix{Float64}
    @pointfield von_mises_stress::Vector{Float64}
    @lthfield defgrad::Matrix{Float64}
    @lthfield weighted_volume::Vector{Float64}
    gradient_weight::Matrix{Float64}
    first_piola_kirchhoff::Matrix{Float64}
end

function init_field(::AbstractNSCMaterial, ::AbstractTimeSolver, system::BondSystem,
                    ::Val{:b_int})
    return zeros(3, get_n_points(system))
end

function init_field(::AbstractNSCMaterial, ::AbstractTimeSolver, system::BondSystem,
                    ::Val{:damage_changed})
    return ones(Bool, get_n_loc_points(system))
end

function init_field(::AbstractNSCMaterial, ::AbstractTimeSolver, system::BondSystem,
                    ::Val{:stress})
    return zeros(9, get_n_loc_points(system))
end

function init_field(::AbstractNSCMaterial, ::AbstractTimeSolver, system::BondSystem,
                    ::Val{:von_mises_stress})
    return zeros(get_n_loc_points(system))
end

function init_field(::AbstractNSCMaterial, ::AbstractTimeSolver, system::BondSystem,
                    ::Val{:defgrad})
    return zeros(9, get_n_points(system))
end

function init_field(::AbstractNSCMaterial, ::AbstractTimeSolver, system::BondSystem,
                    ::Val{:weighted_volume})
    return zeros(get_n_points(system))
end

function init_field(::AbstractNSCMaterial, ::AbstractTimeSolver, system::BondSystem,
                    ::Val{:gradient_weight})
    return zeros(3, get_n_bonds(system))
end

function init_field(::AbstractNSCMaterial, ::AbstractTimeSolver, system::BondSystem,
                    ::Val{:first_piola_kirchhoff})
    return zeros(9, get_n_bonds(system))
end

function initialize!(chunk::BodyChunk{<:BondSystem,<:AbstractNSCMaterial})
    chunk.storage.damage_changed .= true
    return nothing
end

function calc_force_density!(dh::ThreadsBodyDataHandler{<:BondSystem,<:AbstractNSCMaterial},
                             t, Δt)
    (; chunks) = dh
    after_fields = nsc_lth_after_fields(chunks[1].mat)
    pre_fields = filter(x -> !in(x, after_fields), loc_to_halo_fields(chunks[1].storage))
    @threads :static for chunk_id in eachindex(chunks)
        exchange_loc_to_halo!(dh, chunk_id, pre_fields)
    end
    @threads :static for chunk_id in eachindex(chunks)
        calc_weights_and_defgrad!(chunks[chunk_id], t, Δt)
    end
    @threads :static for chunk_id in eachindex(chunks)
        exchange_loc_to_halo!(dh, chunk_id, after_fields)
    end
    @threads :static for chunk_id in eachindex(chunks)
        calc_force_density!(chunks[chunk_id], t, Δt)
    end
    @threads :static for chunk_id in eachindex(chunks)
        exchange_halo_to_loc!(dh, chunk_id)
    end
    return nothing
end

nsc_lth_after_fields(::NSCMaterial) = (:defgrad, :weighted_volume)

function calc_force_density!(dh::MPIBodyDataHandler{<:BondSystem,<:AbstractNSCMaterial}, t,
                             Δt)
    (; chunk) = dh
    exchange_loc_to_halo!(dh, :position)
    calc_weights_and_defgrad!(chunk, t, Δt)
    exchange_loc_to_halo!(dh, (:defgrad, :weighted_volume))
    calc_force_density!(chunk, t, Δt)
    exchange_halo_to_loc!(dh)
    return nothing
end

function calc_weights_and_defgrad!(chunk::BodyChunk{<:BondSystem,<:AbstractNSCMaterial}, t,
                                   Δt)
    (; system, mat, paramsetup, storage) = chunk
    (; dmgmodel) = mat
    storage.n_active_bonds .= 0
    for i in each_point_idx(system)
        calc_failure!(storage, system, mat, dmgmodel, paramsetup, i)
        calc_damage!(storage, system, mat, dmgmodel, paramsetup, i)
        calc_weights_and_defgrad!(storage, system, mat, paramsetup, t, Δt, i)
    end
    return nothing
end

function calc_damage!(storage::AbstractStorage, system::AbstractBondSystem,
                      mat::AbstractNSCMaterial, dmgmodel::AbstractDamageModel,
                      paramsetup::AbstractParameterSetup, i)
    (; n_neighbors) = system
    (; n_active_bonds, damage, damage_changed) = storage
    old_damage = damage[i]
    new_damage = 1 - n_active_bonds[i] / n_neighbors[i]
    # if new_damage > old_damage
    #     damage_changed[i] = true
    # else
    #     damage_changed[i] = true # TODO: somehow this makes an error if set to `false`!!!
    # end
    damage_changed[i] = true
    damage[i] = new_damage
    return nothing
end

function damage_changed(::AbstractNSCMaterial, storage::AbstractStorage, i)
    return storage.damage_changed[i]
end

function calc_weights_and_defgrad!(storage::AbstractStorage, system::AbstractBondSystem,
                                   mat::AbstractNSCMaterial,
                                   paramsetup::AbstractParameterSetup, t, Δt, i)
    params = get_params(paramsetup, i)
    Kinv = nsc_defgrad!(storage, system, mat, params, t, Δt, i)
    if damage_changed(mat, storage, i)
        nsc_weights!(storage, system, mat, params, Kinv, t, Δt, i)
    end
    return nothing
end

function nsc_defgrad!(storage::AbstractStorage, system::AbstractBondSystem,
                      mat::AbstractNSCMaterial, params::AbstractPointParameters, t, Δt, i)
    (; bonds, volume) = system
    (; bond_active, defgrad, weighted_volume) = storage

    K = zero(SMatrix{3,3,Float64,9})
    _F = zero(SMatrix{3,3,Float64,9})
    wi = 0.0
    for bond_id in each_bond_idx(system, i)
        bond = bonds[bond_id]
        j = bond.neighbor
        ΔXij = get_vector_diff(system.position, i, j)
        Δxij = get_vector_diff(storage.position, i, j)
        ωij = kernel(system, bond_id) * bond_active[bond_id]
        temp = ωij * volume[j]
        K += temp * (ΔXij * ΔXij')
        _F += temp * (Δxij * ΔXij')
        wi += temp
    end
    Kinv = inv(K)
    F = _F * Kinv
    update_tensor!(defgrad, i, F)
    weighted_volume[i] = wi

    return Kinv
end

function nsc_weights!(storage::AbstractStorage, system::AbstractBondSystem,
                      mat::AbstractNSCMaterial, params::AbstractPointParameters,
                      Kinv::SMatrix{3,3,T,9}, t, Δt, i) where {T}
    (; bonds, volume) = system
    (; bond_active, gradient_weight) = storage

    for bond_id in each_bond_idx(system, i)
        bond = bonds[bond_id]
        j = bond.neighbor
        ΔXij = get_vector_diff(system.position, i, j)
        ωij = kernel(system, bond_id) * bond_active[bond_id]
        temp = ωij * volume[j]
        Ψ = temp * (Kinv * ΔXij)
        update_vector!(gradient_weight, bond_id, Ψ)
    end

    return nothing
end

function calc_force_density!(chunk::BodyChunk{<:BondSystem,<:AbstractNSCMaterial}, t, Δt)
    (; system, mat, paramsetup, storage) = chunk
    storage.b_int .= 0
    for i in each_point_idx(chunk)
        force_density_point!(storage, system, mat, paramsetup, t, Δt, i)
    end
    nancheck(chunk, t)
    return nothing
end

function force_density_point!(storage::AbstractStorage, system::BondSystem,
                              mat::AbstractNSCMaterial, params::AbstractPointParameters, t,
                              Δt, i)
    ∑P = nsc_stress_integral!(storage, system, mat, params, t, Δt, i)
    nsc_force_density!(storage, system, mat, params, ∑P, t, Δt, i)
    return nothing
end

function nsc_stress_integral!(storage::AbstractStorage, system::AbstractBondSystem,
                              mat::AbstractNSCMaterial, params::AbstractPointParameters, t,
                              Δt, i)
    (; bonds, volume) = system
    (; bond_active, gradient_weight, defgrad, weighted_volume) = storage
    Fi = get_tensor(defgrad, i)
    too_much_damage!(storage, system, mat, Fi, i) && return nothing
    wi = weighted_volume[i]
    ∑P = SMatrix{3,3,Float64,9}(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    for bond_id in each_bond_idx(system, i)
        if bond_active[bond_id]
            bond = bonds[bond_id]
            j, L = bond.neighbor, bond.length
            ΔXij = get_vector_diff(system.position, i, j)
            Δxij = get_vector_diff(storage.position, i, j)
            Fj = get_tensor(defgrad, j)
            Fb = 0.5 * (Fi + Fj)
            ΔXijLL = ΔXij' / (L * L)
            ΔFij = (Δxij - Fb * ΔXij) * ΔXijLL
            Fij = Fb + ΔFij
            Pij = calc_first_piola_kirchhoff!(storage, mat, params, Fij, bond_id)
            Tempij = I - ΔXij * ΔXijLL
            wj = weighted_volume[j]
            ϕ = (wi > 0 && wj > 0) ? (0.5 / wi + 0.5 / wj) : 0.0
            ω̃ij = kernel(system, bond_id) * ϕ
            ∑Pij = ω̃ij * (Pij * Tempij)
            ∑P += ∑Pij
        end
    end
    return ∑P
end

function nsc_force_density!(storage::AbstractStorage, system::AbstractBondSystem,
                            mat::AbstractNSCMaterial, params::AbstractPointParameters,
                            ∑P, t, Δt, i)
    (; bonds, volume) = system
    (; bond_active, gradient_weight, weighted_volume) = storage
    wi = weighted_volume[i]
    for bond_id in each_bond_idx(system, i)
        if bond_active[bond_id]
            bond = bonds[bond_id]
            j, L = bond.neighbor, bond.length
            ΔXij = get_vector_diff(system.position, i, j)
            Pij = get_tensor(storage.first_piola_kirchhoff, bond_id)
            Φij = get_vector(gradient_weight, bond_id)
            wj = weighted_volume[j]
            ϕ = (wi > 0 && wj > 0) ? (0.5 / wi + 0.5 / wj) : 0.0
            ω̃ij = kernel(system, bond_id) * ϕ
            tij = ω̃ij / (L * L) * (Pij * ΔXij) + ∑P * Φij
            update_add_vector!(storage.b_int, i, tij * volume[j])
            update_add_vector!(storage.b_int, j, -tij * volume[i])
        end
    end
    return nothing
end

function calc_first_piola_kirchhoff!(storage::NSCStorage, mat::NSCMaterial,
                                     params::StandardPointParameters, F, bond_id)
    P = first_piola_kirchhoff(mat.constitutive_model, storage, params, F)
    update_tensor!(storage.first_piola_kirchhoff, bond_id, P)
    return P
end

function too_much_damage!(storage::AbstractStorage, system::BondSystem,
                          mat::AbstractNSCMaterial, F, i)
    if storage.damage[i] > mat.maxdmg || containsnan(F)
        # kill all bonds of this point
        storage.bond_active[each_bond_idx(system, i)] .= false
        return true
    end
    return false
end
