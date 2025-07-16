"""
    RKCMaterial(; kernel, model, dmgmodel, maxdmg, reprkernel, regfactor)

A material type used to assign the material of a [`Body`](@ref) with a reproducing kernel
peridynamics (correspondence) formulation with bond-associated quadrature integration at
the center of the bonds.

# Keywords
- `kernel::Function`: Kernel function used for weighting the interactions between points. \\
    (default: `cubic_b_spline_kernel`) \\
    The following kernels can be used:
    - [`const_one_kernel`](@ref)
    - [`linear_kernel`](@ref)
    - [`cubic_b_spline_kernel`](@ref)
    - [`cubic_b_spline_kernel_norm`](@ref)
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
    (default: `1.0`)
- `reprkernel::Symbol`: A kernel function used for the reproducing kernel approximation
    of the moment matrix. This kernel is used to calculate the moment matrix and the
    gradient weights, which are used to approximate the deformation gradient. \\
    (default: `:C1`) \\
    The following kernels can be used:
    - `:C1`: Linear reproducing kernel basis [x, y, z] with first-order accuracy
        for correspondence formulation. This is the default kernel.
    - `:RK1`: First-order reproducing kernel basis [1, x, y, z] with constant term
        for enhanced stability in uniform deformation fields.
    - `:RK2`: Second-order reproducing kernel basis [1, x, y, z, x², y², z²] with
        diagonal quadratic terms for improved accuracy in curved deformation fields.
    - `:PD2`: Second-order peridynamic basis [x, y, z, x², xy, xz, y², yz, z²] with
        full quadratic terms but without constant term.
- `regfactor::Float64`: A regularization factor used to stabilize the inversion of the
    moment matrix. The product of `regfactor * δ^3` is used for the regularization with
    [`invreg`](@ref). The regularization helps to ensure numerical stability by mitigating
    issues such as ill-conditioning or singularity in the matrix operations.\\
    (default: `1e-13`)

# Examples

```julia-repl
julia> mat = RKCMaterial()
RKCMaterial{LinearElastic, typeof(linear_kernel), CriticalStretch}(maxdmg=1.0)
```

---

```julia
RKCMaterial{CM,K,DM}
```

Material type for a reproducing kernel peridynamics (correspondence) formulation with
bond-associated quadrature integration at the center of the bonds.

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
- `reprkernel::Symbol`: A kernel function used for the reproducing kernel approximation. See
    the constructor docs for more informations.
- `regfactor::Float64`: Regularization factor. See the constructor docs for more
    informations.

# Allowed material parameters
When using [`material!`](@ref) on a [`Body`](@ref) with `RKCMaterial`, then the following
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
`RKCMaterial`, the following fields are allowed:
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
struct RKCMaterial{CM,K,DM} <: AbstractRKCMaterial{CM,NoCorrection}
    kernel::K
    constitutive_model::CM
    dmgmodel::DM
    maxdmg::Float64
    reprkernel::Symbol
    regfactor::Float64
    function RKCMaterial(kernel::K, cm::CM, dmgmodel::DM, maxdmg::Real, reprkernel::Symbol,
                         regfactor::Real) where {CM,K,DM}
        return new{CM,K,DM}(kernel, cm, dmgmodel, maxdmg, reprkernel, regfactor)
    end
end

function RKCMaterial(; kernel::Function=cubic_b_spline_kernel,
                     model::AbstractConstitutiveModel=LinearElastic(),
                     dmgmodel::AbstractDamageModel=CriticalStretch(), maxdmg::Real=1.0,
                     reprkernel::Symbol=:C1, regfactor::Real=1e-13)
    get_q_dim(reprkernel) # check if the kernel is implemented
    if !(0 ≤ regfactor ≤ 1)
        msg = "Regularization factor must be in the range 0 ≤ regfactor ≤ 1\n"
        throw(ArgumentError(msg))
    end
    return RKCMaterial(kernel, model, dmgmodel, maxdmg, reprkernel, regfactor)
end

function Base.show(io::IO, @nospecialize(mat::AbstractRKCMaterial))
    print(io, typeof(mat))
    print(io, msg_fields_in_brackets(mat, (:maxdmg,)))
    return nothing
end

function log_material_property(::Val{:reprkernel}, mat; indentation)
    return msg_qty("reproducing kernel", mat.reprkernel; indentation)
end

function log_material_property(::Val{:regfactor}, mat; indentation)
    return msg_qty("regularization factor", mat.regfactor; indentation)
end

@params RKCMaterial StandardPointParameters

@storage RKCMaterial struct RKCStorage
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
    @pointfield update_gradients::Vector{Bool}
    @pointfield stress::Matrix{Float64}
    @pointfield von_mises_stress::Vector{Float64}
    @lthfield defgrad::Matrix{Float64}
    @lthfield weighted_volume::Vector{Float64}
    gradient_weight::Matrix{Float64}
    first_piola_kirchhoff::Matrix{Float64}
end

function init_field(::AbstractRKCMaterial, ::AbstractTimeSolver, system::BondSystem,
                    ::Val{:b_int})
    return zeros(3, get_n_points(system))
end

function init_field(::AbstractRKCMaterial, ::AbstractTimeSolver, system::BondSystem,
                    ::Val{:update_gradients})
    return ones(Bool, get_n_loc_points(system))
end

function init_field(::AbstractRKCMaterial, ::AbstractTimeSolver, system::BondSystem,
                    ::Val{:stress})
    return zeros(9, get_n_loc_points(system))
end

function init_field(::AbstractRKCMaterial, ::AbstractTimeSolver, system::BondSystem,
                    ::Val{:von_mises_stress})
    return zeros(get_n_loc_points(system))
end

function init_field(::AbstractRKCMaterial, ::AbstractTimeSolver, system::BondSystem,
                    ::Val{:defgrad})
    return zeros(9, get_n_points(system))
end

function init_field(::AbstractRKCMaterial, ::AbstractTimeSolver, system::BondSystem,
                    ::Val{:weighted_volume})
    return zeros(get_n_points(system))
end

function init_field(::AbstractRKCMaterial, ::AbstractTimeSolver, system::BondSystem,
                    ::Val{:gradient_weight})
    return zeros(3, get_n_bonds(system))
end

function init_field(::AbstractRKCMaterial, ::AbstractTimeSolver, system::BondSystem,
                    ::Val{:first_piola_kirchhoff})
    return zeros(9, get_n_bonds(system))
end

function initialize!(chunk::BodyChunk{<:BondSystem,<:AbstractRKCMaterial})
    (; system, mat, paramsetup, storage) = chunk
    storage.update_gradients .= true
    for i in each_point_idx(system)
        calc_weights_and_defgrad!(storage, system, mat, paramsetup, 0.0, 0.0, i)
    end
    return nothing
end

function calc_force_density!(dh::ThreadsBodyDataHandler{<:BondSystem,<:AbstractRKCMaterial},
                             t, Δt)
    (; chunks) = dh
    after_fields = rkc_lth_after_fields(chunks[1].mat)
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

rkc_lth_after_fields(::RKCMaterial) = (:defgrad, :weighted_volume)

function calc_force_density!(dh::MPIBodyDataHandler{<:BondSystem,<:AbstractRKCMaterial}, t,
                             Δt)
    (; chunk) = dh
    after_fields = rkc_lth_after_fields(chunk.mat)
    pre_fields = filter(x -> !in(x, after_fields), loc_to_halo_fields(chunk.storage))
    exchange_loc_to_halo!(dh, pre_fields)
    calc_weights_and_defgrad!(chunk, t, Δt)
    exchange_loc_to_halo!(dh, after_fields)
    calc_force_density!(chunk, t, Δt)
    exchange_halo_to_loc!(dh)
    return nothing
end

function calc_weights_and_defgrad!(chunk::BodyChunk{<:BondSystem,<:AbstractRKCMaterial}, t,
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
                      mat::AbstractRKCMaterial, dmgmodel::AbstractDamageModel,
                      paramsetup::AbstractParameterSetup, i)
    (; n_neighbors) = system
    (; n_active_bonds, damage, update_gradients) = storage
    old_damage = damage[i]
    new_damage = 1 - n_active_bonds[i] / n_neighbors[i]
    if new_damage > old_damage
        update_gradients[i] = true
    else
        update_gradients[i] = false
    end
    damage[i] = new_damage
    return nothing
end

function update_gradients(::AbstractRKCMaterial, storage::AbstractStorage, i)
    return storage.update_gradients[i]
end

function calc_weights_and_defgrad!(storage::AbstractStorage, system::AbstractBondSystem,
                                   mat::AbstractRKCMaterial,
                                   paramsetup::AbstractParameterSetup, t, Δt, i)
    params = get_params(paramsetup, i)
    if update_gradients(mat, storage, i)
        rkc_weights!(storage, system, mat, params, t, Δt, i)
    end
    rkc_defgrad!(storage, system, mat, params, t, Δt, i)
    return nothing
end

function rkc_weights!(storage::AbstractStorage, system::AbstractBondSystem,
                      mat::AbstractRKCMaterial, params::AbstractPointParameters, t, Δt, i)
    (; bonds, volume) = system
    (; bond_active, gradient_weight, weighted_volume, update_gradients) = storage
    (; reprkernel, regfactor) = mat
    (; δ) = params

    # get dimenion of the monomial vector and the gradient extraction matrix
    q_dim = get_q_dim(reprkernel)
    Q∇ᵀ = get_gradient_extraction_matrix(reprkernel)

    # calculate moment matrix M
    M = zero(SMatrix{q_dim,q_dim,Float64,q_dim*q_dim})
    wi = 0.0
    for bond_id in each_bond_idx(system, i)
        bond = bonds[bond_id]
        j = bond.neighbor
        ΔXij = get_vector_diff(system.position, i, j)
        Q = get_monomial_vector(reprkernel, ΔXij) ./ δ # scale by horizon
        ωij = kernel(system, bond_id) * bond_active[bond_id]
        temp = ωij * volume[j]
        M += temp * (Q * Q')
        wi += temp
    end
    weighted_volume[i] = wi

    # calculate inverse of moment matrix, must be a full rank matrix!
    Minv = invreg(M, regfactor * δ * δ * δ)

    # calculate gradient weights Φ
    for bond_id in each_bond_idx(system, i)
        bond = bonds[bond_id]
        j = bond.neighbor
        ΔXij = get_vector_diff(system.position, i, j)
        Q = get_monomial_vector(reprkernel, ΔXij) ./ δ # scale by horizon
        ωij = kernel(system, bond_id) * bond_active[bond_id]
        temp = ωij * volume[j] / δ # δ is needed due to scaling by horizon
        MinvQ = Minv * Q
        Φ = temp * (Q∇ᵀ * MinvQ)
        update_vector!(gradient_weight, bond_id, Φ)
    end

    # gradients are evaluated and do not need to be updated anymore
    update_gradients[i] = false

    return nothing
end

@inline get_q_dim(reprkernel::Symbol) = get_q_dim(Val(reprkernel))

function get_q_dim(::Val{reprkernel}) where {reprkernel}
    msg = "Reproducing kernel `$reprkernel` is not implemented!\n"
    return throw(ArgumentError(msg))
end

@inline function get_monomial_vector(rk::Symbol, ΔX::AbstractArray)
    return get_monomial_vector(Val(rk), ΔX)
end

function get_monomial_vector(::Val{reprkernel}, ΔX::AbstractArray) where {reprkernel}
    msg = "Reproducing kernel `$reprkernel` is not implemented!\n"
    return throw(ArgumentError(msg))
end

@inline function get_gradient_extraction_matrix(reprkernel::Symbol)
    return get_gradient_extraction_matrix(Val(reprkernel))
end

function get_gradient_extraction_matrix(::Val{reprkernel}) where {reprkernel}
    msg = "Reproducing kernel `$reprkernel` is not implemented!\n"
    return throw(ArgumentError(msg))
end

@inline get_q_dim(::Val{:C1}) = 3
@inline function get_monomial_vector(::Val{:C1}, ΔX::AbstractArray)
    x, y, z = ΔX[1], ΔX[2], ΔX[3]
    Q = SVector{3,eltype(ΔX)}(x, y, z)
    return Q
end
@inline function get_gradient_extraction_matrix(::Val{:C1})
    Q∇ᵀ = SMatrix{3,3,Int,9}(1, 0, 0, 0, 1, 0, 0, 0, 1)
    return Q∇ᵀ
end

@inline get_q_dim(::Val{:RK1}) = 4
@inline function get_monomial_vector(::Val{:RK1}, ΔX::AbstractArray)
    x, y, z = ΔX[1], ΔX[2], ΔX[3]
    Q = SVector{4,eltype(ΔX)}(1.0, x, y, z)
    return Q
end
@inline function get_gradient_extraction_matrix(::Val{:RK1})
    Q∇ᵀ = SMatrix{3,4,Int,12}(0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1)
    return Q∇ᵀ
end

@inline get_q_dim(::Val{:RK2}) = 7
@inline function get_monomial_vector(::Val{:RK2}, ΔX::AbstractArray)
    x, y, z = ΔX[1], ΔX[2], ΔX[3]
    Q = SVector{7,eltype(ΔX)}(1.0, x, y, z, x * x, y * y, z * z)
    return Q
end
@inline function get_gradient_extraction_matrix(::Val{:RK2})
    Q∇ᵀ = SMatrix{3,7,Int,21}(0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    return Q∇ᵀ
end

@inline get_q_dim(::Val{:PD2}) = 9
@inline function get_monomial_vector(::Val{:PD2}, ΔX::AbstractArray)
    x, y, z = ΔX[1], ΔX[2], ΔX[3]
    Q = SVector{9,eltype(ΔX)}(x, y, z, x * x, x * y, x * z, y * y, y * z, z * z)
    return Q
end
@inline function get_gradient_extraction_matrix(::Val{:PD2})
    Q∇ᵀ = SMatrix{3,9,Int,27}(1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                              0, 0, 0, 0, 0, 0, 0, 0, 0)
    return Q∇ᵀ
end

function rkc_defgrad!(storage::AbstractStorage, system::AbstractBondSystem,
                      mat::AbstractRKCMaterial, params::AbstractPointParameters, t, Δt, i)
    (; bonds) = system
    (; defgrad, gradient_weight) = storage

    F = zero(SMatrix{3,3,Float64,9})
    for bond_id in each_bond_idx(system, i)
        bond = bonds[bond_id]
        j = bond.neighbor
        Δxij = get_vector_diff(storage.position, i, j)
        Φij = get_vector(gradient_weight, bond_id)
        F += Δxij * Φij'
    end
    update_tensor!(defgrad, i, F)

    return nothing
end

function calc_force_density!(chunk::BodyChunk{<:BondSystem,<:AbstractRKCMaterial}, t, Δt)
    (; system, mat, paramsetup, storage) = chunk
    storage.b_int .= 0
    for i in each_point_idx(chunk)
        force_density_point!(storage, system, mat, paramsetup, t, Δt, i)
    end
    nancheck(chunk, t, Δt)
    return nothing
end

function force_density_point!(storage::AbstractStorage, system::AbstractSystem,
                              mat::AbstractRKCMaterial,
                              paramhandler::AbstractParameterHandler, t, Δt, i)
    params = get_params(paramhandler, i)
    force_density_point!(storage, system, mat, params, t, Δt, i)
    return nothing
end

function force_density_point!(storage::AbstractStorage, system::AbstractSystem,
                              mat::AbstractRKCMaterial, params::AbstractPointParameters, t,
                              Δt, i)
    too_much_damage!(storage, system, mat, i) && return nothing
    ∑P = rkc_stress_integral!(storage, system, mat, params, t, Δt, i)
    rkc_force_density!(storage, system, mat, params, ∑P, t, Δt, i)
    return nothing
end

function rkc_stress_integral!(storage::AbstractStorage, system::AbstractBondSystem,
                              mat::AbstractRKCMaterial, params::AbstractPointParameters, t,
                              Δt, i)
    (; bonds, volume) = system
    (; bond_active, defgrad, weighted_volume) = storage
    Fi = get_tensor(defgrad, i)
    wi = weighted_volume[i]
    ∑P = zero(SMatrix{3,3,Float64,9})
    for bond_id in each_bond_idx(system, i)
        if bond_active[bond_id]
            bond = bonds[bond_id]
            j, L = bond.neighbor, bond.length
            ΔXij = get_vector_diff(system.position, i, j)
            Δxij = get_vector_diff(storage.position, i, j)
            Fj = get_tensor(defgrad, j)
            Fij = bond_avg(Fi, Fj, ΔXij, Δxij, L)
            Pij = calc_first_piola_kirchhoff!(storage, mat, params, Fij, bond_id)
            Tempij = I - ΔXij * ΔXij' / (L * L)
            wj = weighted_volume[j]
            ϕ = 0.5 / wi + 0.5 / wj
            ω̃ij = kernel(system, bond_id) * ϕ * volume[j]
            ∑Pij = ω̃ij * (Pij * Tempij)
            ∑P += ∑Pij
        end
    end
    return ∑P
end

function rkc_force_density!(storage::AbstractStorage, system::AbstractBondSystem,
                            mat::AbstractRKCMaterial, params::AbstractPointParameters,
                            ∑P, t, Δt, i)
    (; bonds, volume) = system
    (; bond_active, gradient_weight, first_piola_kirchhoff, weighted_volume,
       b_int) = storage
    wi = weighted_volume[i]
    for bond_id in each_bond_idx(system, i)
        if bond_active[bond_id]
            bond = bonds[bond_id]
            j, L = bond.neighbor, bond.length
            ΔXij = get_vector_diff(system.position, i, j)
            Pij = get_tensor(first_piola_kirchhoff, bond_id)
            Φij = get_vector(gradient_weight, bond_id)
            ϕ = 1 / wi
            ω̃ij = kernel(system, bond_id) * ϕ
            tij = ω̃ij / (L * L) * (Pij * ΔXij) + ∑P * Φij / volume[j]
            update_add_vector!(b_int, i, tij * volume[j])
            update_add_vector!(b_int, j, -tij * volume[i])
        end
    end
    return nothing
end

function calc_first_piola_kirchhoff!(storage::RKCStorage, mat::RKCMaterial,
                                     params::StandardPointParameters, F, bond_id)
    P = first_piola_kirchhoff(mat.constitutive_model, storage, params, F)
    update_tensor!(storage.first_piola_kirchhoff, bond_id, P)
    return P
end

function too_much_damage!(storage::AbstractStorage, system::BondSystem,
                          mat::AbstractRKCMaterial, i)
    if storage.damage[i] > mat.maxdmg
        # kill all bonds of this point
        storage.bond_active[each_bond_idx(system, i)] .= false
        return true
    end
    return false
end

function bond_avg(Fi, Fj, ΔXij, Δxij, L)
    Favg = 0.5 * (Fi + Fj)
    Fcor = (Δxij - Favg * ΔXij) * (ΔXij' / (L * L))
    Fij = Favg + Fcor
    return Fij
end
