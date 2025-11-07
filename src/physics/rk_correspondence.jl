"""
    RKCMaterial(; kernel, model, dmgmodel, monomial, regfactor)

A material type used to assign the material of a [`Body`](@ref) with a reproducing kernel
peridynamics (correspondence) formulation with bond-associated quadrature integration at
the center of the bonds.

# Keywords
- `kernel::Function`: Kernel function used for weighting the interactions between points. \\
    (default: `cubic_b_spline_kernel_norm`) \\
    The following kernels can be used:
    - [`const_one_kernel`](@ref)
    - [`linear_kernel`](@ref)
    - [`cubic_b_spline_kernel_norm`](@ref)
    - [`cubic_b_spline_kernel`](@ref)
- `model::AbstractConstitutiveModel`: Constitutive model defining the material behavior. \\
    (default: `SaintVenantKirchhoff()`) \\
    The following models can be used:
    - [`SaintVenantKirchhoff`](@ref)
    - [`LinearElastic`](@ref)
    - [`NeoHooke`](@ref)
    - [`NeoHookePenalty`](@ref)
- `dmgmodel::AbstractDamageModel`: Damage model defining the damage behavior. \\
    (default: [`CriticalStretch`](@ref))
- `monomial::Symbol`: The monomial vector used for the reproducing kernel approximation
    of the moment matrix. This kernel is used to calculate the moment matrix and the
    gradient weights, which are used to approximate the deformation gradient. \\
    (default: `:C1`) \\
    The following kernels can be used:
    - `:C1`: Linear monomial basis vector [x, y, z] with first-order accuracy, equivalent to
        the standard correspondence formulation. This is the default kernel.
    - `:RK1`: First-order monomial basis vector [1, x, y, z] with constant term
        for enhanced stability in uniform deformation fields.
    - `:RK2`: Second-order monomial basis vector [1, x, y, z, x², y², z²] with
        diagonal quadratic terms for improved accuracy in curved deformation fields.
    - `:PD2`: Second-order monomial basis vector [x, y, z, x², xy, xz, y², yz, z²] with
        full quadratic terms but without constant term.
- `lambda::Float64`: Tikhonov regularization parameter used to stabilize the inversion of
    the moment matrix. This parameter controls the strength of the Tikhonov regularization
    applied during the pseudo-inverse computation with [`invreg`](@ref). Should be a
    non-negative value between 0 and 1, where larger values increase regularization
    strength. See the [`invreg`](@ref) documentation for parameter selection guidelines.\\
    (default: `0`)
- `beta::Float64`: SVD truncation parameter used to stabilize the inversion of the moment
    matrix. This parameter defines the threshold (as a fraction of the largest singular
    value) below which singular values are set to zero during the pseudo-inverse
    computation with [`invreg`](@ref). Should be a positive value between 0 and 1.
    See the [`invreg`](@ref) documentation for parameter selection guidelines.\\
    (default: `sqrt(eps())`)

# Examples

```julia-repl
julia> mat = RKCMaterial()
RKCMaterial{SaintVenantKirchhoff, typeof(linear_kernel), CriticalStretch}()
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
- `monomial::Symbol`: The monomial vector used for the reproducing kernel approximation. See
    the constructor docs for more informations.
- `lambda::Float64`: Tikhonov regularization parameter. See the constructor docs for more
    informations.
- `beta::Float64`: SVD truncation parameter. See the constructor docs for more
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
- `cauchy_stress::Matrix{Float64}`: Cauchy stress tensor of each point
- `von_mises_stress::Vector{Float64}`: Von Mises stress of each point
- `hydrostatic_stress::Vector{Float64}`: Hydrostatic stress of each point.
- `strain_energy_density::Vector{Float64}`: Strain energy density of each point.
"""
struct RKCMaterial{CM,K,DM} <: AbstractRKCMaterial{CM,NoCorrection}
    kernel::K
    constitutive_model::CM
    dmgmodel::DM
    monomial::Symbol
    lambda::Float64
    beta::Float64
    function RKCMaterial(kernel::K, cm::CM, dmgmodel::DM, monomial::Symbol,
                         lambda::Float64, beta::Float64) where {CM,K,DM}
        return new{CM,K,DM}(kernel, cm, dmgmodel, monomial, lambda, beta)
    end
end

function RKCMaterial(; kernel::Function=cubic_b_spline_kernel_norm,
                     model::AbstractConstitutiveModel=SaintVenantKirchhoff(),
                     dmgmodel::AbstractDamageModel=CriticalStretch(),
                     monomial::Symbol=:C1, lambda::Real=0, beta::Real=sqrt(eps()))
    get_q_dim(monomial) # check if the kernel is implemented
    if !(0 ≤ lambda ≤ 1)
        msg = "Regularization factor must be in the range 0 ≤ lambda ≤ 1\n"
        throw(ArgumentError(msg))
    end
    if !(0 ≤ beta ≤ 1)
        msg = "Regularization factor must be in the range 0 ≤ beta ≤ 1\n"
        throw(ArgumentError(msg))
    end
    return RKCMaterial(kernel, model, dmgmodel, monomial, lambda, beta)
end

function Base.show(io::IO, @nospecialize(mat::AbstractRKCMaterial))
    print(io, typeof(mat))
    print(io, msg_fields_in_brackets(mat, ()))
    return nothing
end

function log_material_property(::Val{:monomial}, mat; indentation)
    return msg_qty("monomial type", mat.monomial; indentation)
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
    @pointfield n_active_bonds::Vector{Int}
    @pointfield update_gradients::Vector{Bool}
    @pointfield cauchy_stress::Matrix{Float64}
    @pointfield von_mises_stress::Vector{Float64}
    @pointfield strain_energy_density::Vector{Float64}
    @lthfield defgrad::Matrix{Float64}
    @lthfield weighted_volume::Vector{Float64}
    bond_active::Vector{Bool}
    gradient_weight::Matrix{Float64}
    bond_first_piola_kirchhoff::Matrix{Float64}
    residual::Vector{Float64}
    jacobian::Matrix{Float64}
    displacement_copy::Matrix{Float64}
    b_int_copy::Matrix{Float64}
    temp_force_a::Vector{Float64}
    temp_force_b::Vector{Float64}
    Δu::Vector{Float64}
    affected_points::Vector{Vector{Int}}
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
                    ::Val{:cauchy_stress})
    return zeros(9, get_n_loc_points(system))
end

function init_field(::AbstractRKCMaterial, ::AbstractTimeSolver, system::BondSystem,
                    ::Val{:von_mises_stress})
    return zeros(get_n_loc_points(system))
end

function init_field(::AbstractRKCMaterial, ::AbstractTimeSolver, system::BondSystem,
                    ::Val{:strain_energy_density})
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
                    ::Val{:bond_first_piola_kirchhoff})
    return zeros(9, get_n_bonds(system))
end

# For RKCMaterial, we need two levels of neighbors because:
# - Force at point i uses defgrad of i and all neighbors j
# - Force at neighbor j uses defgrad of j and all neighbors k of j
# So when perturbing i, we need defgrad updates for i, all j, and all k
# WARNING: This only works because for NR there is only one chunk and the system knows
# all points. For multithreading, additional layers of halo points would be needed.
function get_affected_points_rkc(system::AbstractBondSystem, i)
    affected_points = get_affected_points(system, i)
    second_level = Vector{Int}()
    for i in affected_points
        for bond_id in each_bond_idx(system, i)
            bond = system.bonds[bond_id]
            j = bond.neighbor
            if j ∉ affected_points && j ∉ second_level
                push!(second_level, j)
            end
        end
    end
    append!(affected_points, second_level)
    sort!(affected_points)
    return affected_points
end

# Override affected_points initialization for RKCMaterial
function init_field(::AbstractRKCMaterial, ::NewtonRaphson, system::BondSystem,
                    ::Val{:affected_points})
    affected_points = Vector{Vector{Int}}(undef, get_n_loc_points(system))
    for i in each_point_idx(system)
        affected_points[i] = get_affected_points_rkc(system, i)
    end
    return affected_points
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

# Newton-Raphson specific method: calculate forces for perturbed positions
# During Jacobian assembly, we only update deformation gradients, NOT weights or damage.
# Gradient weights and bond failure state are frozen at the current Newton-Raphson iterate.
function calc_force_density!(storage::AbstractStorage, system::AbstractBondSystem,
                             mat::AbstractRKCMaterial, paramsetup::AbstractParameterSetup,
                             idxs::AbstractVector{Int}, t, Δt)
    @inbounds storage.b_int[:, idxs] .= 0.0

    # Update deformation gradient for all affected points
    # Note: idxs contains point i and all its neighbors, so all required defgrads are updated
    for i in idxs
        rkc_defgrad!(storage, system, mat, get_params(paramsetup, i), t, Δt, i)
    end

    # Calculate force density - only for the original perturbed point and its direct neighbors
    # Forces at more distant points don't change due to locality
    for i in idxs
        force_density_point!(storage, system, mat, paramsetup, t, Δt, i)
    end

    nancheck(storage, t, Δt)
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
    (; monomial, lambda, beta) = mat
    (; δ) = params

    # get dimenion of the monomial vector and the gradient extraction matrix
    q_dim = get_q_dim(monomial)
    Q∇ᵀ = get_gradient_extraction_matrix(monomial)

    # calculate moment matrix M
    M = zero(SMatrix{q_dim,q_dim,Float64,q_dim*q_dim})
    wi = 0.0
    for bond_id in each_bond_idx(system, i)
        bond = bonds[bond_id]
        j = bond.neighbor
        ΔXij = get_vector_diff(system.position, i, j)
        Q = get_monomial_vector(monomial, ΔXij ./ δ) # normalize by δ
        ωij = kernel(system, bond_id) * bond_active[bond_id]
        temp = ωij * volume[j]
        M += temp * (Q * Q')
        wi += temp
    end
    weighted_volume[i] = wi

    # calculate regularized inverse of M
    Minv = invreg(M, lambda, beta)

    # calculate gradient weights Φ
    for bond_id in each_bond_idx(system, i)
        bond = bonds[bond_id]
        j = bond.neighbor
        ΔXij = get_vector_diff(system.position, i, j)
        Q = get_monomial_vector(monomial, ΔXij ./ δ) # normalize by δ
        ωij = kernel(system, bond_id) * bond_active[bond_id]
        temp = ωij / δ * volume[j] # note the division by δ here, due to normalization of Q
        MinvQ = Minv * Q
        Φ = temp * (Q∇ᵀ * MinvQ)
        update_vector!(gradient_weight, bond_id, Φ)
    end

    # gradients are evaluated and do not need to be updated anymore
    update_gradients[i] = false

    return nothing
end

@inline get_q_dim(monomial::Symbol) = get_q_dim(Val(monomial))

function get_q_dim(::Val{monomial}) where {monomial}
    msg = "Reproducing kernel `$monomial` is not implemented!\n"
    return throw(ArgumentError(msg))
end

@inline function get_monomial_vector(rk::Symbol, ΔX::AbstractArray)
    return get_monomial_vector(Val(rk), ΔX)
end

function get_monomial_vector(::Val{monomial}, ΔX::AbstractArray) where {monomial}
    msg = "Reproducing kernel `$monomial` is not implemented!\n"
    return throw(ArgumentError(msg))
end

@inline function get_gradient_extraction_matrix(monomial::Symbol)
    return get_gradient_extraction_matrix(Val(monomial))
end

function get_gradient_extraction_matrix(::Val{monomial}) where {monomial}
    msg = "Reproducing kernel `$monomial` is not implemented!\n"
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

    F = SMatrix{3,3,Float64,9}(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)
    for bond_id in each_bond_idx(system, i)
        bond = bonds[bond_id]
        j = bond.neighbor
        ΔXij = get_vector_diff(system.position, i, j)
        Δxij = get_vector_diff(storage.position, i, j)
        Δuij = Δxij - ΔXij
        Φij = get_vector(gradient_weight, bond_id)
        F += Δuij * Φij' # maybe calculating the displacement gradient is more stable?
    end
    update_tensor!(defgrad, i, F)

    return nothing
end

function calc_force_density!(storage::AbstractStorage, system::AbstractBondSystem,
                             mat::AbstractRKCMaterial, paramsetup::AbstractParameterSetup,
                             t, Δt)
    storage.b_int .= 0.0
    for i in each_point_idx(system)
        force_density_point!(storage, system, mat, paramsetup, t, Δt, i)
    end
    nancheck(storage, t, Δt)
    return nothing
end

function force_density_point!(storage::AbstractStorage, system::AbstractSystem,
                              mat::AbstractRKCMaterial, params::AbstractPointParameters, t,
                              Δt, i)
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
    (; bond_active, gradient_weight, bond_first_piola_kirchhoff, weighted_volume,
       b_int) = storage
    wi = weighted_volume[i]
    for bond_id in each_bond_idx(system, i)
        if bond_active[bond_id]
            bond = bonds[bond_id]
            j, L = bond.neighbor, bond.length
            ΔXij = get_vector_diff(system.position, i, j)
            Pij = get_tensor(bond_first_piola_kirchhoff, bond_id)
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
    update_tensor!(storage.bond_first_piola_kirchhoff, bond_id, P)
    return P
end

function bond_avg(Fi, Fj, ΔXij, Δxij, L)
    Favg = 0.5 * (Fi + Fj)
    Fcor = (Δxij - Favg * ΔXij) * (ΔXij' / (L * L))
    Fij = Favg + Fcor
    return Fij
end

function cauchy_stress_point!(storage::AbstractStorage, system::BondSystem,
                              ::AbstractRKCMaterial, ::StandardPointParameters, i)
    (; bonds) = system
    (; bond_active, defgrad, bond_first_piola_kirchhoff, n_active_bonds) = storage
    Fi = get_tensor(defgrad, i)
    σi = zero(SMatrix{3,3,Float64,9})
    for bond_id in each_bond_idx(system, i)
        if bond_active[bond_id]
            bond = bonds[bond_id]
            j, L = bond.neighbor, bond.length
            ΔXij = get_vector_diff(system.position, i, j)
            Δxij = get_vector_diff(storage.position, i, j)
            Fj = get_tensor(defgrad, j)
            Fij = bond_avg(Fi, Fj, ΔXij, Δxij, L)
            Pij = get_tensor(bond_first_piola_kirchhoff, bond_id)
            σij = cauchy_stress(Pij, Fij)
            σi += σij
        end
    end
    update_tensor!(storage.cauchy_stress, i, σi ./ n_active_bonds[i])
    return nothing
end

# calculate the Cauchy stress tensor just when exporting
function export_field(::Val{:cauchy_stress}, mat::AbstractRKCMaterial, system::BondSystem,
                      storage::AbstractStorage, paramsetup::AbstractParameterSetup, t)
    for i in each_point_idx(system)
        params = get_params(paramsetup, i)
        cauchy_stress_point!(storage, system, mat, params, i)
    end
    return storage.cauchy_stress
end

# calculate the von Mises stress from the Cauchy stress tensor just when exporting
function export_field(::Val{:von_mises_stress}, mat::AbstractRKCMaterial,
                      system::BondSystem, storage::AbstractStorage,
                      paramsetup::AbstractParameterSetup, t)
    for i in each_point_idx(system)
        params = get_params(paramsetup, i)
        cauchy_stress_point!(storage, system, mat, params, i)
        σ = get_tensor(storage.cauchy_stress, i)
        storage.von_mises_stress[i] = von_mises_stress(σ)
    end
    return storage.von_mises_stress
end

# calculate the von hydrostatic stress from the Cauchy stress tensor just when exporting,
# use the `von_mises_stress` field to store the hydrostatic stress
function export_field(::Val{:hydrostatic_stress}, mat::AbstractRKCMaterial,
                      system::BondSystem, storage::AbstractStorage,
                      paramsetup::AbstractParameterSetup, t)
    for i in each_point_idx(system)
        params = get_params(paramsetup, i)
        cauchy_stress_point!(storage, system, mat, params, i)
        σ = get_tensor(storage.cauchy_stress, i)
        storage.von_mises_stress[i] = 1/3 * (σ[1,1] + σ[2,2] + σ[3,3])
    end
    return storage.von_mises_stress
end
custom_field(::Type{RKCStorage}, ::Val{:hydrostatic_stress}) = true

# calculate the strain energy density from the deformation gradient just when exporting
function export_field(::Val{:strain_energy_density}, mat::AbstractRKCMaterial,
                      system::BondSystem, storage::AbstractStorage,
                      paramsetup::AbstractParameterSetup, t)
    model = mat.constitutive_model
    for i in each_point_idx(system)
        params = get_params(paramsetup, i)
        F = get_tensor(storage.defgrad, i)
        storage.strain_energy_density[i] = strain_energy_density(model, storage, params, F)
    end
    return storage.strain_energy_density
end
