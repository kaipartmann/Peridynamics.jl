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
- `reprkernel::Symbol`: A kernel function used for the reproducing kernel approximation
    of the moment matrix. This kernel is used to calculate the moment matrix and the
    gradient weights, which are used to approximate the deformation gradient. \\
    (default: `:C1`) \\
    The following kernels can be used:
    - `:C1`: A kernel reproducing the corresponndence formulation of peridynamics
        with a first order accuracy. This is the default kernel.
- `regfactor::Float64`: A regularization factor used to stabilize the inversion of the
    moment matrix. The product of `regfactor * δ^3` is used for the regularization with
    [`invreg`](@ref). The regularization helps to ensure numerical stability by mitigating
    issues such as ill-conditioning or singularity in the matrix operations.\\
    (default: `1e-13`)

# Examples

```julia-repl
julia> mat = RKCMaterial()
RKCMaterial{LinearElastic, typeof(linear_kernel), CriticalStretch}(maxdmg=0.85)
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
                     dmgmodel::AbstractDamageModel=CriticalStretch(), maxdmg::Real=0.85,
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
    q_dim = get_q_dim(reprkernel)

    Q∇ᵀ = get_q_triangle(reprkernel)

    (; δ) = params

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

@inline get_q_triangle(reprkernel::Symbol) = get_q_triangle(Val(reprkernel))

function get_q_triangle(::Val{reprkernel}) where {reprkernel}
    msg = "Reproducing kernel `$reprkernel` is not implemented!\n"
    return throw(ArgumentError(msg))
end

@inline get_q_dim(::Val{:C1}) = 3
@inline function get_monomial_vector(::Val{:C1}, ΔX::AbstractArray)
    x, y, z = ΔX[1], ΔX[2], ΔX[3]
    Q = SVector{3,eltype(ΔX)}(x, y, z)
    return Q
end
@inline function get_q_triangle(::Val{:C1})
    Q∇ᵀ = SMatrix{3,3,Int,9}(1, 0, 0, 0, 1, 0, 0, 0, 1)
    return Q∇ᵀ
end

@inline get_q_dim(::Val{:RK1}) = 4
@inline function get_monomial_vector(::Val{:RK1}, ΔX::AbstractArray)
    x, y, z = ΔX[1], ΔX[2], ΔX[3]
    Q = SVector{4,eltype(ΔX)}(1.0, x, y, z)
    return Q
end
@inline function get_q_triangle(::Val{:RK1})
    Q∇ᵀ = SMatrix{3,4,Int,12}(0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1)
    return Q∇ᵀ
end

@inline get_q_dim(::Val{:RK2}) = 7
@inline function get_monomial_vector(::Val{:RK2}, ΔX::AbstractArray)
    x, y, z = ΔX[1], ΔX[2], ΔX[3]
    Q = SVector{7,eltype(ΔX)}(1.0, x, y, z, x * x, y * y, z * z)
    return Q
end
@inline function get_q_triangle(::Val{:RK2})
    Q∇ᵀ = SMatrix{3,7,Int,21}(0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    return Q∇ᵀ
end

@inline get_q_dim(::Val{:PD2}) = 9
@inline function get_monomial_vector(::Val{:PD2}, ΔX::AbstractArray)
    x, y, z = ΔX[1], ΔX[2], ΔX[3]
    Q = SVector{9,eltype(ΔX)}(x, y, z, x * x, x * y, x * z, y * y, y * z, z * z)
    return Q
end
@inline function get_q_triangle(::Val{:PD2})
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
    nancheck(chunk, t)
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
    too_much_damage!(storage, system, mat, Fi, i) && return zero(SMatrix{3,3,Float64,9})
    wi = weighted_volume[i]
    ∑P = zero(SMatrix{3,3,Float64,9})
    for bond_id in each_bond_idx(system, i)
        if bond_active[bond_id]
            bond = bonds[bond_id]
            j, L = bond.neighbor, bond.length
            ΔXij = get_vector_diff(system.position, i, j)
            Δxij = get_vector_diff(storage.position, i, j)
            Fj = get_tensor(defgrad, j)
            # bond averaging manually
            # Fb = 0.5 * (Fi + Fj)
            # ΔXijLL = ΔXij' / (L * L)
            # ΔFij = (Δxij - Fb * ΔXij) * ΔXijLL
            # Fij = Fb + ΔFij
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
                          mat::AbstractRKCMaterial, F, i)
    if storage.damage[i] > mat.maxdmg || containsnan(F)
        # kill all bonds of this point
        storage.bond_active[each_bond_idx(system, i)] .= false
        return true
    end
    return false
end

function bond_avg(Fi, Fj, ΔXij, Δxij, L)
    # Matrix writing notation:
    # A = [A11 A12 A13
    #      A21 A22 A23
    #      A31 A32 A33]

    #-- read the elements

    # # Fi
    # Fi11, Fi12, Fi13 = Fi[1], Fi[4], Fi[7]
    # Fi21, Fi22, Fi23 = Fi[2], Fi[5], Fi[8]
    # Fi31, Fi32, Fi33 = Fi[3], Fi[6], Fi[9]

    # # Fj
    # Fj11, Fj12, Fj13 = Fj[1], Fj[4], Fj[7]
    # Fj21, Fj22, Fj23 = Fj[2], Fj[5], Fj[8]
    # Fj31, Fj32, Fj33 = Fj[3], Fj[6], Fj[9]

    # # ΔXij
    # ΔX1, ΔX2, ΔX3 = ΔXij[1], ΔXij[2], ΔXij[3]

    # # Δxij
    # Δx1, Δx2, Δx3 = Δxij[1], Δxij[2], Δxij[3]

    # # check the norm
    # # @assert norm(ΔXij) ≈ L "|ΔXij| should be equal to L!"
    # Lsq = L * L

    # #-- calculate the average bond force
    # Favg11 = 0.5 * (Fi11 + Fj11)
    # Favg12 = 0.5 * (Fi12 + Fj12)
    # Favg13 = 0.5 * (Fi13 + Fj13)
    # Favg21 = 0.5 * (Fi21 + Fj21)
    # Favg22 = 0.5 * (Fi22 + Fj22)
    # Favg23 = 0.5 * (Fi23 + Fj23)
    # Favg31 = 0.5 * (Fi31 + Fj31)
    # Favg32 = 0.5 * (Fi32 + Fj32)
    # Favg33 = 0.5 * (Fi33 + Fj33)

    # # calculate the difference vector
    # ΔD1 = Δx1 - Favg11 * ΔX1 - Favg12 * ΔX2 - Favg13 * ΔX3
    # ΔD2 = Δx2 - Favg21 * ΔX1 - Favg22 * ΔX2 - Favg23 * ΔX3
    # ΔD3 = Δx3 - Favg31 * ΔX1 - Favg32 * ΔX2 - Favg33 * ΔX3

    # # calculate the correction term
    # Fcor11 = 1 / Lsq * (ΔD1 * ΔX1)
    # Fcor12 = 1 / Lsq * (ΔD1 * ΔX2)
    # Fcor13 = 1 / Lsq * (ΔD1 * ΔX3)
    # Fcor21 = 1 / Lsq * (ΔD2 * ΔX1)
    # Fcor22 = 1 / Lsq * (ΔD2 * ΔX2)
    # Fcor23 = 1 / Lsq * (ΔD2 * ΔX3)
    # Fcor31 = 1 / Lsq * (ΔD3 * ΔX1)
    # Fcor32 = 1 / Lsq * (ΔD3 * ΔX2)
    # Fcor33 = 1 / Lsq * (ΔD3 * ΔX3)

    # #-- assemble
    # F11 = Favg11 + Fcor11
    # F12 = Favg12 + Fcor12
    # F13 = Favg13 + Fcor13
    # F21 = Favg21 + Fcor21
    # F22 = Favg22 + Fcor22
    # F23 = Favg23 + Fcor23
    # F31 = Favg31 + Fcor31
    # F32 = Favg32 + Fcor32
    # F33 = Favg33 + Fcor33

    # # merge the element type of Fi and Fj
    # T = promote_type(eltype(Fi), eltype(Fj), eltype(ΔXij), eltype(Δxij))
    # Fij = SMatrix{3,3,T,9}(F11, F21, F31,
    #                        F12, F22, F32,
    #                        F13, F23, F33)

    #-- using StaticArrays for simplicity
    Favg = 0.5 * (Fi + Fj)
    Fcor = (Δxij - Favg * ΔXij) * (ΔXij' / (L * L))
    Fij = Favg + Fcor

    return Fij
end
