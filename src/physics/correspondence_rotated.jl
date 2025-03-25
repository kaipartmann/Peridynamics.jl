"""
    CRMaterial(; kernel, model, zem, dmgmodel, maxdmg)

A material type used to assign the material of a [`Body`](@ref) with the local continuum
consistent (correspondence) formulation of non-ordinary state-based peridynamics.

# Keywords
- `kernel::Function`: Kernel function used for weighting the interactions between points. \\
    (default: `linear_kernel`) \\
    The following kernels can be used:
    - [`linear_kernel`](@ref)
    - [`cubic_b_spline_kernel`](@ref)
- `model::AbstractConstitutiveModel`: Constitutive model defining the material behavior. \\
    (default: `LinearElastic()`) \\
    Only the following model can be used:
    - [`LinearElastic`](@ref)
- `zem::AbstractZEMStabilization`: Zero-energy mode stabilization. The
    stabilization algorithm of Silling (2017) is used as default. \\
    (default: [`ZEMSilling`](@ref))
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
julia> mat = CRMaterial()
CRMaterial{LinearElastic, ZEMSilling, typeof(linear_kernel), CriticalStretch}(maxdmg=0.85)
```

---

```julia
CRMaterial{CM,ZEM,K,DM}
```

Material type for the local continuum consistent (correspondence) formulation of
non-ordinary state-based peridynamics.

# Type Parameters
- `CM`: A constitutive model type. See the constructor docs for more informations.
- `ZEM`: A zero-energy mode stabilization type. See the constructor docs for more
         informations.
- `K`: A kernel function type. See the constructor docs for more informations.
- `DM`: A damage model type. See the constructor docs for more informations.

# Fields
- `kernel::Function`: Kernel function used for weighting the interactions between points.
    See the constructor docs for more informations.
- `model::AbstractConstitutiveModel`: Constitutive model defining the material behavior. See
    the constructor docs for more informations.
- `zem::AbstractZEMStabilization`: Zero-energy mode stabilization. See the constructor docs
    for more informations.
- `dmgmodel::AbstractDamageModel`: Damage model defining the damage behavior. See the
    constructor docs for more informations.
- `maxdmg::Float64`: Maximum value of damage a point is allowed to obtain. See the
    constructor docs for more informations.

# Allowed material parameters
When using [`material!`](@ref) on a [`Body`](@ref) with `CRMaterial`, then the following
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
`CRMaterial`, the following fields are allowed:
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
struct CRMaterial{CM,ZEM,K,DM} <: AbstractCorrespondenceMaterial{CM,ZEM}
    kernel::K
    constitutive_model::CM
    zem_stabilization::ZEM
    dmgmodel::DM
    maxdmg::Float64
    function CRMaterial(kernel::K, cm::CM, zem::ZEM, dmgmodel::DM,
                       maxdmg::Real) where {CM,ZEM,K,DM}
        return new{CM,ZEM,K,DM}(kernel, cm, zem, dmgmodel, maxdmg)
    end
end

function CRMaterial(; kernel::Function=linear_kernel,
                    model::AbstractConstitutiveModel=LinearElastic(),
                    zem::AbstractZEMStabilization=ZEMSilling(),
                    dmgmodel::AbstractDamageModel=CriticalStretch(), maxdmg::Real=0.85)
    if !(model isa LinearElastic)
        msg = "only `LinearElastic` is supported as model for `CRMaterial`!\n"
        throw(ArgumentError(msg))
    end
    return CRMaterial(kernel, model, zem, dmgmodel, maxdmg)
end

function Base.show(io::IO, @nospecialize(mat::CRMaterial))
    print(io, typeof(mat))
    print(io, msg_fields_in_brackets(mat, (:maxdmg,)))
    return nothing
end

@params CRMaterial StandardPointParameters

@storage CRMaterial struct CRStorage
    @lthfield position::Matrix{Float64}
    @pointfield displacement::Matrix{Float64}
    @pointfield velocity::Matrix{Float64}
    @lthfield velocity_half::Matrix{Float64}
    @pointfield velocity_half_old::Matrix{Float64}
    @pointfield acceleration::Matrix{Float64}
    @htlfield b_int::Matrix{Float64}
    @pointfield b_int_old::Matrix{Float64}
    @pointfield b_ext::Matrix{Float64}
    @pointfield density_matrix::Matrix{Float64}
    @pointfield damage::Vector{Float64}
    bond_active::Vector{Bool}
    @pointfield n_active_bonds::Vector{Int}
    @pointfield stress::Matrix{Float64}
    @pointfield von_mises_stress::Vector{Float64}
    @pointfield left_stretch::Matrix{Float64}
    @pointfield rotation::Matrix{Float64}
end

function init_field(::CRMaterial, ::AbstractTimeSolver, system::BondSystem,
                    ::Val{:velocity_half})
    return zeros(3, get_n_points(system))
end

function init_field(::CRMaterial, ::AbstractTimeSolver, system::BondSystem, ::Val{:b_int})
    return zeros(3, get_n_points(system))
end

function init_field(::CRMaterial, ::AbstractTimeSolver, system::BondSystem, ::Val{:stress})
    return zeros(9, get_n_loc_points(system))
end

function init_field(::CRMaterial, ::AbstractTimeSolver, system::BondSystem,
                    ::Val{:von_mises_stress})
    return zeros(get_n_loc_points(system))
end

function init_field(::CRMaterial, ::AbstractTimeSolver, system::BondSystem,
                    ::Val{:left_stretch})
    V = zeros(9, get_n_loc_points(system))
    V[[1, 5, 9], :] .= 1.0
    return V
end

function init_field(::CRMaterial, ::AbstractTimeSolver, system::BondSystem,
                    ::Val{:rotation})
    R = zeros(9, get_n_loc_points(system))
    R[[1, 5, 9], :] .= 1.0
    return R
end

function calc_deformation_gradient(storage::CRStorage, system::BondSystem, ::CRMaterial,
                                   ::StandardPointParameters, i)
    (; bonds, volume) = system
    (; bond_active) = storage
    K = zero(SMatrix{3,3,Float64,9})
    _F = zero(SMatrix{3,3,Float64,9})
    _Ḟ = zero(SMatrix{3,3,Float64,9})
    ω0 = 0.0
    for bond_id in each_bond_idx(system, i)
        bond = bonds[bond_id]
        j = bond.neighbor
        ΔXij = get_vector_diff(system.position, i, j)
        Δxij = get_vector_diff(storage.position, i, j)
        Δvij = get_vector_diff(storage.velocity_half, i, j)
        ωij = kernel(system, bond_id) * bond_active[bond_id]
        ω0 += ωij
        temp = ωij * volume[j]
        ΔXijt = ΔXij'
        K += temp * (ΔXij * ΔXijt)
        _F += temp * (Δxij * ΔXijt)
        _Ḟ += temp * (Δvij * ΔXijt)
    end
    Kinv = inv(K)
    F = _F * Kinv
    Ḟ = _Ḟ * Kinv
    return (; F, Ḟ, Kinv, ω0)
end

function calc_first_piola_kirchhoff!(storage::CRStorage, mat::CRMaterial,
                                     params::StandardPointParameters, defgrad_res, Δt, i)

    ### New implementation
    (; F, Ḟ, Kinv) = defgrad_res
    D = init_stress_rotation!(storage, F, Ḟ, Δt, i)
    if iszero(D)
        storage.von_mises_stress[i] = 0.0
        return zero(SMatrix{3,3,Float64,9})
    end
    Δε = D * Δt
    Δθ = tr(Δε)
    Δεᵈᵉᵛ = Δε - Δθ / 3 * I
    σ = get_tensor(storage.stress, i)
    σₙ₊₁ = σ + 2 * params.G * Δεᵈᵉᵛ + params.K * Δθ * I
    update_tensor!(storage.stress, i, σₙ₊₁)
    T = rotate_stress(storage, σₙ₊₁, i)
    storage.von_mises_stress[i] = von_mises_stress(T)
    P = first_piola_kirchhoff(T, F)
    PKinv = P * Kinv
    return PKinv



    # ### SaintVenantKirchhoff
    # (; F, Ḟ, Kinv) = defgrad_res

    # # Compute Green-Lagrange strain tensor
    # E = 0.5 .* (F' * F - I)

    # # Compute strain rate (Lie derivative of Green-Lagrange strain)
    # Ė = symmetrize(F' * Ḟ)  # Time derivative of strain tensor

    # # Retrieve previous stress tensor (Second Piola-Kirchhoff stress)
    # S_prev = get_tensor(storage.stress, i)

    # # Compute stress rate (Lie derivative of the Second Piola-Kirchhoff stress)
    # S_dot = params.λ * tr(Ė) * I + 2 * params.μ * Ė

    # # Incremental stress update
    # S_new = S_prev + S_dot * Δt

    # # Store updated stress tensor
    # update_tensor!(storage.stress, i, S_new)

    # # Convert to First Piola-Kirchhoff stress
    # P = F * S_new

    # # Apply objectivity enforcement using Flanagan & Taylor's rotation algorithm
    # D = init_stress_rotation!(storage, F, Ḟ, Δt, i)

    # if iszero(D)
    #     storage.von_mises_stress[i] = 0.0
    #     return zero(SMatrix{3,3,Float64,9})
    # end

    # # Rotate stress tensor
    # P_rotated = rotate_stress(storage, P, i)

    # # Compute von Mises stress (for monitoring)
    # storage.von_mises_stress[i] = von_mises_stress(P_rotated)

    # # Convert to final stress measure
    # PKinv = P_rotated * Kinv




    # ### Neo-NeoHooke
    # (; F, Ḟ, Kinv) = defgrad_res

    # # Compute deformation tensors
    # J = det(F)                # Jacobian determinant

    # # Prevent issues with inversion if J is too small
    # if J < 1e-8
    #     @warn "Jacobian determinant J = $J is too small, setting J = 1e-8"
    #     J = 1e-8  # Prevent singularities
    # end

    # C = F' * F                # Right Cauchy-Green strain tensor
    # Cinv = safe_inverse(C)    # Use a stabilized inverse function
    # I1 = tr(C)                # First invariant of C
    # μ, κ = params.μ, params.K # Shear and bulk modulus

    # # Compute strain rate (Lie derivative of C)
    # Ċ = symmetrize(F' * Ḟ)

    # # Ensure the trace computation is numerically stable
    # tr_Ċ = tr(Ċ)
    # if tr_Ċ < 0 && abs(tr_Ċ) < 1e-12  # Numerical noise handling
    #     tr_Ċ = 0.0
    # end

    # # Compute stress rate (Lie derivative of Second Piola-Kirchhoff stress)
    # S_dot = μ * J^(-2/3) * (I - (1/3) * I1 * Cinv) * Ċ - (2/3) * μ * J^(-2/3) * tr_Ċ * Cinv +
    #         κ * (J - 1) * J * Cinv * tr_Ċ

    # # Retrieve previous stress tensor (Second Piola-Kirchhoff stress)
    # S_prev = get_tensor(storage.stress, i)

    # # Incremental stress update
    # S_new = S_prev + S_dot * Δt

    # # Store updated stress tensor
    # update_tensor!(storage.stress, i, S_new)

    # # Convert to First Piola-Kirchhoff stress
    # P = F * S_new

    # # Apply objectivity enforcement using Flanagan & Taylor's rotation algorithm
    # D = init_stress_rotation!(storage, F, Ḟ, Δt, i)

    # if iszero(D)
    #     storage.von_mises_stress[i] = 0.0
    #     return zero(SMatrix{3,3,Float64,9})
    # end

    # # Rotate stress tensor
    # P_rotated = rotate_stress(storage, P, i)

    # # Compute von Mises stress (for monitoring)
    # storage.von_mises_stress[i] = von_mises_stress(P_rotated)

    # # Convert to final stress measure
    # PKinv = P_rotated * Kinv

    # return PKinv
end

# # Safe inverse function to prevent singularity issues
# @inline function safe_inverse(A::SMatrix{3,3,T,9}) where T
#     try
#         return inv(A)
#     catch e
#         @warn "Matrix inversion failed, returning identity matrix as fallback"
#         return I
#     end
# end

# Helper function to ensure symmetry
# @inline function symmetrize(A::SMatrix{3,3,T,9}) where T
#     return 0.5 * (A + A')
# end
