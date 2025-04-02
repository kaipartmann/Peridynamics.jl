"""
    NSCRMaterial(; kernel, model, zem, dmgmodel, maxdmg)

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
julia> mat = NSCRMaterial()
NSCRMaterial{LinearElastic, typeof(linear_kernel), CriticalStretch}(maxdmg=0.85)
```

---

```julia
NSCRMaterial{CM,K,DM}
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
When using [`material!`](@ref) on a [`Body`](@ref) with `NSCRMaterial`, then the following
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
`NSCRMaterial`, the following fields are allowed:
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
struct NSCRMaterial{CM,K,DM} <: AbstractNSCMaterial{CM,NoCorrection}
    kernel::K
    constitutive_model::CM
    dmgmodel::DM
    maxdmg::Float64
    function NSCRMaterial(kernel::K, cm::CM, dmgmodel::DM, maxdmg::Real) where {CM,K,DM}
        return new{CM,K,DM}(kernel, cm, dmgmodel, maxdmg)
    end
end

function NSCRMaterial(; kernel::Function=cubic_b_spline_kernel,
                      model::AbstractConstitutiveModel=LinearElastic(),
                      dmgmodel::AbstractDamageModel=CriticalStretch(), maxdmg::Real=0.85)
    if kernel !== cubic_b_spline_kernel
        msg = "The kernel $kernel is not recommended for use with NSCRMaterial."
        msg *= "Use the `cubic_b_spline_kernel` function instead.\n"
        throw(ArgumentError(msg))
    end
    return NSCRMaterial(kernel, model, dmgmodel, maxdmg)
end

function Base.show(io::IO, @nospecialize(mat::NSCRMaterial))
    print(io, typeof(mat))
    print(io, msg_fields_in_brackets(mat, (:maxdmg,)))
    return nothing
end

@params NSCRMaterial StandardPointParameters

@storage NSCRMaterial struct NSCRStorage
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
    @pointfield update_gradients::Vector{Bool}
    @pointfield stress::Matrix{Float64}
    @pointfield von_mises_stress::Vector{Float64}
    @lthfield defgrad::Matrix{Float64}
    @lthfield defgrad_dot::Matrix{Float64}
    @lthfield weighted_volume::Vector{Float64}
    gradient_weight::Matrix{Float64}
    first_piola_kirchhoff::Matrix{Float64}
    left_stretch::Matrix{Float64}
    rotation::Matrix{Float64}
    bond_stress::Matrix{Float64}
end

function init_field(::NSCRMaterial, ::AbstractTimeSolver, system::BondSystem,
                    ::Val{:velocity_half})
    return zeros(3, get_n_points(system))
end

function init_field(::NSCRMaterial, ::AbstractTimeSolver, system::BondSystem,
                    ::Val{:defgrad_dot})
    return zeros(9, get_n_points(system))
end

function init_field(::NSCRMaterial, ::AbstractTimeSolver, system::BondSystem,
                    ::Val{:left_stretch})
    V = zeros(9, get_n_bonds(system))
    V[[1, 5, 9], :] .= 1.0
    return V
end

function init_field(::NSCRMaterial, ::AbstractTimeSolver, system::BondSystem,
                    ::Val{:rotation})
    R = zeros(9, get_n_bonds(system))
    R[[1, 5, 9], :] .= 1.0
    return R
end

function init_field(::NSCRMaterial, ::AbstractTimeSolver, system::BondSystem,
                    ::Val{:bond_stress})
    return zeros(9, get_n_bonds(system))
end

nsc_lth_after_fields(::NSCRMaterial) = (:defgrad, :defgrad_dot, :weighted_volume)

function nsc_defgrad!(storage::NSCRStorage, system::BondSystem, mat::NSCRMaterial,
                      params::StandardPointParameters, t, Δt, i)
    (; bonds, volume) = system
    (; bond_active, defgrad, defgrad_dot, weighted_volume) = storage

    K = zero(SMatrix{3,3,Float64,9})
    _F = zero(SMatrix{3,3,Float64,9})
    _Ḟ = zero(SMatrix{3,3,Float64,9})
    wi = 0.0
    for bond_id in each_bond_idx(system, i)
        bond = bonds[bond_id]
        j = bond.neighbor
        ΔXij = get_vector_diff(system.position, i, j)
        Δxij = get_vector_diff(storage.position, i, j)
        Δvij = get_vector_diff(storage.velocity_half, i, j)
        ωij = kernel(system, bond_id) * bond_active[bond_id]
        temp = ωij * volume[j]
        K += temp * (ΔXij * ΔXij')
        _F += temp * (Δxij * ΔXij')
        _Ḟ += temp * (Δvij * ΔXij')
        wi += temp
    end
    Kinv = inv(K)
    F = _F * Kinv
    Ḟ = _Ḟ * Kinv
    update_tensor!(defgrad, i, F)
    update_tensor!(defgrad_dot, i, Ḟ)
    weighted_volume[i] = wi

    return Kinv
end

function nsc_stress_integral!(storage::NSCRStorage, system::BondSystem, mat::NSCRMaterial,
                              params::StandardPointParameters, t, Δt, i)
    (; bonds) = system
    (; bond_active, defgrad, defgrad_dot, weighted_volume) = storage
    Fi = get_tensor(defgrad, i)
    Ḟi = get_tensor(defgrad_dot, i)
    too_much_damage!(storage, system, mat, Fi, i) && return nothing
    wi = weighted_volume[i]
    ∑P = SMatrix{3,3,Float64,9}(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    for bond_id in each_bond_idx(system, i)
        if bond_active[bond_id]
            bond = bonds[bond_id]
            j, L = bond.neighbor, bond.length
            ΔXij = get_vector_diff(system.position, i, j)
            Δxij = get_vector_diff(storage.position, i, j)
            Δvij = get_vector_diff(storage.velocity_half, i, j)
            Fj = get_tensor(defgrad, j)
            Ḟj = get_tensor(defgrad_dot, j)
            Fb = 0.5 * (Fi + Fj)
            Ḟb = 0.5 * (Ḟi + Ḟj)
            ΔXijLL = ΔXij' / (L * L)
            ΔFij = (Δxij - Fb * ΔXij) * ΔXijLL
            ΔḞij = (Δvij - Ḟb * ΔXij) * ΔXijLL
            Fij = Fb + ΔFij
            Ḟij = Ḟb + ΔḞij
            Pij = calc_first_piola_kirchhoff!(storage, mat, params, Fij, Ḟij, Δt, bond_id)
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

function calc_first_piola_kirchhoff!(storage::NSCRStorage, mat::NSCRMaterial,
                                     params::StandardPointParameters, F::SMatrix{3,3,FT,9},
                                     Ḟ::SMatrix{3,3,FT,9}, Δt, bond_id) where {FT}
    D = init_stress_rotation!(storage, F, Ḟ, Δt, bond_id)
    iszero(D) && return zero(SMatrix{3,3,FT,9})
    Δε = D * Δt
    Δθ = tr(Δε)
    Δεᵈᵉᵛ = Δε - Δθ / 3 * I
    σ = get_tensor(storage.bond_stress, bond_id)
    σₙ₊₁ = σ + 2 * params.G * Δεᵈᵉᵛ + params.K * Δθ * I
    update_tensor!(storage.bond_stress, bond_id, σₙ₊₁)
    T = rotate_stress(storage, σₙ₊₁, bond_id)
    P = first_piola_kirchhoff(T, F)
    update_tensor!(storage.first_piola_kirchhoff, bond_id, P)
    return P
end
