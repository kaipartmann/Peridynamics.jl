"""
    CRMaterial(; kernel, model, zem, dmgmodel, maxdmg)

A material type used to assign the material of a [`Body`](@ref) with the local continuum
consistent (correspondence) formulation of non-ordinary state-based peridynamics.

# Keywords
- `kernel::Function`: Kernel function used for weighting the interactions between points. \\
    (default: `linear_kernel`) \\
    The following kernels can be used:
    - [`const_one_kernel`](@ref)
    - [`linear_kernel`](@ref)
    - [`cubic_b_spline_kernel`](@ref)
- `model::AbstractConstitutiveModel`: Constitutive model defining the material behavior. \\
    (default: `LinearElastic()`) \\
    Only the following model can be used:
    - [`LinearElastic`](@ref)
- `zem::AbstractZEMStabilization`: Algorithm of zero-energy mode stabilization. \\
    (default: [`ZEMSilling`](@ref)) \\
    The following algorithms can be used:
    - [`ZEMSilling`](@ref)
    - [`ZEMWan`](@ref)
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
    zem::ZEM
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

@params CRMaterial CPointParameters

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
    zem_stiffness_rotated::MArray{NTuple{4,3},Float64,4,81}
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

function init_field(::CRMaterial, ::AbstractTimeSolver, system::BondSystem,
                    ::Val{:zem_stiffness_rotated})
    return zero(MArray{NTuple{4,3},Float64,4,81})
end

function calc_deformation_gradient(storage::CRStorage, system::BondSystem, ::CRMaterial,
                                   ::CPointParameters, i)
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
                                     params::CPointParameters, defgrad_res, Δt, i)
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
end

function calc_zem_stiffness_tensor!(storage::CRStorage, system::BondSystem, mat::CRMaterial,
                                    params::CPointParameters, zem::ZEMWan, defgrad_res, i)
    (; Kinv) = defgrad_res
    (; rotation, zem_stiffness_rotated) = storage
    R = get_tensor(rotation, i)
    C_1 = calc_rotated_zem_stiffness_tensor!(zem_stiffness_rotated, params.C, Kinv, R)
    return C_1
end
