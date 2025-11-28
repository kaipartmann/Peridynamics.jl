"""
    CRMaterial(; kernel, model, zem, dmgmodel, maxdmg)

The same as the [`CMaterial`](@ref) but with rotation of the stress tensor for large
deformation simulations, therefore not all models are supported.

Supported models:
- `SaintVenantKirchhoff`
- `LinearElastic`

Please take a look at the [`CMaterial`](@ref) docs for more information about the
material!
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
                    model::AbstractConstitutiveModel=SaintVenantKirchhoff(),
                    zem::AbstractZEMStabilization=ZEMSilling(),
                    dmgmodel::AbstractDamageModel=CriticalStretch(), maxdmg::Real=0.85)
    if !(typeof(model) <: Union{SaintVenantKirchhoff,LinearElastic})
        msg = "model `$(typeof(model))` is currently not supported for `CRMaterial`!\n"
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
    @pointfield n_active_bonds::Vector{Int}
    @pointfield unrotated_stress::Matrix{Float64}
    @pointfield defgrad::Matrix{Float64}
    @pointfield cauchy_stress::Matrix{Float64}
    @pointfield von_mises_stress::Vector{Float64}
    @pointfield strain_energy_density::Vector{Float64}
    @pointfield left_stretch::Matrix{Float64}
    @pointfield rotation::Matrix{Float64}
    bond_active::Vector{Bool}
    zem_stiffness_rotated::MArray{NTuple{4,3},Float64,4,81}
end

function init_field(::CRMaterial, ::AbstractTimeSolver, system::BondSystem,
                    ::Val{:velocity_half})
    return zeros(3, get_n_points(system))
end

function init_field(::CRMaterial, ::AbstractTimeSolver, system::BondSystem, ::Val{:b_int})
    return zeros(3, get_n_points(system))
end

function init_field(::CRMaterial, ::AbstractTimeSolver, system::BondSystem,
                    ::Val{:unrotated_stress})
    return zeros(9, get_n_loc_points(system))
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

function calc_deformation_gradient!(storage::CRStorage, system::BondSystem, ::CRMaterial,
                                    ::CPointParameters, i)
    (; bonds, volume) = system
    (; bond_active) = storage
    K = zero(SMatrix{3,3,Float64,9})
    _F = zero(SMatrix{3,3,Float64,9})
    _Ḟ = zero(SMatrix{3,3,Float64,9})
    for bond_id in each_bond_idx(system, i)
        bond = bonds[bond_id]
        j = bond.neighbor
        ΔXij = get_vector_diff(system.position, i, j)
        Δxij = get_vector_diff(storage.position, i, j)
        Δvij = get_vector_diff(storage.velocity_half, i, j)
        ωij = kernel(system, bond_id) * bond_active[bond_id]
        temp = ωij * volume[j]
        ΔXijt = ΔXij'
        K += temp * (ΔXij * ΔXijt)
        _F += temp * (Δxij * ΔXijt)
        _Ḟ += temp * (Δvij * ΔXijt)
    end
    Kinv = inv(K)
    F = _F * Kinv
    Ḟ = _Ḟ * Kinv
    Peridynamics.update_tensor!(storage.defgrad, i, F)
    return (; F, Ḟ, Kinv)
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
    σ = get_tensor(storage.unrotated_stress, i)
    σₙ₊₁ = σ + 2 * params.G * Δεᵈᵉᵛ + params.K * Δθ * I
    update_tensor!(storage.unrotated_stress, i, σₙ₊₁)
    T = rotate_stress(storage, σₙ₊₁, i)
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

# calculate the von Mises stress from the Cauchy stress tensor just when exporting
function export_field(::Val{:cauchy_stress}, ::CRMaterial, system::BondSystem,
                      storage::AbstractStorage, ::AbstractParameterSetup, t)
    for i in each_point_idx(system)
        σ = get_tensor(storage.unrotated_stress, i)
        T = rotate_stress(storage::AbstractStorage, σ, i)
        Peridynamics.update_tensor!(storage.cauchy_stress, i, T)
    end
    return storage.cauchy_stress
end

# calculate the von Mises stress from the Cauchy stress tensor just when exporting
function export_field(::Val{:von_mises_stress}, ::CRMaterial, system::BondSystem,
                      storage::AbstractStorage, ::AbstractParameterSetup, t)
    for i in each_point_idx(system)
        σ = get_tensor(storage.unrotated_stress, i)
        T = rotate_stress(storage::AbstractStorage, σ, i)
        storage.von_mises_stress[i] = von_mises_stress(T)
    end
    return storage.von_mises_stress
end

# calculate the von hydrostatic stress from the Cauchy stress tensor just when exporting,
# use the `von_mises_stress` field to store the hydrostatic stress
function export_field(::Val{:hydrostatic_stress}, ::CRMaterial,
                      system::BondSystem, storage::CRStorage, ::AbstractParameterSetup, t)
    for i in each_point_idx(system)
        σ = get_tensor(storage.unrotated_stress, i)
        T = rotate_stress(storage::AbstractStorage, σ, i)
        storage.von_mises_stress[i] = 1/3 * (T[1,1] + T[2,2] + T[3,3])
    end
    return storage.von_mises_stress
end
custom_field(::Type{CRStorage}, ::Val{:hydrostatic_stress}) = true

# calculate the strain energy density from the deformation gradient just when exporting
function export_field(::Val{:strain_energy_density}, mat::CRMaterial, system::BondSystem,
                      storage::AbstractStorage, paramsetup::AbstractParameterSetup, t)
    model = mat.constitutive_model
    for i in each_point_idx(system)
        params = get_params(paramsetup, i)
        F = get_tensor(storage.defgrad, i)
        storage.strain_energy_density[i] = strain_energy_density(model, storage, params, F)
    end
    return storage.strain_energy_density
end
