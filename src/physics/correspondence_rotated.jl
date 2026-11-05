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
    @inherit VelocityVerletFields DynamicRelaxationFields
    @inherit BondFracFields
    @lth velocity_half::PointVector
    @htl b_int::PointVector
    unrotated_stress::PointTensor
    defgrad::PointTensor
    cauchy_stress::PointTensor
    von_mises_stress::PointScalar
    strain_energy_density::PointScalar
    left_stretch::PointTensor = I
    rotation::PointTensor = I
    zem_stiffness_rotated::MArray{NTuple{4,3},Float64,4,81}
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
custom_field(::Type{<:CRStorage}, ::Val{:hydrostatic_stress}) = true

# calculate the strain energy density from the deformation gradient just when exporting
function export_field(::Val{:strain_energy_density}, mat::CRMaterial, system::BondSystem,
                      storage::AbstractStorage, paramsetup::AbstractParameterSetup, t)
    model = mat.constitutive_model
    for i in each_point_idx(system)
        params = get_params(paramsetup, i)
        F = get_tensor(storage.defgrad, i)
        storage.strain_energy_density[i] = strain_energy_density(model, storage, params, F,
                                                                 i)
    end
    return storage.strain_energy_density
end
