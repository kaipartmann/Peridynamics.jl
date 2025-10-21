@doc raw"""
    LinearElastic

Linear elastic constitutive model that can be specified when using a [`CMaterial`](@ref) and
[`BACMaterial`](@ref).
The first Piola-Kirchhoff stress ``\boldsymbol{P}`` is given by
```math
\boldsymbol{P} = \mathbb{C} : \boldsymbol{E} \; ,
```
with the elastic stiffness tensor  ``\mathbb{C}`` and the Green-Lagrange strain tensor
``\boldsymbol{E}`` with
```math
\boldsymbol{E} = \frac{1}{2} \left( \boldsymbol{F}^{\top} \boldsymbol{F} - \boldsymbol{I}
                             \right) \; .
```
"""
struct LinearElastic <: AbstractConstitutiveModel end

function first_piola_kirchhoff(::LinearElastic, storage::AbstractStorage,
                               params::AbstractPointParameters, F::SMatrix{3,3,T,9}) where T
    E = 0.5 .* (F' * F - I)
    Evoigt = SVector{6,Float64}(E[1,1], E[2,2], E[3,3], 2 * E[2,3], 2 * E[3,1], 2 * E[1,2])
    Cvoigt = get_hooke_matrix_voigt(params.nu, params.λ, params.μ)
    Svoigt = Cvoigt * Evoigt
    # Convert second Piola-Kirchhoff stress from Voigt to tensor form
    S = SMatrix{3,3,Float64,9}(Svoigt[1], Svoigt[6], Svoigt[5],
                               Svoigt[6], Svoigt[2], Svoigt[4],
                               Svoigt[5], Svoigt[4], Svoigt[3])
    # First Piola-Kirchhoff stress: P = F * S
    P = F * S
    return P
end

function strain_energy_density(::LinearElastic, storage::AbstractStorage,
                               params::AbstractPointParameters, F::SMatrix{3,3,T,9}) where T
    E = 0.5 .* (F' * F - I)
    # For energy density, use the standard form: Ψ = 0.5 * λ * tr(E)^2 + μ * tr(E*E)
    # This is equivalent to the Saint-Venant-Kirchhoff model
    Ψ = 0.5 * params.λ * tr(E)^2 + params.μ * tr(E * E)
    return Ψ
end

function get_hooke_matrix_voigt(nu, λ, μ)
    a = (1 - nu) * λ / nu
    Cvoigt = SMatrix{6,6,Float64,36}(
        a, λ, λ, 0, 0, 0,
        λ, a, λ, 0, 0, 0,
        λ, λ, a, 0, 0, 0,
        0, 0, 0, μ, 0, 0,
        0, 0, 0, 0, μ, 0,
        0, 0, 0, 0, 0, μ
    )
    return Cvoigt
end

function get_hooke_matrix(nu, λ, μ)
    Cvoigt = get_hooke_matrix_voigt(nu, λ, μ)
    C = zero(MArray{NTuple{4,3},Float64,4,81})
    voigt_map = @SVector [(1,1), (2,2), (3,3), (2,3), (1,3), (1,2)]
    for I in 1:6, J in 1:6
        i1, i2 = voigt_map[I]
        j1, j2 = voigt_map[J]
        C[i1, i2, j1, j2] = Cvoigt[I, J]
        C[i2, i1, j1, j2] = Cvoigt[I, J]
        C[i1, i2, j2, j1] = Cvoigt[I, J]
        C[i2, i1, j2, j1] = Cvoigt[I, J]
    end
    return C
end


@doc raw"""
    NeoHooke

Neo-Hookean constitutive model that can be specified when using a [`CMaterial`](@ref) and
[`BACMaterial`](@ref).
The first Piola-Kirchhoff stress ``\boldsymbol{P}`` is given by
```math
\begin{aligned}
\boldsymbol{C} &= \boldsymbol{F}^{\top} \boldsymbol{F} \; , \\
\boldsymbol{S} &= \mu \left( \boldsymbol{I} - \boldsymbol{C}^{-1} \right)
    + \lambda \log(J) \boldsymbol{C}^{-1} \; , \\
\boldsymbol{P} &= \boldsymbol{F} \, \boldsymbol{S} \; ,
\end{aligned}
```
with the deformation gradient ``\boldsymbol{F}``, the right Cauchy-Green deformation tensor
``\boldsymbol{C}``, the Jacobian ``J = \mathrm{det}(\boldsymbol{F})``, the second
Piola-Kirchhoff stress ``\boldsymbol{S}``, and the first and second Lamé parameters
``\lambda`` and ``\mu``.
"""
struct NeoHooke <: AbstractConstitutiveModel end

function first_piola_kirchhoff(::NeoHooke, storage::AbstractStorage,
                               params::AbstractPointParameters, F::SMatrix{3,3,T,9}) where T
    J = det(F)
    Cinv = inv(F' * F)
    S = params.μ * (I - Cinv) + params.λ * log(J) * Cinv
    P = F * S
    return P
end

function strain_energy_density(::NeoHooke, storage::AbstractStorage,
                               params::AbstractPointParameters, F::SMatrix{3,3,T,9}) where T
    J = det(F)
    C = F' * F
    I₁ = tr(C)
    Ψ = 0.5 * params.μ * (I₁ - 3) - params.μ * log(J) + 0.5 * params.λ * log(J)^2
    return Ψ
end

@doc raw"""
    MooneyRivlin

Mooney-Rivlin constitutive model that can be specified when using a [`CMaterial`](@ref) and
[`BACMaterial`](@ref).
The first Piola-Kirchhoff stress ``\boldsymbol{P}`` is given by
```math
\begin{aligned}
\boldsymbol{C} &= \boldsymbol{F}^{\top} \boldsymbol{F} \; , \\
\boldsymbol{S} &= G \left( \boldsymbol{I} - \frac{1}{3} \mathrm{tr}(\boldsymbol{C})
                           \boldsymbol{C}^{-1} \right) \cdot J^{-\frac{2}{3}}
                + \frac{K}{4} \left( J^2 - J^{-2} \right) \boldsymbol{C}^{-1} \; , \\
\boldsymbol{P} &= \boldsymbol{F} \, \boldsymbol{S} \; ,
\end{aligned}
```
with the deformation gradient ``\boldsymbol{F}``, the right Cauchy-Green deformation tensor
``\boldsymbol{C}``, the Jacobian ``J = \mathrm{det}(\boldsymbol{F})``, the second
Piola-Kirchhoff stress ``\boldsymbol{S}``, the shear modulus ``G``, and the
bulk modulus ``K``.

# Error handling
If the Jacobian ``J`` is smaller than the machine precision `eps()` or a `NaN`, the first
Piola-Kirchhoff stress tensor is defined as ``\boldsymbol{P} = \boldsymbol{0}``.
"""
struct MooneyRivlin <: AbstractConstitutiveModel end

function first_piola_kirchhoff(::MooneyRivlin, storage::AbstractStorage,
                               params::AbstractPointParameters, F::SMatrix{3,3,T,9}) where T
    J = det(F)
    J < eps() && return zero(SMatrix{3,3,T,9})
    isnan(J) && return zero(SMatrix{3,3,T,9})
    C = F' * F
    Cinv = inv(C)
    S = params.G .* (I - 1 / 3 .* tr(C) .* Cinv) .* J^(-2 / 3) .+
        params.K / 4 .* (J^2 - J^(-2)) .* Cinv
    P = F * S
    return P
end

function strain_energy_density(::MooneyRivlin, storage::AbstractStorage,
                               params::AbstractPointParameters, F::SMatrix{3,3,T,9}) where T
    J = det(F)
    J < eps() && return zero(T)
    isnan(J) && return zero(T)
    C = F' * F
    I₁ = tr(C)
    Ψ = 0.5 * params.G * (I₁ * J^(-2 / 3) - 3) + params.K / 8 * (J^2 + J^(-2) - 2)
    return Ψ
end

@doc raw"""
    SaintVenantKirchhoff

Saint-Venant-Kirchhoff constitutive model that can be specified when using a
[`CMaterial`](@ref) and [`BACMaterial`](@ref).
The first Piola-Kirchhoff stress ``\boldsymbol{P}`` is given by
```math
\begin{aligned}
\boldsymbol{E} &= \frac{1}{2} \left( \boldsymbol{F}^{\top} \boldsymbol{F} - \boldsymbol{I}
                              \right) \; , \\
\boldsymbol{S} &= \lambda \, \mathrm{tr}(\boldsymbol{E}) \, \boldsymbol{I}
                + 2 \mu \boldsymbol{E} \; , \\
\boldsymbol{P} &= \boldsymbol{F} \, \boldsymbol{S} \; ,
\end{aligned}
```
with the deformation gradient ``\boldsymbol{F}``, the Green-Lagrange strain tensor
``\boldsymbol{E}``, the second Piola-Kirchhoff stress ``\boldsymbol{S}``, and the first and
second Lamé parameters ``\lambda`` and ``\mu``.
"""
struct SaintVenantKirchhoff <: AbstractConstitutiveModel end

function first_piola_kirchhoff(::SaintVenantKirchhoff, storage::AbstractStorage,
                               params::AbstractPointParameters, F::SMatrix{3,3,T,9}) where T
    E = 0.5 .* (F' * F - I)
    S = params.λ * tr(E) * I + 2 * params.μ * E
    P = F * S
    return P
end

function strain_energy_density(::SaintVenantKirchhoff, storage::AbstractStorage,
                               params::AbstractPointParameters, F::SMatrix{3,3,T,9}) where T
    E = 0.5 .* (F' * F - I)
    Ψ = 0.5 * params.λ * tr(E)^2 + params.μ * tr(E * E)
    return Ψ
end
