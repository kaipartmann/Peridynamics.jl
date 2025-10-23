@doc raw"""
    SaintVenantKirchhoff

Saint-Venant-Kirchhoff constitutive model that can be specified when using a
[`CMaterial`](@ref) and [`BACMaterial`](@ref).

The strain energy density ``\Psi`` is given by
```math
\Psi = \frac{1}{2} \lambda \, \mathrm{tr}(\boldsymbol{E})^2 + \mu \, \mathrm{tr}(\boldsymbol{E} \cdot \boldsymbol{E}) \; ,
```
with the first and second Lamé parameters ``\lambda`` and ``\mu``, and the Green-Lagrange
strain tensor
```math
\boldsymbol{E} = \frac{1}{2} \left( \boldsymbol{F}^{\top} \boldsymbol{F} - \boldsymbol{I}
                              \right) \; .
```

The first Piola-Kirchhoff stress ``\boldsymbol{P}`` is given by
```math
\begin{aligned}
\boldsymbol{S} &= \lambda \, \mathrm{tr}(\boldsymbol{E}) \, \boldsymbol{I}
                + 2 \mu \boldsymbol{E} \; , \\
\boldsymbol{P} &= \boldsymbol{F} \, \boldsymbol{S} \; ,
\end{aligned}
```
with the deformation gradient ``\boldsymbol{F}`` and the second Piola-Kirchhoff stress
``\boldsymbol{S}``.

!!! note
    This model is equivalent to the [`LinearElastic`](@ref) model, both using the same
    strain energy density function.
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

@doc raw"""
    LinearElastic

Linear elastic constitutive model that can be specified when using a [`CMaterial`](@ref) and
[`BACMaterial`](@ref).

The strain energy density ``\Psi`` is given by
```math
\Psi = \frac{1}{2} \lambda \, \mathrm{tr}(\boldsymbol{E})^2 + \mu \, \mathrm{tr}(\boldsymbol{E} \cdot \boldsymbol{E}) \; ,
```
with the first and second Lamé parameters ``\lambda`` and ``\mu``, and the Green-Lagrange
strain tensor ``\boldsymbol{E}``
```math
\boldsymbol{E} = \frac{1}{2} \left( \boldsymbol{F}^{\top} \boldsymbol{F} - \boldsymbol{I}
                             \right) \; .
```

The first Piola-Kirchhoff stress ``\boldsymbol{P}`` is given by
```math
\begin{aligned}
\boldsymbol{S} &= \mathbb{C} : \boldsymbol{E} \; , \\
\boldsymbol{P} &= \boldsymbol{F} \, \boldsymbol{S} \; ,
\end{aligned}
```
with the deformation gradient ``\boldsymbol{F}``, the elastic stiffness tensor ``\mathbb{C}``,
and the second Piola-Kirchhoff stress ``\boldsymbol{S}``.

!!! note
    This model is equivalent to the Saint-Venant-Kirchhoff model, but uses a
    different implementation based on the elastic stiffness tensor.
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

Compressible Neo-Hookean hyperelastic constitutive model that can be specified when using
a [`CMaterial`](@ref) and [`BACMaterial`](@ref).

The strain energy density ``\Psi`` is given by
```math
\Psi = \frac{1}{2} \mu \left( I_1 - 3 \right) - \mu \log(J) + \frac{1}{2} \lambda \log(J)^2 \; ,
```
with the first invariant ``I_1 = \mathrm{tr}(\boldsymbol{C})`` of the right Cauchy-Green
deformation tensor ``\boldsymbol{C} = \boldsymbol{F}^{\top} \boldsymbol{F}``, the Jacobian
``J = \mathrm{det}(\boldsymbol{F})``, and the first and second Lamé parameters ``\lambda``
and ``\mu``.

The first Piola-Kirchhoff stress ``\boldsymbol{P}`` is given by
```math
\begin{aligned}
\boldsymbol{C} &= \boldsymbol{F}^{\top} \boldsymbol{F} \; , \\
\boldsymbol{S} &= \mu \left( \boldsymbol{I} - \boldsymbol{C}^{-1} \right)
    + \lambda \log(J) \boldsymbol{C}^{-1} \; , \\
\boldsymbol{P} &= \boldsymbol{F} \, \boldsymbol{S} \; ,
\end{aligned}
```
with the deformation gradient ``\boldsymbol{F}`` and the second Piola-Kirchhoff stress
``\boldsymbol{S}``.

# Reference
Treloar, L. R. G. (1943). "The elasticity of a network of long-chain molecules—II."
*Transactions of the Faraday Society*, 39, 241–246.
DOI: [10.1039/TF9433900241](https://doi.org/10.1039/TF9433900241)
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
    NeoHookePenalty

Compressible Neo-Hookean hyperelastic model with a penalty-type volumetric formulation,
suitable for modeling rubber-like and biological materials. Can be specified when using a
[`CMaterial`](@ref) and [`BACMaterial`](@ref).

The strain energy density ``\Psi`` is given by
```math
\Psi = \frac{1}{2} G \left( \bar{I}_1 - 3 \right) + \frac{K}{8} \left( J^2 + J^{-2} - 2 \right) \; ,
```
with the modified first invariant ``\bar{I}_1 = I_1 J^{-2/3}`` where
``I_1 = \mathrm{tr}(\boldsymbol{C})`` is the first invariant of the right Cauchy-Green
deformation tensor ``\boldsymbol{C} = \boldsymbol{F}^{\top} \boldsymbol{F}``, the Jacobian
``J = \mathrm{det}(\boldsymbol{F})``, the shear modulus ``G``, and the bulk modulus ``K``.

The first Piola-Kirchhoff stress ``\boldsymbol{P}`` is given by
```math
\begin{aligned}
\boldsymbol{C} &= \boldsymbol{F}^{\top} \boldsymbol{F} \; , \\
\boldsymbol{S} &= G \left( \boldsymbol{I} - \frac{1}{3} \mathrm{tr}(\boldsymbol{C})
                           \boldsymbol{C}^{-1} \right) J^{-\frac{2}{3}}
                + \frac{K}{4} \left( J^2 - J^{-2} \right) \boldsymbol{C}^{-1} \; , \\
\boldsymbol{P} &= \boldsymbol{F} \, \boldsymbol{S} \; ,
\end{aligned}
```
with the deformation gradient ``\boldsymbol{F}`` and the second Piola-Kirchhoff stress
``\boldsymbol{S}``.

!!! note "Penalty-type volumetric formulation"
    This model uses a penalty-type volumetric term ``\Psi_{\mathrm{vol}} = \frac{K}{8}(J^2 + J^{-2} - 2)``,
    which is computationally efficient and widely used in commercial finite element codes
    for nearly-incompressible materials. The term penalizes volume changes from the
    reference configuration (``J = 1``).

    This differs from other volumetric formulations such as:
    - Standard Neo-Hookean (logarithmic): ``\Psi_{\mathrm{vol}} = -\mu \ln J + \frac{\lambda}{2}\ln^2 J``
    - Simo-Miehe (polyconvex): ``\Psi_{\mathrm{vol}} = \frac{K}{4}(J^2 - 1 - 2\ln J)``

# Error handling
If the Jacobian ``J`` is smaller than the machine precision `eps()` or a `NaN`, the strain
energy density and first Piola-Kirchhoff stress tensor are defined as zero:
``\Psi = 0`` and ``\boldsymbol{P} = \boldsymbol{0}``.
"""
struct NeoHookePenalty <: AbstractConstitutiveModel end

function first_piola_kirchhoff(::NeoHookePenalty, storage::AbstractStorage,
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

function strain_energy_density(::NeoHookePenalty, storage::AbstractStorage,
                               params::AbstractPointParameters, F::SMatrix{3,3,T,9}) where T
    J = det(F)
    J < eps() && return zero(T)
    isnan(J) && return zero(T)
    C = F' * F
    I₁ = tr(C)
    Ψ = 0.5 * params.G * (I₁ * J^(-2 / 3) - 3) + params.K / 8 * (J^2 + J^(-2) - 2)
    return Ψ
end
