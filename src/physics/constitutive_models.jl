struct LinearElastic <: AbstractConstitutiveModel end

function first_piola_kirchhoff(::LinearElastic, storage::AbstractStorage,
                               params::AbstractPointParameters, F::SMatrix{3,3,T,9}) where T
    E = 0.5 .* (F' * F - I)
    Evoigt = SVector{6,Float64}(E[1,1], E[2,2], E[3,3], 2 * E[2,3], 2 * E[3,1], 2 * E[1,2])
    Cvoigt = get_hooke_matrix(params.nu, params.λ, params.μ)
    Pvoigt = Cvoigt * Evoigt
    P = SMatrix{3,3,Float64,9}(Pvoigt[1], Pvoigt[6], Pvoigt[5],
                               Pvoigt[6], Pvoigt[2], Pvoigt[4],
                               Pvoigt[5], Pvoigt[4], Pvoigt[3])
    return P
end

function get_hooke_matrix(nu, λ, μ)
    a = (1 - nu) * λ / nu
    CVoigt = SMatrix{6,6,Float64,36}(a, λ, λ, 0, 0, 0, λ, a, λ, 0, 0, 0, λ, λ, a, 0, 0, 0,
                                     0, 0, 0, μ, 0, 0, 0, 0, 0, 0, μ, 0, 0, 0, 0, 0, 0, μ)
    return CVoigt
end

struct NeoHooke <: AbstractConstitutiveModel end

function first_piola_kirchhoff(::NeoHooke, storage::AbstractStorage,
                               params::AbstractPointParameters, F::SMatrix{3,3,T,9}) where T
    J = det(F)
    Cinv = inv(F' * F)
    S = params.μ * (I - Cinv) + params.λ * log(J) * Cinv
    P = F * S
    return P
end

struct NeoHookeNonlinear <: AbstractConstitutiveModel end

function first_piola_kirchhoff(::NeoHookeNonlinear, storage::AbstractStorage,
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

struct SaintVenantKirchhoff <: AbstractConstitutiveModel end

function first_piola_kirchhoff(::SaintVenantKirchhoff, storage::AbstractStorage,
                               params::AbstractPointParameters, F::SMatrix{3,3,T,9}) where T
    E = 0.5 .* (F' * F - I)
    S = params.λ * tr(E) * I + 2 * params.μ * E
    P = F * S
    return P
end
