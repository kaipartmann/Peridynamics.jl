struct NeoHookeNonlinear <: AbstractConstitutiveModel end

function cauchy_stress(::NeoHookeNonlinear, storage::AbstractStorage,
                       params::AbstractPointParameters, F::SMatrix{3,3})
    J = det(F)
    J < eps() && return zero(SMatrix{3,3})
    C = F' * F
    Cinv = inv(C)
    S = params.G .* (I - 1 / 3 .* tr(C) .* Cinv) .* J^(-2 / 3) .+
        params.K / 4 .* (J^2 - J^(-2)) .* Cinv
    P = F * S
    σ = 1/J .* P * F'
    return σ
end

struct SaintVenantKirchhoff <: AbstractConstitutiveModel end

function cauchy_stress(::SaintVenantKirchhoff, storage::AbstractStorage,
                       params::AbstractPointParameters, F::SMatrix{3,3})
    J = det(F)
    J < eps() && return zero(SMatrix{3,3})
    E = 0.5 .* (F' * F - I)
    S = params.λ * tr(E) * I + 2 * params.μ * E
    P = F * S
    σ = 1/J .* P * F'
    return σ
end
