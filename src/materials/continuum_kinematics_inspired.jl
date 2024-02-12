struct CKIMaterial <: AbstractMaterial
    δ::Float64
    rho::Float64
    E::Float64
    nu::Float64
    G::Float64
    K::Float64
    C1::Float64
    C2::Float64
    C3::Float64
    Gc::Float64
    εc::Float64
end

function CKIMaterial(; horizon::Real, rho::Real, E::Real, nu::Real,
                     Gc::Real=-1, epsilon_c::Real=-1, C1::Real=0,
                     C2::Real=0, C3::Real=0)
    K = E / (3 * (1 - 2 * nu))
    G = E / (2 * (1 + nu))
    λ = E * nu / ((1 + nu) * (1 - 2nu))
    μ = G
    if (Gc !== -1) && (epsilon_c == -1)
        epsilon_c = sqrt(5.0 * Gc / (9.0 * K * horizon))
    elseif (Gc == -1) && (epsilon_c !== -1)
        Gc = 9.0 / 5.0 * K * horizon * epsilon_c^2
    elseif (Gc !== -1) && (epsilon_c !== -1)
        msg = "duplicate definition: define either Gc or epsilon_c, not both!"
        throw(ArgumentError(msg))
    elseif (Gc == -1) && (epsilon_c == -1)
        throw(ArgumentError("either Gc or epsilon_c have to be defined!"))
    end
    if C1 == 0 && C2 == 0 && C3 == 0
        C1 = 30 / π * μ / horizon^4
        C2 = 0.0
        C3 = 32 / π^4 * (λ - μ) / horizon^12
    else
        msg = "CPD parameters choosen manually!\n"
        msg *= "Be careful when adjusting CPD parameters to avoid unexpected outcomes!"
        @warn msg
    end
    return CKIMaterial(horizon, rho, E, nu, G, K, C1, C2, C3, Gc, epsilon_c)
end
