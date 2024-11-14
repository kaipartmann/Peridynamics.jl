# peridynamic kernels / influence functions for arbitrary bond systems

@inline linear_kernel(δ, L) = δ / L

@inline function cubic_b_spline(δ, L)
    ξ = L / δ
    if 0 < ξ ≤ 0.5
        return 2/3 - 4 * ξ^2 + 4 * ξ^3
    elseif 0.5 < ξ ≤ 1
        return 4/3 - 4 * ξ + 4 * ξ^2 - 4/3 * ξ^3
    else
        return 0
    end
end
