@doc raw"""
    linear_kernel(δ, L)

A linear kernel function ``\omega`` (also called influence function) used for weighting the
bonds in a family. The kernel function is defined as
```math
\omega = \frac{\delta}{L} \; ,
```
with the horizon ``\delta`` and the initial bond length ``L``.
"""
@inline linear_kernel(δ, L) = δ / L

@doc raw"""
    cubic_b_spline_kernel(δ, L)

A cubic B-spline kernel function ``\omega`` used for weighting the bonds in a family.
The kernel function is defined as
```math
\begin{aligned}
\xi &= \frac{L}{\delta} \; , \\
\omega &= \left\{
    \begin{array}{ll}
        \frac{2}{3} - 4 \xi^2 + 4 \xi^3 & \quad \text{if} \; 0 < \xi \leq 0.5 \; , \\[3pt]
        \frac{4}{3} - 4 \xi + 4 \xi^2 - \frac{4}{3} \xi^3 & \quad \text{if} \; 0.5 < \xi \leq 1 \; , \\[3pt]
        0 & \quad \text{else} \; ,
    \end{array}
    \right.
\end{aligned}
```
with the horizon ``\delta`` and the initial bond length ``L``.
"""
@inline function cubic_b_spline_kernel(δ, L)
    ξ = L / δ
    if 0 < ξ ≤ 0.5
        return 2/3 - 4 * ξ^2 + 4 * ξ^3
    elseif 0.5 < ξ ≤ 1
        return 4/3 - 4 * ξ + 4 * ξ^2 - 4/3 * ξ^3
    else
        return 0
    end
end
