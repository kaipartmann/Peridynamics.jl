@doc raw"""
    const_one_kernel(δ, L)

A kernel function ``\omega`` (also called influence function) used for weighting the
bonds in a family. The kernel function is simply defined as a constant value 1:
```math
\omega = 1
```
"""
@inline const_one_kernel(δ, L) = 1

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

@doc raw"""
    cubic_b_spline_kernel_norm(δ, L)

A cubic B-spline kernel function ``\omega`` used for weighting the bonds in a family.
The kernel function is defined as
```math
\begin{aligned}
\xi &= \frac{L}{\delta} \; , \\
\omega &= \frac{8}{\pi \, \delta^3} \cdot \left\{
    \begin{array}{ll}
        1 - 6 \xi^2 + 6 \xi^3 & \quad \text{if} \; 0 < \xi \leq 0.5 \; , \\[3pt]
        2 (1 - \xi)^3 & \quad \text{if} \; 0.5 < \xi \leq 1 \; , \\[3pt]
        0 & \quad \text{else} \; ,
    \end{array}
    \right.
\end{aligned}
```
with the horizon ``\delta`` and the initial bond length ``L``. This kernel is properly
normalized to satisfy the condition ``\int_{\mathcal{H}(X)} \omega(\Delta X) dV' = 1``.
"""
@inline function cubic_b_spline_kernel_norm(δ, L)
    ξ = L / δ
    C = 8.0/(π * δ^3)
    if 0 < ξ ≤ 0.5
        return C * (1 - 6 * ξ^2 + 6 * ξ^3)
    elseif 0.5 < ξ ≤ 1
        return C * 2 * (1 - ξ)^3
    else
        return 0
    end
end
