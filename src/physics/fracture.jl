"""
    CriticalStretch

A damage model based on the stretch of the bond. The bond is considered to be broken
if the stretch exceeds a critical value.
The critical value can be defined via the fracture energy `Gc` or the critical stretch `εc`
using the [`material!`](@ref) function. The damage model is defined globally for the whole
body as part of the material.
"""
struct CriticalStretch <: AbstractDamageModel end

@inline fracture_kwargs() = (:Gc, :epsilon_c)

"""
    failure_permit!(body, set_name, fail_permit)

$(internal_api_warning())

Set the failure permission for points of the set `set_name` of a `body`.

# Arguments

- `body::AbstractBody`: [`Body`](@ref) where the failure permission will be set.
- `set_name::Symbol`: The name of a point set of this body.
- `fail_permit::Bool`: If `true`, failure is allowed, and if `false` then no bonds of this
    point are allowed to break during the simulation.

!!! danger "Overwriting failure permission with `material!` and `failure_permit!`"
    The function `material!` calls `failure_permit!`, so if it is used afterwards,
    previously set failure permissions might be overwritten!

# Throws

- Error if the body does not contain a set with `set_name`.
"""
function failure_permit! end

function failure_permit!(body::AbstractBody, set_name::Symbol, fail_permit::Bool)
    check_if_set_is_defined(body.point_sets, set_name)
    body.fail_permit[body.point_sets[set_name]] .= fail_permit
    return nothing
end

"""
    no_failure!(body::AbstractBody, set_name::Symbol)
    no_failure!(body::AbstractBody)

Disallow failure for all points of the point set `set_name` of the `body`.
If no `set_name` is specified, failure is prohibited for the whole `body`.

# Arguments

- `body::AbstractBody`: [`Body`](@ref) for which failure is prohibited.
- `set_name::Symbol`: The name of a point set of this body.

!!! danger "Overwriting failure permission with `material!` and `no_failure!`"
    The function `material!` sets failure permissions due to the provided input parameters,
    so if it is used afterwards, previously set failure prohibitions might be overwritten!

# Throws

- Error if the body does not contain a set with `set_name`.

# Examples

```julia-repl
julia> no_failure!(body)

julia> body
1000-point Body{BBMaterial{NoCorrection}}:
  1 point set(s):
    1000-point set `all_points`
  1000 points with failure prohibited
```
"""
function no_failure! end

function no_failure!(body::AbstractBody, set_name::Symbol)
    check_if_set_is_defined(body.point_sets, set_name)
    points_without_material = findfirst(x -> x == 0,
                                        body.params_map[body.point_sets[set_name]])
    if points_without_material !== nothing
        msg = string("Not all points of point set \":", set_name,
                     "\" have material parameters!\n")
        msg *= "Please use the `material!` function to define them before prohibiting "
        msg *= "failure, otherwise failure permissions might be overwritten!\n"
        throw(ArgumentError(msg))
    end
    failure_permit!(body, set_name, false)
    return nothing
end

function no_failure!(body::AbstractBody)
    no_failure!(body, :all_points)
    return nothing
end


"""
    get_frac_params(::AbstractDamageModel, p::Dict{Symbol,Any}, δ::Float64, K::Float64)

$(internal_api_warning())

Read or calculate the necessary fracture parameters for the specified damage model from the
dictionary created with [`material!`](@ref). This function has to be defined when creating a
new damage model. Otherwise, a default method returns a empty named tuple `(; )`.
"""
function get_frac_params end

function get_frac_params(::CriticalStretch, p::Dict{Symbol,Any}, δ::Float64, K::Float64)
    local Gc::Float64
    local εc::Float64

    if haskey(p, :Gc) && !haskey(p, :epsilon_c)
        Gc = float(p[:Gc])
        εc = sqrt(5.0 * Gc / (9.0 * K * δ))
    elseif !haskey(p, :Gc) && haskey(p, :epsilon_c)
        εc = float(p[:epsilon_c])
        Gc = 9.0 / 5.0 * K * δ * εc^2
    elseif haskey(p, :Gc) && haskey(p, :epsilon_c)
        msg = "insufficient keywords for calculation of fracture parameters!\n"
        msg *= "Define either Gc or epsilon_c, not both!\n"
        throw(ArgumentError(msg))
    else
        Gc = 0.0;
        εc = 0.0;
    end

    return (; Gc, εc)
end

function get_frac_params(::AbstractDamageModel, p, δ, K)
    return (; )
end

"""
    set_failure_permissions!(body, set_name, params)

$(internal_api_warning())

Grant or prohibit failure permission depending on the submitted fracture parameters by
calling [`failure_permit!`](@ref).

If fracture parameters are found, failure is allowed. If no fracture parameters are found,
failure is not allowed.
"""
function set_failure_permissions!(body::AbstractBody, set_name::Symbol,
                                  params::AbstractPointParameters)
    if has_fracture(body.mat, params)
        failure_permit!(body, set_name, true)
    else
        failure_permit!(body, set_name, false)
    end
    return nothing
end

"""
    has_fracture(mat, params)

$(internal_api_warning())

Return `true` if at least one fracture parameter is set `!=0` in `params` and the system
therefore is supposed to have failure allowed or return `false` if not.
"""
function has_fracture(mat::AbstractMaterial, params::AbstractPointParameters)
    return has_fracture(mat.dmgmodel, params)
end

function has_fracture(::CriticalStretch, params::AbstractPointParameters)
    if isapprox(params.Gc, 0; atol=eps()) || isapprox(params.εc, 0; atol=eps())
        return false
    else
        return true
    end
end

function required_fields_fracture(::Type{Material}) where {Material<:AbstractMaterial}
    fields = (req_point_data_fields_fracture(Material)...,
              req_bond_data_fields_fracture(Material)...,
              req_data_fields_fracture(Material)...)
    return fields
end

function req_point_data_fields_fracture(::Type{Material}) where {Material<:AbstractMaterial}
    return ()
end

function req_bond_data_fields_fracture(::Type{Material}) where {Material<:AbstractMaterial}
    return ()
end

function req_data_fields_fracture(::Type{Material}) where {Material<:AbstractMaterial}
    return ()
end

function log_param_property(::Val{:Gc}, param; indentation)
    return msg_qty("critical energy release rate", param.Gc; indentation)
end

function log_param_property(::Val{:εc}, param; indentation)
    return msg_qty("critical stretch", param.εc; indentation)
end
#-------------------------------------------------------------------------------------------
# EquivalentStress damage model
#-------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------
# Predefined equivalent stress functions
#-------------------------------------------------------------------------------------------

"""
    equiv_von_mises_stress(σ::SMatrix{3,3})

Calculate the von Mises equivalent stress from a 3×3 Cauchy stress tensor.
This is the default stress function for [`EquivalentStress`](@ref).

```math
\\sigma_{\\text{vm}} = \\sqrt{\\sigma_x^2 + \\sigma_y^2 + \\sigma_z^2 - \\sigma_x \\sigma_y
    - \\sigma_x \\sigma_z - \\sigma_y \\sigma_z + 3(\\tau_{xy}^2 + \\tau_{xz}^2 + \\tau_{yz}^2)}
```
"""
@inline function equiv_von_mises_stress(σ::SMatrix{3,3,T,9}) where {T}
    σx, σy, σz = σ[1,1], σ[2,2], σ[3,3]
    τxy, τxz, τyz = σ[1,2], σ[1,3], σ[2,3]
    a = σx * σx + σy * σy + σz * σz
    b = -σx * σy - σx * σz - σy * σz
    c = 3 * (τxy * τxy + τxz * τxz + τyz * τyz)
    d = a + b + c
    d < 0 && return zero(T)
    return √d
end

"""
    equiv_hydrostatic_stress(σ::SMatrix{3,3})

Calculate the hydrostatic (mean) stress from a 3×3 Cauchy stress tensor.
Useful for pressure-dependent failure criteria.

```math
\\sigma_{\\text{h}} = \\frac{1}{3}(\\sigma_{11} + \\sigma_{22} + \\sigma_{33}) = \\frac{1}{3}\\text{tr}(\\boldsymbol{\\sigma})
```
"""
@inline function equiv_hydrostatic_stress(σ::SMatrix{3,3,T,9}) where {T}
    return (σ[1,1] + σ[2,2] + σ[3,3]) / 3
end

"""
    equiv_first_principal_stress(σ::SMatrix{3,3})

Calculate the first (maximum) principal stress from a 3×3 Cauchy stress tensor.
This is commonly used in phase-field fracture simulations as a failure criterion
for brittle materials under tension.

The principal stresses are the eigenvalues of the stress tensor, and the first
principal stress is the largest one.
"""
@inline function equiv_first_principal_stress(σ::SMatrix{3,3,T,9}) where {T}
    # Calculate eigenvalues of the symmetric stress tensor
    # For a 3x3 symmetric matrix, we can use the analytical solution
    λ = eigvals(Symmetric(σ))
    return maximum(λ)
end

"""
    equiv_max_shear_stress(σ::SMatrix{3,3})

Calculate the maximum shear stress (Tresca criterion) from a 3×3 Cauchy stress tensor.

```math
\\tau_{\\text{max}} = \\frac{\\sigma_1 - \\sigma_3}{2}
```

where ``\\sigma_1`` and ``\\sigma_3`` are the maximum and minimum principal stresses.
"""
@inline function equiv_max_shear_stress(σ::SMatrix{3,3,T,9}) where {T}
    λ = eigvals(Symmetric(σ))
    σ1, σ3 = maximum(λ), minimum(λ)
    return (σ1 - σ3) / 2
end

"""
    equiv_tensile_principal_stress(σ::SMatrix{3,3})

Calculate an equivalent tensile stress based on the spectral decomposition of the stress
tensor. This follows the phase-field fracture approach by Miehe et al. (2010) where only
the tensile (positive) part of the stress drives fracture.

This function computes an equivalent von Mises-like stress using only the positive
principal stresses:

```math
\\sigma_{\\text{tensile}} = \\sqrt{\\langle\\sigma_1\\rangle^2 + \\langle\\sigma_2\\rangle^2 + \\langle\\sigma_3\\rangle^2}
```

where ``\\langle\\sigma_i\\rangle = \\max(0, \\sigma_i)`` is the positive part (Macaulay bracket).

This ensures that compressive states do not drive damage evolution, which is physically
correct for brittle materials.

# Reference
Miehe, C., Hofacker, M., & Welschinger, F. (2010). A phase field model for rate-independent
crack propagation: Robust algorithmic implementation based on operator splits.
*Computer Methods in Applied Mechanics and Engineering*, 199(45-48), 2765-2778.
"""
@inline function equiv_tensile_principal_stress(σ::SMatrix{3,3,T,9}) where {T}
    λ = eigvals(Symmetric(σ))
    # Use Macaulay brackets: <σᵢ> = max(0, σᵢ)
    σ1_pos = max(zero(T), λ[1])
    σ2_pos = max(zero(T), λ[2])
    σ3_pos = max(zero(T), λ[3])
    # Equivalent tensile stress magnitude
    return sqrt(σ1_pos^2 + σ2_pos^2 + σ3_pos^2)
end

#-------------------------------------------------------------------------------------------

"""
    EquivalentStress(; stress_func=equiv_von_mises_stress, softening=true)

A flexible stress-based damage model designed for correspondence formulations, particularly
[`RKCMaterial`](@ref) and [`RKCRMaterial`](@ref). This model allows specifying a custom
function to compute the equivalent stress from the Cauchy stress tensor, enabling various
failure criteria such as von Mises, hydrostatic, first principal stress, and more.

Unlike [`CriticalStretch`](@ref), which abruptly removes bonds, this model can implement
**gradual bond softening** to maintain stability of the moment matrix during dynamic
fracture simulations.

# Motivation

In correspondence-based peridynamics formulations, the deformation gradient is calculated
from the moment matrix, which depends on the neighborhood structure. When bonds are abruptly
switched off, the moment matrix can become ill-conditioned or singular, leading to numerical
instabilities. This issue is particularly severe during wave-induced crack propagation.

This damage model addresses the issue by:
1. Computing an equivalent stress at each bond using a user-specified function
2. Implementing gradual softening of bond contributions via a damage variable `d ∈ [0, 1]`
3. Degrading the bond's influence on the moment matrix and force density smoothly by the
   factor `(1 - d)`, avoiding sudden changes

# Keywords

- `stress_func::Function`: A function that takes a 3×3 Cauchy stress tensor (`SMatrix{3,3}`)
  and returns a scalar equivalent stress value. Predefined options:
  - [`equiv_von_mises_stress`](@ref) (default): von Mises equivalent stress
  - [`equiv_hydrostatic_stress`](@ref): Hydrostatic (mean) stress
  - [`equiv_first_principal_stress`](@ref): First (maximum) principal stress
  - [`equiv_max_shear_stress`](@ref): Maximum shear stress (Tresca)

- `softening::Bool`: If `true` (default), enables gradual bond softening. If `false`,
  bonds are immediately switched off when the critical stress is exceeded.

# Material parameters

When using this damage model, set the critical equivalent stress `sigma_c` with the
[`material!`](@ref) function:

```julia
material!(body; horizon=δ, rho=ρ, E=E, nu=ν, sigma_c=σ_c)
```

Alternatively, you can specify the critical energy release rate `Gc`:
```julia
material!(body; horizon=δ, rho=ρ, E=E, nu=ν, Gc=Gc)
```

The relationship between `Gc` and `sigma_c` is derived from fracture mechanics:
```math
\\sigma_c = \\sqrt{\\frac{5 G_c E}{9 \\delta}}
```

# Softening behavior

The bond damage evolves according to a linear softening law. Once the equivalent stress
at a bond exceeds the critical stress `σ_c`, the bond damage increases based on the
excess stress. The bond contribution to the moment matrix and force density is scaled
by the phase-field degradation function `g(d) = (1 - d)²`, where `d` is the bond damage:
- `d = 0`: Bond is fully intact, `g(d) = 1`
- `d = 1`: Bond is fully broken, `g(d) = 0`

The quadratic degradation function is derived from phase-field fracture theory and provides
smoother energy release compared to linear degradation.

# Tensile/Compressive splitting

When `tensile_split=true`, the degradation is applied only to the tensile (positive) part
of the stress, following the spectral decomposition approach from Miehe et al. (2010).
This ensures that compressive stress components remain undegraded, which is physically
correct for brittle materials where cracks cannot propagate under pure compression:

```math
\\boldsymbol{P}_{\\text{eff}} = g(d) \\boldsymbol{P}^+ + \\boldsymbol{P}^-
```

Note: The tensile split adds computational cost due to eigendecomposition. For simulations
where tensile stress dominates, setting `tensile_split=false` (default) is faster.

# Examples

Using von Mises stress (default):
```julia
body = Body(RKCMaterial(dmgmodel=EquivalentStress()), pos, vol)
material!(body; horizon=0.01, rho=7850, E=210e9, nu=0.3, sigma_c=500e6)
```

Using first principal stress with tensile splitting (phase-field style):
```julia
body = Body(RKCMaterial(dmgmodel=EquivalentStress(
    stress_func=equiv_tensile_principal_stress,
    tensile_split=true
)), pos, vol)
material!(body; horizon=0.01, rho=7850, E=210e9, nu=0.3, sigma_c=300e6)
```

Custom stress function:
```julia
# Example: tensile failure only (positive first principal stress)
my_stress_func(σ) = max(0.0, equiv_first_principal_stress(σ))
body = Body(RKCMaterial(dmgmodel=EquivalentStress(stress_func=my_stress_func)), pos, vol)
```

# See also

- [`CriticalStretch`](@ref): Stretch-based damage model (simpler but less stable)
- [`equiv_von_mises_stress`](@ref): Von Mises equivalent stress function
- [`equiv_first_principal_stress`](@ref): First principal stress function
- [`equiv_tensile_principal_stress`](@ref): Tensile principal stress function (phase-field style)
- [`equiv_hydrostatic_stress`](@ref): Hydrostatic stress function
- [`equiv_max_shear_stress`](@ref): Maximum shear stress function
"""
struct EquivalentStress{F<:Function} <: AbstractDamageModel
    stress_func::F
    softening::Bool
    tensile_split::Bool
    viscous_damping::Float64   # Local viscous damping coefficient (0.0 = no damping)
    max_damage_rate::Float64   # Maximum damage increment per step (rate limiting)
    residual_stiffness::Float64 # Minimum stiffness factor (prevents full loss of resistance)
end

function EquivalentStress(; stress_func::Function=equiv_von_mises_stress,
                          softening::Bool=true, tensile_split::Bool=false,
                          viscous_damping::Real=0.05,
                          max_damage_rate::Real=0.1,
                          residual_stiffness::Real=0.01)
    return EquivalentStress(stress_func, softening, tensile_split,
                            Float64(viscous_damping), Float64(max_damage_rate),
                            Float64(residual_stiffness))
end

@inline fracture_kwargs_equiv_stress() = (:Gc, :sigma_c)

function get_frac_params(::EquivalentStress, p::Dict{Symbol,Any}, δ::Float64,
                         K::Float64)
    # Get Young's modulus from the parameter dictionary
    E = haskey(p, :E) ? float(p[:E]) : 3 * K * (1 - 2 * p[:nu])

    local Gc::Float64
    local σc::Float64

    if haskey(p, :Gc) && !haskey(p, :sigma_c)
        Gc = float(p[:Gc])
        # Relationship: σ_c = sqrt(5 * Gc * E / (9 * δ))
        # This is derived from peridynamic fracture energy considerations
        σc = sqrt(5.0 * Gc * E / (9.0 * δ))
    elseif !haskey(p, :Gc) && haskey(p, :sigma_c)
        σc = float(p[:sigma_c])
        # Inverse relationship: Gc = 9 * δ * σ_c^2 / (5 * E)
        Gc = 9.0 * δ * σc^2 / (5.0 * E)
    elseif haskey(p, :Gc) && haskey(p, :sigma_c)
        msg = "insufficient keywords for calculation of fracture parameters!\n"
        msg *= "Define either Gc or sigma_c, not both!\n"
        throw(ArgumentError(msg))
    else
        Gc = 0.0
        σc = 0.0
    end

    return (; Gc, σc)
end

function has_fracture(::EquivalentStress, params::AbstractPointParameters)
    if isapprox(params.Gc, 0; atol=eps()) || isapprox(params.σc, 0; atol=eps())
        return false
    else
        return true
    end
end

function log_param_property(::Val{:σc}, param; indentation)
    return msg_qty("critical equivalent stress", param.σc; indentation)
end
