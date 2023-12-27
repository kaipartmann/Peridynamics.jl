function Base.show(io::IO, ::MIME"text/plain", bc::T) where {T <: AbstractBC}
    msg = string(length(bc.point_id_set), "-points ", T)
    if bc.dim == 1
        msg *= " in x-direction (dim=1)"
    elseif bc.dim == 2
        msg *= " in y-direction (dim=2)"
    elseif bc.dim == 3
        msg *= " in z-direction (dim=3)"
    end
    print(io, msg)
    return nothing
end

function Base.show(io::IO, ::MIME"text/plain", ic::T) where {T <: AbstractIC}
    msg = string(length(ic.point_id_set), "-points ", T)
    if ic.dim == 1
        msg *= " in x-direction (dim=1)"
    elseif ic.dim == 2
        msg *= " in y-direction (dim=2)"
    elseif ic.dim == 3
        msg *= " in z-direction (dim=3)"
    end
    print(io, msg)
    return nothing
end

@doc raw"""
    VelocityBC <: AbstractBC

Velocity boundary condition. The value of the velocity is calculated every time step with
the function `fun` and applied to the dimension `dim` of the points specified by
`point_id_set`.

# Fields
- `fun::Function`: function `f(t)` with current time `t` as argument, that calculates the
value of the velocity for each timestep
- `point_id_set::Vector{Int}`: point-id set with all points for which the boundary condition
gets applied every timestep
- `dim::Int`: dimension on which the boundary condition gets applied to. Possible values:
- x-direction: `dim=1`
- y-direction: `dim=2`
- z-direction: `dim=3`

# Examples

The constant velocity $v = 0.1$ in $y$-direction gets applied to the first $10$ points:
```julia-repl
julia> VelocityBC(t -> 0.1, 1:10, 2)
VelocityBC(var"#7#8"(), [1, 2, 3, 4, 5, 6, 7, 8, 9, 10], 2)
```
"""
struct VelocityBC{F <: Function} <: AbstractBC
    fun::F
    point_id_set::Vector{Int}
    dim::Int
end

@doc raw"""
    ForceDensityBC <: AbstractBC

Force density boundary condition. The value of the force density is calculated every time
step with the function `fun` and applied to the dimension `dim` of the points specified by
`point_id_set`.

# Fields
- `fun::Function`: function `f(t)` with current time `t` as argument, that calculates the
value of the force density for each timestep
- `point_id_set::Vector{Int}`: point-id set with all points for which the boundary condition
gets applied every timestep
- `dim::Int`: dimension on which the boundary condition gets applied to. Possible values:
- x-direction: `dim=1`
- y-direction: `dim=2`
- z-direction: `dim=3`
"""
struct ForceDensityBC{F <: Function} <: AbstractBC
    fun::F
    point_id_set::Vector{Int}
    dim::Int
end

@doc raw"""
    VelocityIC <: AbstractIC

Velocity initial condition. The value `val` of the velocity is applied as initial condition
to the dimension `dim` of the points specified by `point_id_set`.

# Fields
- `val::Float64`: value of the velocity
- `point_id_set::Vector{Int}`: point-id set with all points for which the initial condition
gets applied to
- `dim::Int`: dimension on which the initial condition gets applied to. Possible values:
- x-direction: `dim=1`
- y-direction: `dim=2`
- z-direction: `dim=3`
"""
struct VelocityIC <: AbstractIC
    val::Float64
    point_id_set::Vector{Int}
    dim::Int
end

@doc raw"""
    PosDepVelBC <: AbstractBC

Position dependent velocity boundary condition. The value of the force density is calculated
every time step with the function `fun` and applied to the dimension `dim` of the points
specified by `point_id_set`.

# Fields
- `fun::Function`: function `f(x, y, z, t)` with x-, y-, z-position of the point and current
time `t` as arguments, calculates the value of the force density for each timestep
dependent of the point position
- `point_id_set::Vector{Int}`: point-id set with all points for which the boundary condition
gets applied to
- `dim::Int`: dimension on which the initial condition gets applied to. Possible values:
- x-direction: `dim=1`
- y-direction: `dim=2`
- z-direction: `dim=3`
"""
struct PosDepVelBC{F <: Function} <: AbstractBC
    fun::F
    point_id_set::Vector{Int}
    dim::Int
end

function apply_bcs!(body::AbstractPDBody, bcs::Vector{<:AbstractBC}, t::Float64)
    @sync for bc in bcs
        Threads.@spawn apply_boundarycondition!(body, bc, t)
    end
    return nothing
end

function apply_boundarycondition!(body::AbstractPDBody, bc::VelocityBC, t::Float64)
    value = bc.fun(t)
    dim = bc.dim
    @simd for i in bc.point_id_set
        @inbounds body.velocity_half[dim, i] = value
    end
    return nothing
end

function apply_boundarycondition!(body::AbstractPDBody, bc::ForceDensityBC, t::Float64)
    value = bc.fun(t)
    dim = bc.dim
    @simd for i in bc.point_id_set
        @inbounds body.b_ext[dim, i] = value
    end
    return nothing
end

function apply_boundarycondition!(body::AbstractPDBody, bc::PosDepVelBC, t::Float64)
    dim = bc.dim
    @simd for i in bc.point_id_set
        @inbounds body.velocity_half[dim, i] = bc.fun(body.position[1, i], body.position[2, i],
                                            body.position[3, i], t)
    end
    return nothing
end

function apply_ics!(body::AbstractPDBody, ics::Vector{<:AbstractIC})
    @sync for ic in ics
        Threads.@spawn apply_initialcondition!(body, ic)
    end
    return nothing
end

function apply_initialcondition!(body::AbstractPDBody, ic::VelocityIC)
    dim = ic.dim
    @simd for i in ic.point_id_set
        @inbounds body.velocity[dim, i] = ic.val
    end
    return nothing
end
