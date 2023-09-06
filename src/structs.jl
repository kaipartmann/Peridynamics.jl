
"""
    PreCrack(point_id_set_a::Vector{Int}, point_id_set_b::Vector{Int})

Definition of an preexisting crack in the model. Points in `point_id_set_a` cannot have
interactions with points in `point_id_set_b`.

# Fields
- `point_id_set_a::Vector{Int}`: first point-id set
- `point_id_set_b::Vector{Int}`: second point-id set
"""
struct PreCrack
    point_id_set_a::Vector{Int}
    point_id_set_b::Vector{Int}
end

function Base.show(io::IO, ::MIME"text/plain", precrack::PreCrack)
    println(io, typeof(precrack), ":")
    println(io, " ", length(precrack.point_id_set_a), " points in set a")
    print(io, " ", length(precrack.point_id_set_b), " points in set b")
    return nothing
end

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
struct VelocityBC <: AbstractBC
    fun::Function
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
struct ForceDensityBC <: AbstractBC
    fun::Function
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
struct PosDepVelBC <: AbstractBC
    fun::Function
    point_id_set::Vector{Int}
    dim::Int
end

"""
    ExportSettings

Export settings.

# Fields
- `path::String`: path where results will be saved
- `exportfreq::Int`: export frequency, will export every `exportfreq`-th timestep
- `resultfile_prefix::String`: prefix of the result-filename
- `logfile::String`: name of logfile
- `exportflag::Bool`: disable export for a simulation where saved results are not needed

---
```julia
ExportSettings([path::String, freq::Int])
```

Create `ExportSettings` only by `path` and `freq`. If no arguments are specified, the
`exportflag` will be set to `false` and export disabled.

# Arguments
- `path::String`: path where results will be saved
- `freq::Int`: export frequency
"""
mutable struct ExportSettings
    path::String
    exportfreq::Int
    resultfile_prefix::String
    logfile::String
    exportflag::Bool
end

ExportSettings() = ExportSettings("", 0, "", "", false)
ExportSettings(path::String, freq::Int) = ExportSettings(path, freq, "", "", true)

function Base.show(io::IO, ::MIME"text/plain", es::ExportSettings)
    print(io, typeof(es))
    return nothing
end
