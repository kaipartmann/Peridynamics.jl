const JOB_KWARGS = (:path, :freq, :fields)

"""
    Job(spatial_setup, time_solver; kwargs...)

Job that contains all the information required for a peridynamic simulation

# Arguments

- `spatial_setup::AbstractSpatialSetup`: Body or Multibody setup for the simulation
- `time_solver::AbstractTimeSolver`: Method for calculating discrete time steps

# Keywords

- `path::String`: Storage path for results
- `freq::Int`: Frequency of time steps that are exported
- `fields::NTuple{N,Symbol}`: Exported fields
             Possible export fields depend on the selected material model.
             See material type documentation.
             Default export fields: `(:displacement, :damage)`

# Throws

- Error if keyword is not allowed

# Example

```julia-repl
julia> b = Body(BBMaterial(), pos, vol)

julia> vv = VelocityVerlet(steps=2000)

julia> job = Job(b, vv;
           path=joinpath(@__DIR__, "results", "mode_I"),
           fields=(:displacement, :velocity, :acceleration, :damage))
Job{Body{BBMaterial, Peridynamics.BBPointParameters}, VelocityVerlet}(Body{BBMaterial,
Peridynamics.BBPointParameters}(BBMaterial(), 12500, [-0.49 -0.47 …
```

Note that for the keyword `fields` a NTuple is expected! If you want to export only one
field, insert `,` after the field:
```julia-repl
julia> job = Job(b, vv;
           path=joinpath(@__DIR__, "results", "mode_I"),
           fields=(:displacement,))
```

---

!!! warning "Internal use only"
    Please note that the fields are intended for internal use only. They are *not* part of
    the public API of Peridynamics.jl, and thus can be altered (or removed) at any time
    without it being considered a breaking change.

```julia
Job{SpatialSetup,TimeSolver,Options}
```

# Type Parameters

- `SpatialSetup <: AbstractSpatialSetup`: Type of the spatial setup
- `TimeSolver <: AbstractTimeSolver`: Type of the time solver
- `Options <: AbstractJobOptions`: Type of the export options

# Fields

- `spatial_setup::AbstractSpatialSetup`: Body or Multibody setup for the simulation
- `time_solver::AbstractTimeSolver`: Method for calculating discrete time steps
- `options::AbstractJobOptions`: Options for simulation data export
"""
struct Job{S<:AbstractSpatialSetup,T<:AbstractTimeSolver,O<:AbstractJobOptions}
    spatial_setup::S
    time_solver::T
    options::O

    function Job(spatial_setup::S, time_solver::T, options::O) where {S,T,O}
        pre_submission_check(spatial_setup)
        return new{S,T,O}(spatial_setup, time_solver, options)
    end
end

function Job(spatial_setup::S, time_solver::T; kwargs...) where {S,T}
    o = Dict{Symbol,Any}(kwargs)
    check_kwargs(o, JOB_KWARGS)
    options = get_export_options(storage_type(spatial_setup, time_solver), o)
    return Job(spatial_setup, time_solver, options)
end

const DEFAULT_EXPORT_FIELDS = (:displacement, :damage)

struct JobOptions{F} <: AbstractJobOptions
    export_allowed::Bool
    root::String
    vtk::String
    logfile::String
    freq::Int
    fields::F
end

function JobOptions(root::String, freq::Int, fields)
    vtk = joinpath(root, "vtk")
    logfile = joinpath(root, "logfile.log")
    return JobOptions(true, root, vtk, logfile, freq, fields)
end

function JobOptions()
    return JobOptions(false, "", "", "", 0, Vector{Symbol}())
end

function get_export_options(::Type{S}, o::Dict{Symbol,Any}) where {S<:AbstractStorage}
    local root::String
    local freq::Int

    if haskey(o, :path) && haskey(o, :freq)
        root = normpath(string(o[:path]))
        freq = Int(o[:freq])
    elseif haskey(o, :path) && !haskey(o, :freq)
        root = normpath(string(o[:path]))
        freq = 10
    elseif !haskey(o, :path) && haskey(o, :freq)
        msg = "if `freq` is specified, the keyword `path` is also needed!\n"
        throw(ArgumentError(msg))
    else
        root = ""
        freq = 0
    end
    freq < 0 && throw(ArgumentError("`freq` should be larger than zero!\n"))

    fields = get_export_fields(S, o)

    if isempty(root)
        eo = JobOptions()
    else
        eo = JobOptions(root, freq, fields)
    end

    return eo
end
