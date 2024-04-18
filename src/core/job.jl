const JOB_KWARGS = (:path, :freq, :fields)

"""
    Job{S<:AbstractSpatialSetup,T<:AbstractTimeSolver,O<:AbstractOptions}

Job that contains all the information required for a peridynamic simulation

# Fields

- `spatial_setup<:AbstractSpatialSetup`: Body or Multibody setup for the simulation
- `time_solver<:AbstractTimeSolver`: method for calculating discrete time steps
- `options<:AbstractOptions`: options for simulation data export

---

Constructors:

    Job(spatial_setup::AbstractSpatialSetup, time_solver::AbstractTimeSolver,
        options::AbstractOptions)
    Job(spatial_setup::AbstractSpatialSetup, time_solver::AbstractTimeSolver; kwargs...)

# Arguments

- `spatial_setup<:AbstractSpatialSetup`: Body or Multibody setup for the simulation
- `time_solver<:AbstractTimeSolver`: method for calculating discrete time steps
- `options<:AbstractOptions`: options for simulation data export

# Keywords

- `path::String`: storage path for results
- `freq::Int`: frequency of time steps that are exported
- `fields`: exported fields
             Possible export fields depend on selected material model.
             See material type documentation.
             default export fields: `:displacement`, `:damage`

# Throws

- error if keyword is not allowed

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
"""
struct Job{S<:AbstractSpatialSetup,T<:AbstractTimeSolver,O<:AbstractOptions}
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

struct ExportOptions <: AbstractOptions
    exportflag::Bool
    root::String
    vtk::String
    logfile::String
    freq::Int
    fields::Vector{Symbol}
end

function ExportOptions(root::String, freq::Int, fields::Vector{Symbol})
    vtk = joinpath(root, "vtk")
    logfile = joinpath(root, "logfile.log")
    return ExportOptions(true, root, vtk, logfile, freq, fields)
end

function ExportOptions()
    return ExportOptions(false, "", "", "", 0, Vector{Symbol}())
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
        eo = ExportOptions()
    else
        eo = ExportOptions(root, freq, fields)
    end

    return eo
end
