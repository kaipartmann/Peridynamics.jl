"""
    Job{S<:SpatialSetup,T<:AbstractTimeSolver}

Job that contains all the information required for a peridynamic simulation

# Fields

- `spatial_setup<:SpatialSetup`:
- `time_solver<:AbstractTimeSolver`:
- `options::ExportOptions`:


TODO @kaipartmann + kwargs

---

Constructors:

    Job(spatial_setup<:SpatialSetup, time_solver<:AbstractTimeSolver; kwargs...)

# Arguments

- `spatial_setup<:SpatialSetup`:
- `time_solver<:AbstractTimeSolver`:

# Keywords

-

# Throws

- error if keyword is not allowed

# Example

```julia-repl
julia> b = Body(BBMaterial(), pos, vol)
julia> vv = VelocityVerlet(steps=2000)
julia> job = Job(b, vv;
           path=joinpath(@__DIR__, "results", "mode_I_2"),
           write=(:displacement, :damage))
Job{Body{BBMaterial, Peridynamics.BBPointParameters}, VelocityVerlet}(Body{BBMaterial,
Peridynamics.BBPointParameters}(BBMaterial(), 12500, [-0.49 -0.47 …
```
"""
struct Job{S<:SpatialSetup,T<:AbstractTimeSolver}
    spatial_setup::S
    time_solver::T
    options::ExportOptions

    function Job(spatial_setup::S, time_solver::T, options::ExportOptions) where {S,T}
        pre_submission_check(spatial_setup)
        return new{S,T}(spatial_setup, time_solver, options)
    end
end

function Job(spatial_setup::S, time_solver::T; kwargs...) where {S,T}
    o = Dict{Symbol,Any}(kwargs)
    check_kwargs(o, JOB_KWARGS)
    options = get_export_options(storage_type(spatial_setup, time_solver), o)
    return Job(spatial_setup, time_solver, options)
end
