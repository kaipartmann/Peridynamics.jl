const JOB_KWARGS = (:path, :freq, :fields)

"""
    Job(spatial_setup, time_solver; kwargs...)

A type that contains all the information necessary for a peridynamic simulation. You can
[`submit`](@ref) a `Job` to start the simulation.

# Arguments
- `spatial_setup`: A [`Body`](@ref) or [`MultibodySetup`](@ref).
- `time_solver`: [`VelocityVerlet`](@ref) or [`DynamicRelaxation`](@ref).

# Keywords
- `path::String`: Path to store results. If it does not exist yet it will be created during
    the simulation. (optional)
- `freq::Int`: Output frequency of result files. A output file will be written every
    `freq`-th time step. (default: `10`)
- `fields`: Fields that should be exported to output files. Allowed keywords depend on the
    selected material model. Please look at the documentation of the material you specified
    when creating the body. (default: `(:displacement, :damage)`)

    If `spatial_setup` is a **`Body`**, the `fields` keyword can be of the form:
    - `fields::Symbol`: A symbol specifying a single output field.
    - `fields::NTuple{N,Symbol} where N`: A Tuple specifying multiple output fields.
    - `fields::Vector{Symbol}`: A Vector specifying multiple output fields.

    If `spatial_setup` is a **`MultibodySetup`**, the `fields` keyword can also be specified
    for every body separately:
    - `fields::Dict{Symbol,T}`: A Dictionary containing the fields separately for every
        body. `T` is here every possible type of the `fields` keyword that can be used for a
        single body.

!!! note "No file export"
    If no keyword is specified when creating a `Job`, then no files will be exported.

# Example
```julia-repl
julia> job = Job(multibody_setup, verlet_solver; path="my_results/sim1")
Job:
  spatial_setup  25880-point MultibodySetup
  time_solver    VelocityVerlet(n_steps=2000, safety_factor=0.7)
  options        export_allowed=true, freq=10
```
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
    options = get_job_options(spatial_setup, time_solver, o)
    return Job(spatial_setup, time_solver, options)
end

function Base.show(io::IO, @nospecialize(job::Job))
    n_points = Peridynamics.n_points(job.spatial_setup)
    if job.spatial_setup isa AbstractMultibodySetup
        job_descr = "-point multibody Job with "
    else
        job_descr = "-point Job with "
    end
    solver = typeof(job.time_solver)
    print(io, n_points, job_descr, solver, " solver")
    return nothing
end

function Base.show(io::IO, ::MIME"text/plain", @nospecialize(job::Job))
    if get(io, :compact, false)
        show(io, job)
    else
        println(io, "Job:")
        print(io, msg_fields(job))
    end
    return nothing
end

"""
    JobOptions{F,V}

$(internal_api_warning())

A type that contains the options of a job

# Type Parameters

- `F`: Type for fields of simulation
- `V`: Type for basename of vtk-files

# Fields

- `export_allowed::Bool`: Specify if data is exported for the job
- `root::String`: Path of the folder where all data is saved
- `vtk::String`: Path of the folder where vtk-files are saved
- `logfile::String`: Complete path of the logfile
- `freq::Int`: Frequency of exported time steps
- `fields::F`: Exported fields of the job
- `vtk_filebase::V`: Basename of exported vtk-files
"""
struct JobOptions{F,V} <: AbstractJobOptions
    export_allowed::Bool
    root::String
    vtk::String
    logfile::String
    freq::Int
    fields::F
    vtk_filebase::V
end

function Base.show(io::IO, @nospecialize(options::JobOptions))
    print(io, msg_fields_inline(options, (:export_allowed, :freq,)))
    return nothing
end

function Base.show(io::IO, ::MIME"text/plain", @nospecialize(options::JobOptions))
    if get(io, :compact, false)
        show(io, options)
    else
        println(io, "Job options:")
        print(io, msg_fields(options, (:export_allowed, :root, :freq, :fields)))
    end
    return nothing
end

function JobOptions(root::String, freq::Int, fields, vtk_filebase)
    vtk = joinpath(root, "vtk")
    logfile = joinpath(root, "logfile.log")
    return JobOptions(true, root, vtk, logfile, freq, fields, vtk_filebase)
end

function JobOptions(::AbstractBody)
    return JobOptions(false, "", "", "", 0, Vector{Symbol}(), "")
end

function JobOptions(::AbstractMultibodySetup)
    return JobOptions(false, "", "", "", 0, Dict{Symbol,Vector{Symbol}}(),
                      Dict{Symbol,String}())
end

function get_job_options(spatial_setup, solver, o)
    root, freq = get_root_and_freq(o)

    fields = get_export_fields(spatial_setup, solver, o)
    vtk_filebase = get_vtk_filebase(spatial_setup, root)

    if isempty(root)
        options = JobOptions(spatial_setup)
    else
        options = JobOptions(root, freq, fields, vtk_filebase)
    end

    return options
end

function get_root_and_freq(o::Dict{Symbol,Any})
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

    return root, freq
end
