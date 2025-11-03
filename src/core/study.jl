"""
    Study(jobcreator::Function, setups::Vector{<:NamedTuple}; root::String)

$(internal_api_warning())
$(experimental_api_warning())

A structure for managing parameter studies with multiple peridynamic simulations. The
`Study` type coordinates the execution of multiple simulation jobs with different parameter
configurations, tracks their status, and logs all results to a central logfile.

# Arguments
- `jobcreator::Function`: A function with signature `jobcreator(setup::NamedTuple, root::String)`
    that creates and returns a [`Job`](@ref) object. The function receives:
    - `setup`: A `NamedTuple` containing the parameters for one simulation
    - `root`: The root directory path for the study (to construct individual job paths)
- `setups::Vector{<:NamedTuple}`: A vector of parameter configurations. Each element must be
    a `NamedTuple` with the same field names.

# Keywords
- `root::String`: Root directory path where all simulation results and the study logfile
    will be stored. This directory and all job subdirectories will be created during
    [`submit!`](@ref).

# Fields
- `jobcreator::Function`: The job creation function
- `setups`: Vector of parameter configurations
- `jobs`: Vector of created [`Job`](@ref) objects
- `jobpaths::Vector{String}`: Paths to individual job directories (must be unique)
- `root::String`: Root directory for the study
- `logfile::String`: Path to the central study logfile
- `sim_success::Vector{Bool}`: Status flags indicating successful completion of each job

# Throws
- `ArgumentError`: If `setups` is empty
- `ArgumentError`: If setups don't have consistent field names
- `ArgumentError`: If job paths are not unique
- `ArgumentError`: If `jobcreator` fails to create a valid job

# Example
```julia
function create_job(setup::NamedTuple, root::String)
    body = Body(BBMaterial(), uniform_box(1.0, 1.0, 1.0, 0.1))
    material!(body, horizon=0.3, E=setup.E, rho=1000, Gc=100)
    velocity_ic!(body, :all_points, :x, setup.velocity)
    solver = VelocityVerlet(steps=1000)
    path = joinpath(root, "sim_E\$(setup.E)_v\$(setup.velocity)")
    return Job(body, solver; path=path, freq=10)
end

setups = [
    (; E=1e9, velocity=1.0),
    (; E=2e9, velocity=1.0),
    (; E=1e9, velocity=2.0),
]

study = Study(create_job, setups; root="my_parameter_study")
```

See also: [`submit!`](@ref), [`Job`](@ref)
"""
struct Study{F,S,J}
    jobcreator::F
    setups::S
    jobs::J
    jobpaths::Vector{String}
    root::String
    logfile::String
    sim_success::Vector{Bool}

    function Study(jobcreator::F, setups::S; root::String) where {F,S}
        check_setups(setups)
        jobs = [jobcreator(setup, root) for setup in setups]
        J = typeof(jobs)
        jobpaths = [job.options.root for job in jobs]
        check_jobpaths_unique(jobpaths)
        sim_success = fill(false, length(jobs))
        logfile = joinpath(root, "study_log.log")
        new{F,S,J}(jobcreator, setups, jobs, jobpaths, root, logfile, sim_success)
    end
end

function Base.show(io::IO, ::MIME"text/plain", @nospecialize(study::Study))
    n_jobs = length(study.jobs)
    print(io, "Study with $n_jobs jobs:\n")
    for (i, jobpath) in enumerate(study.jobpaths)
        status = study.sim_success[i] ? "✓" : "✗"
        println(io, "  [$(status)] Job: `$(jobpath)`")
    end
    return nothing
end

function check_setups(setups::Vector{<:NamedTuple})
    if isempty(setups)
        throw(ArgumentError("setups vector cannot be empty!\n"))
    end

    # Check that all setups have the same field names
    first_keys = keys(setups[1])
    for (i, setup) in enumerate(setups[2:end])
        if keys(setup) != first_keys
            msg = "All setups must have the same field names!\n"
            msg *= "  First setup has fields: $first_keys\n"
            msg *= "  Setup $(i+1) has fields: $(keys(setup))\n"
            throw(ArgumentError(msg))
        end
    end

    return nothing
end

function check_jobpaths_unique(jobpaths::Vector{String})
    jobpaths_set = Set{String}(jobpaths)
    if length(jobpaths) != length(jobpaths_set)
        throw(ArgumentError("job paths must be unique for each job in the study!\n"))
    end
    return nothing
end

"""
    submit!(study::Study; kwargs...)

$(internal_api_warning())
$(experimental_api_warning())

Submit and execute all simulation jobs in a parameter study. Jobs are run sequentially,
and each job can utilize MPI or multithreading as configured. If a job fails, the error
is logged and execution continues with the remaining jobs.

This function creates the study root directory and logfile, then submits each job using
[`submit`](@ref). All simulation results and individual job logs are stored in their
respective job directories, while a central study logfile tracks the overall progress and
status of all simulations.

# Arguments
- `study::Study`: The study containing all jobs to execute

# Keywords
- `quiet::Bool`: If `true`, suppress terminal output for individual jobs (default: `false`)
- Additional keywords are passed to [`submit`](@ref) for each job

# Behavior
- Jobs execute sequentially (one after another)
- Each job inherits MPI/threading configuration from the Julia session
- Failed jobs don't stop execution of remaining jobs
- `study.sim_success` is updated with the status of each job
- A study logfile is created at `study.logfile` with:
  - Study header with timestamp
  - Parameters and status for each simulation
  - Execution time for successful jobs
- Individual job logs are written to their respective directories

# Example
```julia
study = Study(create_job, setups; root="my_study")
submit!(study; quiet=true)

# Check which simulations succeeded
successful_indices = findall(study.sim_success)
println("Successful simulations: ", successful_indices)

# Read the study logfile
logcontent = read(study.logfile, String)
```

See also: [`Study`](@ref), [`submit`](@ref)
"""
function submit!(study::Study; kwargs...)
    # Create root directory if it doesn't exist
    if !isdir(study.root)
        mkpath(study.root)
    end

    open(study.logfile, "w+") do io
        write(io, get_logfile_head())
        write(io, peridynamics_banner(color=false))
        write(io, "\nSIMULATION STUDY LOGFILE\n\n")
    end
    for (i, job) in enumerate(study.jobs)
        success = false
        simtime = @elapsed begin
            try
                submit(job; kwargs...)
                success = true
            catch err
                # Try to log the error to the job's logfile, but if that fails
                # (e.g., invalid path), just continue
                try
                    log_it(job.options, "\nERROR: Simulation failed with error!\n")
                    log_it(job.options, sprint(showerror, err, catch_backtrace()))
                    log_it(job.options, "\n")
                catch log_err
                    # If logging fails, just continue - error will be recorded in study log
                end
            end
        end
        study.sim_success[i] = success
        open(study.logfile, "a") do io
            msg = "Simulation `$(study.jobpaths[i])`:\n"
            for (key, value) in pairs(study.setups[i])
                msg *= "  $(key): $(value)\n"
            end
            if success
                msg *= @sprintf("  status: completed ✓ (%.2f seconds)\n", simtime)
            else
                msg *= "  status: failed ✗\n"
            end
            msg *= "\n"
            write(io, msg)
        end
    end
    n_successful = count(study.sim_success)
    n_jobs = length(study.jobs)
    msg = "\nStudy completed ($n_successful / $n_jobs simulations successful)\n"
    print_log(stdout, msg)
    return nothing
end

"""
    process_each_job(f::Function, study::Study, default_result::NamedTuple)

$(internal_api_warning())
$(experimental_api_warning())

Apply a processing function to each successfully completed job in a parameter study.
This function iterates through all jobs in the study and applies the user-defined
processing function `f` to jobs that completed successfully. Failed jobs or jobs where
the processing function errors will use the `default_result` instead.

# Arguments
- `f::Function`: A processing function with signature `f(job::Job, setup::NamedTuple)`
    that returns a `NamedTuple` containing the processed results. The function receives:
    - `job`: The [`Job`](@ref) object for the simulation
    - `setup`: The parameter configuration `NamedTuple` for this job
- `study::Study`: The study containing the jobs to process
- `default_result::NamedTuple`: The default result to use when a job failed or when
    the processing function throws an error

# Returns
- `Vector{<:NamedTuple}`: A vector of results with the same length as the number of
    jobs in the study. Each element is either the result from applying `f` or the
    `default_result` for failed/errored cases.

# Behavior
- Only processes jobs where `study.sim_success[i] == true`
- Failed jobs automatically receive `default_result` (warning logged)
- If processing function `f` throws an error, that job receives `default_result` (error logged)
- Processing is sequential (one job at a time)
- All jobs in the study will have a corresponding entry in the results vector

# Example
```julia
# After running a parameter study
study = Study(create_job, setups; root="my_study")
submit!(study)

# Define a function to extract maximum displacement from results
function extract_max_displacement(job::Job, setup::NamedTuple)
    # Read results from job output directory
    results_file = joinpath(job.options.root, "results_step_0010.jld2")
    data = load_results(results_file)
    max_u = maximum(norm, eachcol(data.displacement))
    return (; E=setup.E, velocity=setup.velocity, max_displacement=max_u)
end

default = (; E=0.0, velocity=0.0, max_displacement=NaN)
results = process_each_job(extract_max_displacement, study, default)

# Filter successful results
successful_results = [r for r in results if !isnan(r.max_displacement)]
```

See also: [`Study`](@ref), [`submit!`](@ref)
"""
function process_each_job(f::F, study::Study, default_result::NamedTuple) where {F}
    results = fill(default_result, length(study.jobs))
    for (i, job) in enumerate(study.jobs)
        if study.sim_success[i]
            setup = study.setups[i]
            res = try
                f(job, setup)
            catch err
                @error "error processing job $(job.options.root)" error=err
                default_result
            end
            results[i] = res
        else
            @warn "skipping processing for failed job $(job.options.root)"
        end
    end
    return results
end
