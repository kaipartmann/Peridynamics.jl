"""
    Study(jobcreator::Function, setups::Vector{<:NamedTuple}; kwargs...)

A structure for managing parameter studies with multiple peridynamic simulations. The
`Study` type coordinates the execution of multiple simulation jobs with different parameter
configurations, tracks their status, and logs all results to a central logfile.

If a logfile already exists at the study root (from a previous interrupted run), the
constructor automatically reads it and initializes the `sim_success` field based on the
recorded completion status. This enables seamless resumption of interrupted studies.

# Arguments
- `jobcreator::Function`: A function with signature
    `jobcreator(setup::NamedTuple, root::String)`
    that creates and returns a [`Job`](@ref) object. The function receives:
    - `setup`: A `NamedTuple` containing the parameters for one simulation
    - `root`: The root directory path for the study (to construct individual job paths)
- `setups::Vector{<:NamedTuple}`: A vector of parameter configurations. Each element must be
    a `NamedTuple` with the same field names.

# Keywords
- `root::String`: Root directory path where all simulation results and the study logfile
    will be stored. This directory and all job subdirectories will be created during
    [`submit!`](@ref).
- `logfile_name::String`: Name of the study logfile. The file will be created in the `root`
    directory.\\
    (default: `"study_log.log"`)

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
# some function that creates a Job from a parameter setup
function create_job(setup::NamedTuple, root::String)
    body = Body(BBMaterial(), uniform_box(1.0, 1.0, 1.0, 0.1))
    material!(body, horizon=0.3, E=setup.E, rho=1000, Gc=100)
    velocity_ic!(body, :all_points, :x, setup.velocity)
    solver = VelocityVerlet(steps=1000)
    # create a unique path for this job based on parameters
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

    function Study(jobcreator::F, setups::S; root::AbstractString,
                   logfile_name::AbstractString="study_log.log") where {F,S}
        check_setups(setups)
        jobs = [jobcreator(setup, root) for setup in setups]
        J = typeof(jobs)
        jobpaths = [job.options.root for job in jobs]
        check_jobpaths_unique(jobpaths)
        sim_success = fill(false, length(jobs))
        logfile = joinpath(root, logfile_name)
        study = new{F,S,J}(jobcreator, setups, jobs, jobpaths, root, logfile, sim_success)
        # If a logfile already exists from a previous run, initialize sim_success
        # from the logfile so processing or resuming works across interrupted runs.
        isfile(logfile) && update_sim_success_from_log!(study)
        return study
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

Submit and execute all simulation jobs in a parameter study. Jobs are run sequentially,
and each job can utilize MPI or multithreading as configured. If a job fails, the error
is logged and execution continues with the remaining jobs.

This function creates the study root directory and logfile, then submits each job using
[`submit`](@ref). All simulation results and individual job logs are stored in their
respective job directories, while a central study logfile tracks the overall progress and
status of all simulations.

## Resuming Interrupted Runs

If the study logfile already exists (from a previous run), `submit!` will:
- Read the logfile and refresh `study.sim_success` to detect previously completed jobs
- Skip jobs that already completed successfully (marked with "skipped" in the logfile)
- Only execute jobs that failed or were not yet started
- Append a "RESUMED" marker with timestamp to the logfile

This allows you to safely resume a study after an interruption (crash, timeout, manual
stop) without re-running successful simulations. Simply call `submit!(study)` again with
the same study configuration.

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

# If the process was interrupted, simply call submit! again to resume:
# It will skip completed jobs and only run remaining ones
study = Study(create_job, setups; root="my_study")  # Reads existing logfile
submit!(study; quiet=true)  # Resumes from where it left off

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
    isdir(study.root) || mkpath(study.root)

    # If logfile exists, refresh sim_success and append a resume marker. Otherwise
    # create a new logfile with header.
    if isfile(study.logfile)
        update_sim_success_from_log!(study)
        @mpiroot :wait open(study.logfile, "a") do io
            datetime = Dates.format(Dates.now(), "yyyy-mm-dd, HH:MM:SS")
            write(io, "\n\n--- RESUMED: $datetime ---\n\n")
        end
    else
        @mpiroot :wait open(study.logfile, "w+") do io
            write(io, get_logfile_head())
            write(io, peridynamics_banner(color=false))
            write(io, "\nSIMULATION STUDY LOGFILE\n\n")
        end
    end

    # Now submit each job
    n_jobs = length(study.jobs)
    for (i, job) in enumerate(study.jobs)
        # If this job already completed in a previous run, skip execution
        @mpiroot :wait open(study.logfile, "a") do io
            msg = @sprintf "(%d/%d) Simulation `%s`:\n" i n_jobs job.options.root
            for (key, value) in pairs(study.setups[i])
                msg *= "  $(key): $(value)\n"
            end
            write(io, msg)
        end
        if study.sim_success[i]
            # Write a short note about skipping to the study logfile
            @mpiroot :wait open(study.logfile, "a") do io
                write(io, "  status: skipped (completed in a previous run)\n\n")
            end
            continue
        end
        simtime = @elapsed begin
            success = try
                # remove the vtk files from previous failed runs if any
                @mpiroot :wait if isdir(job.options.vtk)
                    rm(job.options.vtk; force=true, recursive=true)
                end
                submit(job; kwargs...)
                true
            catch err
                @mpiroot :wait if !(quiet())
                    printstyled(stderr, "\nERROR:"; color=:red, bold=true)
                    msg = " job `$(job.options.root)` failed with error!\n"
                    msg *= sprint(showerror, err, catch_backtrace())
                    println(stderr, msg)
                end
                false
            end
        end
        study.sim_success[i] = success
        @mpiroot :wait open(study.logfile, "a") do io
            if success
                msg = @sprintf("  status: completed ✓ (%.2f seconds)\n\n", simtime)
            else
                msg = "  status: failed ✗\n\n"
            end
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
    update_sim_success_from_log!(study::Study)

$(internal_api_warning())

Read the `study.logfile` and update `study.sim_success` flags according to the
last recorded status for each job. This allows resuming processing or submission
after an interrupted run.
"""
function update_sim_success_from_log!(study::Study)
    isfile(study.logfile) || return nothing

    # Read logfile line by line
    lines = readlines(study.logfile)

    for (i, path) in enumerate(study.jobpaths)
        # Search for the Simulation line for this job path
        sim_line_pattern = "Simulation `$(path)`:"
        sim_line_idxs = findall(line -> occursin(sim_line_pattern, line), lines)
        if isempty(sim_line_idxs)
            # No record for this job in logfile
            study.sim_success[i] = false
            continue
        end
        sim_line_idx = sim_line_idxs[end]  # Take the last occurrence (restarts also logged)

        # Look for status line in the next few lines after the simulation line
        found_status = false
        for j in (sim_line_idx + 1):length(lines)
            if occursin("status: completed", lines[j])
                study.sim_success[i] = true
                found_status = true
                break
            elseif occursin("status: failed", lines[j])
                study.sim_success[i] = false
                found_status = true
                break
            elseif occursin("status: skipped", lines[j])
                study.sim_success[i] = true
                found_status = true
                break
            elseif occursin("Simulation `", lines[j])
                # Reached next simulation entry without finding status
                break
            end
        end

        if !found_status
            # No status found after simulation line
            study.sim_success[i] = false
        end
    end
    return nothing
end

"""
    process_each_job(f::Function, study::Study, default_result::NamedTuple)

Apply a processing function to each successfully completed job in a parameter study.
This function iterates through all jobs in the study and applies the user-defined
processing function `f` to jobs that completed successfully. Failed jobs or jobs where
the processing function errors will use the `default_result` instead.

Before processing, this function automatically refreshes `study.sim_success` from the
logfile if it exists. This ensures that jobs completed in a previous interrupted run
are properly detected and processed, even if the current Julia session doesn't reflect
their completion status yet.

# Arguments
- `f::Function`: A processing function with signature `f(job::Job, setup::NamedTuple)`
    that returns a `NamedTuple` containing the processed results. The function receives:
    - `job`: The [`Job`](@ref) object for the simulation
    - `setup`: The parameter configuration `NamedTuple` for this job
- `study::Study`: The study containing the jobs to process
- `default_result::NamedTuple`: The default result to use when a job failed or when
    the processing function throws an error

# Keywords
- `only_root::Bool`: If `true`, the processing function `f` runs only on the MPI root rank.
    Non-root ranks will use `default_result` for all jobs. This is useful when processing
    involves file I/O or expensive computations that should not be duplicated across ranks.
    Note that the returned results vector will contain meaningful data only on the root rank
    unless you manually broadcast results afterward. (default: `false`)

# Returns
- `Vector{<:NamedTuple}`: A vector of results with the same length as the number of
    jobs in the study. Each element is either the result from applying `f` or the
    `default_result` for failed/errored cases.

# Behavior
- Refreshes job completion status from logfile (if exists) before processing
- Only processes jobs where `study.sim_success[i] == true`
- Failed jobs automatically receive `default_result` (warning logged)
- If processing function `f` throws an error, that job receives `default_result` (error
    logged with MPI barrier synchronization)
- Processing is sequential (one job at a time)
- All jobs in the study will have a corresponding entry in the results vector

# MPI Behavior
- By default (`only_root=false`), the processing function `f` is called on **all MPI ranks**
- With `only_root=true`, `f` runs only on root rank; non-root ranks use `default_result`
- Error handling includes automatic MPI barrier synchronization across all ranks

!!! danger "Processing on MPI only at the root process"
    When using `only_root=true`, ensure that the processing function `f` does not contain
    any MPI calls or operations that require synchronization across ranks, as only the root
    process will execute it. This envolves also the `barrier` keyword in the
    [`process_each_export`](@ref) function if used within `f`, which should be always set to
    `false` in this case to avoid deadlocks!

# Examples
```julia
# After running a parameter study (possibly interrupted and resumed)
study = Study(create_job, setups; root="my_study")
# Status is automatically refreshed from logfile in constructor

# Define a function to extract maximum displacement from results
function extract_max_displacement(job::Job, setup::NamedTuple)
    # Process files and calculate maximum displacement
    process_each_export(job; serial=true) do r0, r, id
        # Read and process data
    end
    max_u = ...
    return (; E=setup.E, velocity=setup.velocity, max_displacement=max_u)
end

default = (; E=0.0, velocity=0.0, max_displacement=NaN)
# Use only_root=true since we're doing file I/O
results = process_each_job(extract_max_displacement, study, default; only_root=true)

# Filter successful results (only root rank will have meaningful data)
if mpi_isroot()
    successful_results = [r for r in results if !isnan(r.max_displacement)]
end
```

See also: [`Study`](@ref), [`submit!`](@ref), [`process_each_export`](@ref),
[`mpi_isroot`](@ref)
"""
function process_each_job(f::F, study::Study, default_result::NamedTuple;
                         only_root::Bool=false) where {F}
    # Refresh sim_success from logfile if it exists (in case study was interrupted/resumed)
    isfile(study.logfile) && update_sim_success_from_log!(study)

    results = fill(default_result, length(study.jobs))
    for (i, job) in enumerate(study.jobs)
        if study.sim_success[i]
            setup = study.setups[i]
            res = try
                # Only execute on root if only_root=true
                if only_root && !mpi_isroot()
                    default_result
                else
                    f(job, setup)
                end
            catch err
                @mpiroot :wait begin
                    msg = "job `$(job.options.root)` failed with error!\n"
                    msg *= sprint(showerror, err, catch_backtrace())
                    # print error to stderr
                    printstyled(stderr, "\nERROR: "; color=:red, bold=true)
                    println(stderr, msg)
                    # also write error to a separate logfile in the job directory
                    t = Dates.format(Dates.now(), "yyyy-mm-dd_HH-MM-SS")
                    err_logfile = joinpath(job.options.root, "$(t)_proc_error.log")
                    open(err_logfile, "w+") do io
                        write(io, get_logfile_head())
                        write(io, peridynamics_banner(color=false))
                        write(io, "ERROR: " * msg)
                    end
                end
                default_result
            end
            results[i] = res
        else
            @warn "skipping processing for failed job $(job.options.root)"
        end
    end
    return results
end
