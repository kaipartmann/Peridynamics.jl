"""
    Study

A structure for managing parameter studies with multiple peridynamic simulations.

# Fields
- `create_job::Function`: A function with signature `create_job(setup::NamedTuple)` that
    creates a [`Job`](@ref) from a setup configuration.
- `setups::Vector{NamedTuple}`: A vector of setup configurations. Each setup must be a
    `NamedTuple` with the same field names.
- `jobs::Vector{Job}`: A vector of jobs created from the setups.
- `submission_status::Vector{Bool}`: Status vector indicating whether each job was
    submitted successfully (`true`) or encountered an error (`false`).
- `postproc_status::Vector{Bool}`: Status vector indicating whether post-processing for
    each job was successful (`true`) or encountered an error (`false`).
- `results::Vector{Vector{NamedTuple}}`: Storage for results from post-processing. Each
    element corresponds to a simulation and contains a vector of `NamedTuple`s returned by
    the processing function for each time step.

# Example
```julia
function create_job(setup::NamedTuple)
    # Create body, solver, etc. using setup parameters
    body = ...
    solver = VelocityVerlet(steps=setup.n_steps)
    job = Job(body, solver; path=setup.path, freq=setup.freq)
    return job
end

setups = [
    (; n_steps=1000, path="sim1", freq=10),
    (; n_steps=2000, path="sim2", freq=20),
]

study = Study(create_job, setups)
```

See also: [`submit!`](@ref), [`postproc!`](@ref)
"""
struct Study
    create_job::Function
    setups::Vector{NamedTuple}
    jobs::Vector{Job}
    submission_status::Vector{Bool}
    postproc_status::Vector{Bool}
    results::Vector{Vector{NamedTuple}}

    function Study(create_job::Function, setups::Vector{<:NamedTuple})
        check_setups(setups)
        n = length(setups)
        jobs = Vector{Job}(undef, n)
        submission_status = fill(false, n)
        postproc_status = fill(false, n)
        results = [NamedTuple[] for _ in 1:n]

        for (i, setup) in enumerate(setups)
            try
                job = create_job(setup)
                if !(job isa Job)
                    msg = "create_job function must return a Job object!\n"
                    msg *= "For setup $i, got $(typeof(job)) instead.\n"
                    throw(ArgumentError(msg))
                end
                jobs[i] = job
            catch err
                msg = "Error creating job for setup $i:\n"
                msg *= "  Setup: $setup\n"
                msg *= "  Error: $err\n"
                throw(ArgumentError(msg))
            end
        end

        return new(create_job, setups, jobs, submission_status, postproc_status, results)
    end
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

function Base.show(io::IO, @nospecialize(study::Study))
    n_jobs = length(study.jobs)
    n_submitted = count(study.submission_status)
    n_postproc = count(study.postproc_status)
    print(io, "Study with $n_jobs simulations ($n_submitted submitted, $n_postproc post-processed)")
    return nothing
end

function Base.show(io::IO, ::MIME"text/plain", @nospecialize(study::Study))
    if get(io, :compact, false)
        show(io, study)
    else
        n_jobs = length(study.jobs)
        n_submitted = count(study.submission_status)
        n_postproc = count(study.postproc_status)
        println(io, "Study:")
        println(io, "  Number of simulations:     $n_jobs")
        println(io, "  Successfully submitted:    $n_submitted")
        println(io, "  Successfully post-processed: $n_postproc")
        if n_jobs > 0
            println(io, "  Setup parameters:          $(keys(study.setups[1]))")
        end
    end
    return nothing
end

"""
    submit!(study::Study; kwargs...)

Submit all jobs in the study for simulation. Jobs are executed **sequentially**, with each
job running to completion before the next begins. Each individual job can utilize MPI or
multithreading as configured through the [`submit`](@ref) function.

Each job is run independently, and if one job encounters an error, the remaining jobs will
still be executed. The submission status for each job is stored in `study.submission_status`.

!!! note "Execution model"
    Jobs in a study run sequentially, not in parallel. This is by design because:
    - Each simulation typically uses all available computational resources (MPI ranks/threads)
    - Running multiple MPI jobs simultaneously from one process is problematic
    - For true parallel execution of parameter studies, submit separate jobs to your HPC scheduler

# Arguments
- `study::Study`: The study containing the jobs to submit.

# Keywords
- `quiet::Bool`: If `true`, suppress output for each individual job. (default: `false`)

# Returns
- `nothing`

# Example
```julia
study = Study(create_job, setups)
submit!(study)

# Check which jobs completed successfully
successful_jobs = findall(study.submission_status)
```

See also: [`Study`](@ref), [`submit`](@ref)
"""
function submit!(study::Study; quiet::Bool=false)
    n_jobs = length(study.jobs)

    println("Starting parameter study with $n_jobs simulations...")

    for (i, job) in enumerate(study.jobs)
        println("\n" * "="^80)
        println("Simulation $i of $n_jobs")
        println("Setup: $(study.setups[i])")
        println("="^80)

        try
            submit(job; quiet=quiet)
            study.submission_status[i] = true
            println("✓ Simulation $i completed successfully")
        catch err
            study.submission_status[i] = false
            println("✗ Simulation $i encountered an error:")
            println("  Error type: $(typeof(err))")
            println("  Error message: $err")
            if err isa Exception
                println("  Stacktrace:")
                for (exc, bt) in Base.catch_stack()
                    showerror(stdout, exc, bt)
                    println()
                end
            end
            println("Continuing with remaining simulations...")
        end
    end

    println("\n" * "="^80)
    n_successful = count(study.submission_status)
    println("Parameter study completed: $n_successful of $n_jobs simulations successful")
    println("="^80)

    return nothing
end

"""
    postproc!(proc_func::Function, study::Study; kwargs...)

Apply a post-processing function to all successfully submitted jobs in the study. Results
are stored in `study.results` as a vector of vectors of `NamedTuple`s.

Post-processing can be performed in parallel (multithreaded or MPI) for each individual
simulation by setting `serial=false` (default). The parallelization happens within each
simulation's time steps, not across different simulations.

# Arguments
- `proc_func::Function`: A function with signature `proc_func(r0, r, id)` where:
    - `r0`: Reference results (from [`read_vtk`](@ref) of the initial export)
    - `r`: Current time step results (from [`read_vtk`](@ref))
    - `id`: File ID / time step number
  The function should return either a `NamedTuple` with results or `nothing`.
- `study::Study`: The study to post-process.

# Keywords
- `serial::Bool`: If `true`, process results serially on a single thread. (default: `false`)

# Returns
- `nothing`: Results are stored in `study.results`. For simulation `i`, access results via
    `study.results[i]`, which contains a vector of `NamedTuple`s (one per time step).

# Example
```julia
function proc_func(r0, r, id)
    # Calculate some quantity
    max_displacement = maximum(r[:displacement])
    return (; time_step=id, max_disp=max_displacement)
end

study = Study(create_job, setups)
submit!(study)
postproc!(proc_func, study)

# Access results for first simulation
results_sim1 = study.results[1]
```

See also: [`Study`](@ref), [`process_each_export`](@ref)
"""
function postproc!(proc_func::Function, study::Study; serial::Bool=false)
    check_process_function(proc_func)

    n_jobs = length(study.jobs)
    n_successful = count(study.submission_status)

    if n_successful == 0
        @warn "No successfully submitted jobs to post-process!"
        return nothing
    end

    println("\nStarting post-processing for $n_successful successful simulations...")

    # Clear previous results and determine if we should collect
    collect_results = Ref(false)
    first_result_type_checked = Ref(false)

    for (i, job) in enumerate(study.jobs)
        # Clear previous results for this simulation
        empty!(study.results[i])

        # Skip jobs that weren't submitted successfully
        if !study.submission_status[i]
            study.postproc_status[i] = false
            continue
        end

        println("\nPost-processing simulation $i of $n_jobs...")

        try
            function wrapper_func(r0, r, id)
                result = proc_func(r0, r, id)

                # On first result, check if we should collect
                if !first_result_type_checked[]
                    if result isa NamedTuple
                        collect_results[] = true
                    elseif !isnothing(result)
                        @warn "proc_func returned $(typeof(result)) instead of NamedTuple or nothing. Results will not be collected."
                    end
                    first_result_type_checked[] = true
                end

                # Collect if appropriate
                if collect_results[] && result isa NamedTuple
                    push!(study.results[i], result)
                end

                return nothing
            end

            process_each_export(wrapper_func, job; serial=serial)

            study.postproc_status[i] = true
            println("✓ Post-processing simulation $i completed successfully")

        catch err
            study.postproc_status[i] = false
            println("✗ Post-processing simulation $i encountered an error:")
            println("  Error type: $(typeof(err))")
            println("  Error message: $err")
            println("Continuing with remaining simulations...")
        end
    end

    println("\n" * "="^80)
    n_postproc_successful = count(study.postproc_status)
    println("Post-processing completed: $n_postproc_successful of $n_successful simulations successful")
    println("="^80)

    return nothing
end
