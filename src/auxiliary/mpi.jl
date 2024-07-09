const MPI_INITIALIZED = Ref(false)
const MPI_RUN = Ref(false)
const MPI_RUN_FORCED = Ref(false)
const MPI_ISROOT = Ref(false)
const MPI_PROGRESS_BARS = Ref(false)
const TO = TimerOutput()

@inline mpi_comm() = MPI.COMM_WORLD
@inline mpi_rank() = MPI.Comm_rank(MPI.COMM_WORLD)
@inline mpi_nranks() = MPI.Comm_size(MPI.COMM_WORLD)
@inline mpi_run() = MPI_RUN[]
@inline mpi_chunk_id() = mpi_rank() + 1

"""
    mpi_isroot()

Helper function that returns a bool indicating if a process is the MPI root process.
It can be safely used even for multithreading simulations, as it is always `true` if the
package is started in a normal Julia environment which is not started by MPI.
"""
@inline mpi_isroot() = MPI_ISROOT[]
@inline mpi_progress_bars() = MPI_PROGRESS_BARS[]
@inline set_mpi_progress_bars!(b::Bool) = (MPI_PROGRESS_BARS[] = b; return nothing)

"""
    force_mpi_run!()

Helper function to force the usage of the MPI backend.
After this function is called, all following simulations will use MPI.
"""
@inline function force_mpi_run!()
    MPI_RUN[] = true
    MPI_RUN_FORCED[] = true
    return nothing
end

"""
    force_threads_run!()

Helper function to force the usage of the multithreading backend.
After this function is called, all following simulations will use multithreading.
"""
@inline function force_threads_run!()
    MPI_RUN[] = false
    MPI_RUN_FORCED[] = true
    return nothing
end

@inline function set_mpi_run!(b::Bool)
    MPI_RUN_FORCED[] && return nothing
    MPI_RUN[] = b
    return nothing
end

function init_mpi()
    if !MPI_INITIALIZED[]
        MPI.Init(finalize_atexit=true)
        MPI_INITIALIZED[] = true
    end
    MPI_ISROOT[] = mpi_rank() == 0
    set_mpi_run!(mpi_run_initial_check())
    return nothing
end

function mpi_run_initial_check()
    MPI_RUN_FORCED[] && return MPI_RUN[]
    haskey(ENV, "MPI_LOCALRANKID") && return true
    mpi_nranks() > 1 && return true
    nthreads() > 1 && return false
    for env_key in keys(ENV)
        if contains(env_key, "SLURM")
            msg = "will fall back to multithreading, due to "
            msg *= "`nranks` = `nthreads` = 1\n"
            msg *= "If you want to run MPI with only 1 rank, use `force_mpi_run!()`"
            @warn msg
            return false
        end
    end
    return false
end

function log_mpi_timers(options::AbstractJobOptions)
    options.export_allowed || return nothing
    file = joinpath(options.root, @sprintf("timers_rank_%d.log", mpi_rank()))
    open(file, "w+") do io
        show(IOContext(io, :displaysize => (24,150)), TO)
        write(io, "\n")
    end
    return nothing
end

"""
    enable_mpi_timers!()

Helper function to enable timers defined with the `TimerOutputs` for simulations with the
MPI backend. The results of the timers then will be exported into the specified path of a
[`Job`](@ref). By default, not timers will be used with MPI simulations. It can be safely
used with multithreading.

See also [`disable_mpi_timers!`](@ref).
"""
function enable_mpi_timers!()
    Core.eval(Peridynamics, :(TimerOutputs.enable_debug_timings(Peridynamics)))
    return nothing
end

"""
    disable_mpi_timers!()

Helper function to disable timers defined with the `TimerOutputs` for simulations with the
MPI backend. It is mainly used to reset the behaviour after a call of
[`enable_mpi_timers!`](@ref). It can be safely used with multithreading.
"""
function disable_mpi_timers!()
    Core.eval(Peridynamics, :(TimerOutputs.disable_debug_timings(Peridynamics)))
    return nothing
end

"""
    enable_mpi_progress_bars!()

Helper function to enable progress bars with MPI simulations on a personal computer.
After this function is called, progress bars are beeing shown with MPI simulations like with
multithreading simulations. This behavior can be reset to default with
[`reset_mpi_progress_bars!`](@ref).

!!! warning "Progress bars and output files"
    Progress bars are by default disabled with MPI simulations, because they can really mess
    up with the output files produced by a HPC system.
    Therefore, a warning is shown as a reminder to reset this behaviour before
    submitting a job to a cluster!
"""
function enable_mpi_progress_bars!()
    mpi_run() || return nothing
    set_mpi_progress_bars!(true)
    return nothing
end

"""
    reset_mpi_progress_bars!()

After this function is called, progress bars are again disabled on MPI simulations
(standard setting). This will reset the behavior after a call of
[`enable_mpi_progress_bars!`](@ref).
"""
function reset_mpi_progress_bars!()
    mpi_run() || return nothing
    set_mpi_progress_bars!(false)
    return nothing
end

"""
    @mpitime expression

Time the `expression` if the mpi rank is zero. Lowers to:
```julia
if mpi_isroot()
    @time expression
else
    expression
end
```

See also: [`mpi_isroot`](@ref).
"""
macro mpitime(expr)
    return quote
        if Peridynamics.mpi_isroot()
            Base.@time $(esc(expr))
        else
            $(esc(expr))
        end
    end
end

"""
    @mpiroot [option] expression

Run the code if the mpi rank is zero. Lowers to something similar as:
```julia
if mpi_isroot()
    expression
end
```

# Options
- `:wait`: All MPI ranks will wait until the root rank finishes evaluating `expression`.

See also: [`mpi_isroot`](@ref).
"""
macro mpiroot end

macro mpiroot(expr)
    return quote
        if Peridynamics.mpi_isroot()
            $(esc(expr))
        end
    end
end

macro mpiroot(option, expr)
    if !(option isa QuoteNode) || option.value !== :wait
        msg = "argument `$option` is not a valid option input!\n"
        throw(ArgumentError(msg))
    end
    return quote
        if Peridynamics.mpi_isroot()
            $(esc(expr))
        end
        Peridynamics.mpi_run() && MPI.Barrier(mpi_comm())
    end
end
