const QUIET = Ref(false)

@inline quiet() = QUIET[]
@inline set_quiet!(b::Bool) = (QUIET[] = b; return nothing)

is_logging(io) = isa(io, Base.TTY) == false || (get(ENV, "CI", nothing) == "true")

@inline progress_enabled() = !is_logging(stderr) && !quiet()

function init_logs(options::ExportOptions)
    options.exportflag || return nothing
    mpi_isroot() || return nothing
    mkpath(options.vtk)
    init_logfile(options)
    return nothing
end

function init_logfile(options::ExportOptions)
    open(options.logfile, "w+") do io
        write(io, "--- LOGFILE ---\n")
        if mpi_run()
            write(io, @sprintf("MPI SIMULATION WITH %d RANKS\n", mpi_nranks()))
        else
            write(io, @sprintf("MULTITHREADING SIMULATION WITH %d THREADS\n", nthreads()))
        end
        write(io, get_git_info())
    end
    return nothing
end

function get_git_info()
    repo = LibGit2.GitRepo(pkgdir(@__MODULE__))
    head_name = LibGit2.headname(repo)
    head_oid = LibGit2.head_oid(repo)
    dirty = LibGit2.isdirty(repo)
    tag_list = LibGit2.tag_list(repo)
    info = "GIT-INFO: "
    if head_name == "main"
        info *= "commit: $head_oid"
    else
        info *= "branch: $head_name, commit: $head_oid"
    end
    for tag in tag_list
        tag_hash = LibGit2.target(LibGit2.GitObject(repo, tag))
        if tag_hash == head_oid
            info *= ", tag: $tag"
            break
        end
    end
    if dirty
        info *= ", local changes detected"
    end
    info *= "\n"
    return info
end
