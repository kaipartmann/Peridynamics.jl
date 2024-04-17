const QUIET = Ref(false)

@inline quiet() = QUIET[]
@inline set_quiet!(b::Bool) = (QUIET[] = b; return nothing)

is_logging(io) = isa(io, Base.TTY) == false || (get(ENV, "CI", nothing) == "true")

@inline progress_enabled() = !is_logging(stderr) && !quiet()

function init_logs(options::AbstractOptions)
    options.exportflag || return nothing
    mpi_isroot() || return nothing
    mkpath(options.vtk)
    init_logfile(options)
    return nothing
end

function init_logfile(options::AbstractOptions)
    open(options.logfile, "w+") do io
        write(io, get_logfile_head())
        write(io, peridynamics_banner())
        write(io, get_run_info())
    end
    return nothing
end

function add_to_logfile(options::AbstractOptions, msg::AbstractString)
    open(options.logfile, "a") do io
        write(io, msg)
    end
    return nothing
end

function get_logfile_head()
    msg = "LOGFILE CREATED ON "
    msg *= Dates.format(Dates.now(), "yyyy-mm-dd, HH:MM:SS")
    msg *= "\n"
    msg *= get_git_info()
    msg *= "\n"
    return msg
end

function get_git_info()
    git_info = try
        _get_git_info()
    catch
        ""
    end
    return git_info
end

function _get_git_info()
    repo = LibGit2.GitRepo(pkgdir(@__MODULE__))
    head_name = LibGit2.headname(repo)
    head_oid = LibGit2.head_oid(repo)
    git_info = "GIT-INFO: "
    if head_name == "main"
        git_info *= "commit: $head_oid"
    else
        git_info *= "branch: $head_name, commit: $head_oid"
    end
    for tag in LibGit2.tag_list(repo)
        if LibGit2.target(LibGit2.GitObject(repo, tag)) == head_oid
            git_info *= "\n          tag: $tag"
            break
        end
    end
    if LibGit2.isdirty(repo)
        git_info *= "\n          local changes detected!"
    end
    git_info *= "\n"
    return git_info
end

function peridynamics_banner()
    msg = raw"""
    ______         _     _                             _
    | ___ \       (_)   | |                           (_)
    | |_/ /__ _ __ _  __| |_   _ _ __   __ _ _ __ ___  _  ___ ___
    |  __/ _ \ '__| |/ _` | | | | '_ \ / _` | '_ ` _ \| |/ __/ __|
    | | |  __/ |  | | (_| | |_| | | | | (_| | | | | | | | (__\__ \
    \_|  \___|_|  |_|\__,_|\__, |_| |_|\__,_|_| |_| |_|_|\___|___/
                            __/ |
                           |___/   """
    msg *= "Copyright (c) $(Dates.format(Dates.now(), "yyyy")) Kai Partmann\n\n"
    return msg
end

function get_run_info()
    msg = ""
    if mpi_run()
        msg *= @sprintf("MPI SIMULATION WITH %d RANKS\n", mpi_nranks())
    else
        msg *= @sprintf("MULTITHREADING SIMULATION WITH %d THREADS\n", nthreads())
    end
    return msg
end
