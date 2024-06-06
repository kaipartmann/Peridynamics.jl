const QUIET = Ref(false)
const PROGRESS_BARS = Ref(true)

@inline quiet() = QUIET[]
@inline set_quiet!(b::Bool) = (QUIET[] = b; return nothing)

@inline progress_bars() = PROGRESS_BARS[]
@inline set_progress_bars!(b::Bool) = (PROGRESS_BARS[] = b; return nothing)

function set_progress_bars!()
    is_logging = isa(stderr, Base.TTY) == false || (get(ENV, "CI", nothing) == "true")
    progress_bars_enabled = !is_logging && !quiet()
    set_progress_bars!(progress_bars_enabled)
    mpi_run() || return nothing
    if !progress_bars_enabled && !quiet() && mpi_progress_bars()
        set_progress_bars!(mpi_progress_bars())
        msg = "progress bar settings overwritten manually!\n"
        msg *= "The use of progress bars with MPI can lead to a mess in output files!"
        @warn msg
    end
    return nothing
end

function init_logs(options::AbstractOptions)
    mpi_isroot() || return nothing
    print_log(peridynamics_banner_color())
    print_log(get_run_info())
    set_progress_bars!()
    options.exportflag || return nothing
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
    options.exportflag || return nothing
    mpi_isroot() || return nothing
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

function peridynamics_banner_color(; indentation::Int=10)
    indent = indentation > 0 ? " "^indentation : ""
    msg = indent
    msg *= " _____         \u001b[31m_\u001b[0m     _                             "
    msg *= "\u001b[32m_\u001b[0m\n"
    msg *= indent
    msg *= "| ___ \\       \u001b[31m(_)\u001b[0m   | |"
    msg *= "                           \u001b[32m(_)\u001b[0m\n"
    msg *= indent
    msg *= "| |_/ /__ _ __ _  __| |_   _ _ __   __ _ _ __ ___  _  ___ ___\n"
    msg *= indent
    msg *= "|  __/ _ \\ '__| |/ _` | | | | '_ \\ / _` | '_ ` _ \\| |/ __/ __|\n"
    msg *= indent
    msg *= "| | |  __/ |  | | (_| | |_| | | | | (_| | | | | | | | (__\\__ \\\n"
    msg *= indent
    msg *= "\\_|  \\___|_|  |_|\\__,_|\\__, |_| |_|\\__,_|_| |_| |_|_|\\___|___/\n"
    msg *= indent
    msg *= "                        __/ |\n"
    msg *= indent
    msg *= "                       |___/   "
    msg *= "Copyright (c) $(Dates.format(Dates.now(), "yyyy")) Kai Partmann\n\n"
    return msg
end

function peridynamics_banner(; indentation::Int=10)
    indent = indentation > 0 ? " "^indentation : ""
    msg = indent
    msg *= " _____         _     _                             _\n"
    msg *= indent
    msg *= "| ___ \\       (_)   | |                           (_)\n"
    msg *= indent
    msg *= "| |_/ /__ _ __ _  __| |_   _ _ __   __ _ _ __ ___  _  ___ ___\n"
    msg *= indent
    msg *= "|  __/ _ \\ '__| |/ _` | | | | '_ \\ / _` | '_ ` _ \\| |/ __/ __|\n"
    msg *= indent
    msg *= "| | |  __/ |  | | (_| | |_| | | | | (_| | | | | | | | (__\\__ \\\n"
    msg *= indent
    msg *= "\\_|  \\___|_|  |_|\\__,_|\\__, |_| |_|\\__,_|_| |_| |_|_|\\___|___/\n"
    msg *= indent
    msg *= "                        __/ |\n"
    msg *= indent
    msg *= "                       |___/   "
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

function print_log(msg::AbstractString)
    mpi_isroot() || return nothing
    quiet() && return nothing
    print(msg)
    return nothing
end

function log_it(options::AbstractOptions, msg::AbstractString)
    print_log(msg)
    add_to_logfile(options, msg)
    return nothing
end

function log_qty(descr::AbstractString, qty::AbstractString; linewidth::Int=82,
                 indentation::Int=2, filler::Char='.')
    len_filling = linewidth - indentation - length(descr) - length(qty)
    if len_filling > 1
        n_fillers = len_filling - 2
        filling = " " * filler^n_fillers * " "
    else
        filling = " "
    end
    msg = " "^indentation * descr * filling * qty * "\n"
    return msg
end

function log_qty(descr::AbstractString, qty::Real; kwargs...)
    return log_qty(descr, @sprintf("%.7g", qty); kwargs...)
end

function log_qty(descr::AbstractString, qty::Any; kwargs...)
    return log_qty(descr, string(qty); kwargs...)
end

function log_create_data_handler_start()
    mpi_isroot() || return nothing
    quiet() && return nothing
    progress_bars() || return nothing
    print("DATA HANDLER CREATION ... ⏳")
    return nothing
end

function log_create_data_handler_end()
    mpi_isroot() || return nothing
    quiet() && return nothing
    progress_bars() || return nothing
    println("\rDATA HANDLER CREATION COMPLETED ✓")
    return nothing
end
