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

function init_logs(options::AbstractJobOptions)
    mpi_isroot() || return nothing
    print_log(peridynamics_banner())
    print_log(get_run_info())
    set_progress_bars!()
    options.export_allowed || return nothing
    mkpath(options.vtk)
    init_logfile(options)
    return nothing
end

function init_logfile(options::AbstractJobOptions)
    open(options.logfile, "w+") do io
        write(io, get_logfile_head())
        write(io, peridynamics_banner(color=false))
        write(io, get_run_info())
    end
    return nothing
end

function add_to_logfile(options::AbstractJobOptions, msg::AbstractString)
    options.export_allowed || return nothing
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
    msg *= get_version_info()
    msg *= get_git_info()
    msg *= "\n"
    return msg
end

function get_version_info()
    version_info = "VERSION: "
    version_info *= string(pkgversion(@__MODULE__))
    version_info *= "\n"
    return version_info
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

function peridynamics_banner(; color::Bool=true, indentation::Int=10)
    indent = indentation > 0 ? " "^indentation : ""
    if color
        c = Base.text_colors
        tx, d1, d2, d3, d4 = c[:normal], c[:blue], c[:red], c[:green], c[:magenta]
    else
        tx, d1, d2, d3, d4 = "", "", "", "", ""
    end
    msg = indent
    msg *= "                                                     $(d3)_$(tx)\n"
    msg *= indent
    msg *= " _____         $(d1)_$(tx)     _" * " "^29 * "$(d2)_$(d3)(_)$(d4)_$(tx)\n"
    msg *= indent
    msg *= "| ___ \\       $(d1)(_)$(tx)   | |" * " "^27 * "$(d2)(_) $(d4)(_)$(tx)\n"
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

function print_log(io::IO, msg::AbstractString)
    mpi_isroot() || return nothing
    quiet() && return nothing
    print(io, msg)
    return nothing
end

print_log(msg::AbstractString) = print_log(stdout, msg)

function log_it(options::AbstractJobOptions, msg::AbstractString)
    print_log(msg)
    add_to_logfile(options, msg)
    return nothing
end

function msg_qty(descr_raw::AbstractString, qty_raw::AbstractString;
                 linewidth::Union{Int,Nothing}=nothing,
                 leftwidth::Union{Int,Nothing}=nothing, indentation::Int=2,
                 filler::Char='.', separator::AbstractString="",
                 delimiter::AbstractString=" ", newline::Bool=true)
    descr, qty = strip(descr_raw), strip(qty_raw)
    len_filling = get_len_filling(linewidth, leftwidth, indentation, descr, separator, qty)
    msg = _msg_qty(descr, qty, len_filling, indentation, filler, separator, delimiter,
                   newline)
    return msg
end

@inline default_linewidth() = 82

function get_len_filling(::Nothing, ::Nothing, indentation::Int, descr::AbstractString,
                         separator::AbstractString, qty::AbstractString)
    n = default_linewidth() - indentation - length(descr) - length(separator) - length(qty)
    return n
end

function get_len_filling(linewidth::Int, ::Nothing, indentation::Int, descr::AbstractString,
                         separator::AbstractString, qty::AbstractString)
    return linewidth - indentation - length(descr) - length(separator) - length(qty)
end

function get_len_filling(::Nothing, leftwidth::Int, indentation::Int, descr::AbstractString,
                         separator::AbstractString, ::AbstractString)
    return leftwidth - indentation - length(descr) - length(separator)
end

function get_len_filling(::Int, ::Int, ::Int, ::AbstractString, ::AbstractString,
                         ::AbstractString)
    msg = "insufficient keywords for calculation of filling length!\n"
    msg *= "Specify either `linewidth` or `leftwidth`, not both!\n"
    return throw(ArgumentError(msg))
end

function _msg_qty(descr, qty, len_filling, indentation, filler, separator, delimiter,
                  newline)
    if len_filling > 1
        n_fillers = len_filling - 2
        filling = delimiter * filler^n_fillers * delimiter
    else
        filling = delimiter
    end
    lineend = newline ? "\n" : ""
    msg = " "^indentation * descr * separator * filling * qty * lineend
    return msg
end

function msg_qty(descr::AbstractString, qty::Real; kwargs...)
    return msg_qty(descr, @sprintf("%.7g", qty); kwargs...)
end

function msg_qty(descr::Any, qty::Any; kwargs...)
    return msg_qty(string(descr), string(qty); kwargs...)
end

function msg_fields_in_brackets(obj::T) where {T}
    fields = fieldnames(T)
    values = Tuple(getfield(obj, field) for field in fields)
    return _msg_fields_in_brackets(fields, values)
end

function msg_fields_in_brackets(obj, fields)
    values = Tuple(getfield(obj, field) for field in fields)
    return _msg_fields_in_brackets(fields, values)
end

function _msg_fields_in_brackets(fields, values)
    n_fields = length(fields)
    msg = "("
    for i in eachindex(fields, values)
        field, value = fields[i], values[i]
        msg *= msg_qty(field, value; linewidth=0, indentation=0, separator="=",
                       delimiter="", newline=false)
        if i == n_fields
            msg *= ")"
        else
            msg *= ", "
        end
    end
    return msg
end

function msg_fields(obj::T) where {T}
    fields = fieldnames(T)
    values = Tuple(getfield(obj, field) for field in fields)
    return _msg_fields(fields, values)
end

function msg_fields(obj, fields)
    values = Tuple(getfield(obj, field) for field in fields)
    return _msg_fields(fields, values)
end

function _msg_fields(fields, values)
    n_fields = length(fields)
    fields_str = string.(fields)
    leftwidth = maximum(x -> length(x), fields_str) + 4
    msg = ""
    for i in eachindex(fields, values)
        descr = fields_str[i]
        qty = values[i]
        if i == n_fields
            newline = false
        else
            newline = true
        end
        msg *= msg_qty(descr, qty; leftwidth=leftwidth, newline=newline, filler=' ')
    end
    return msg
end

function log_create_data_handler_start(io::IO=stdout)
    mpi_isroot() || return nothing
    quiet() && return nothing
    progress_bars() || return nothing
    print(io, "DATA HANDLER CREATION ... ⏳")
    return nothing
end

function log_create_data_handler_end(io::IO=stdout)
    mpi_isroot() || return nothing
    quiet() && return nothing
    progress_bars() || return nothing
    println(io, "\rDATA HANDLER CREATION COMPLETED ✔")
    return nothing
end

function log_simulation_duration(options::AbstractJobOptions, duration::Float64)
    msg = @sprintf("SIMULATION COMPLETED AFTER %g SECONDS ✔\n", duration)
    log_it(options, msg)
    return nothing
end
