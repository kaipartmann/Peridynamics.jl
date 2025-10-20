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
    if mpi_progress_bars()
        set_progress_bars!(true)
        msg = "progress bar settings overwritten manually!\n"
        msg *= "The use of progress bars with MPI can lead to a mess in output files!"
        @warn msg
    else
        set_progress_bars!(false)
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

function msg_path(descr_raw::AbstractString, path_raw::AbstractString;
                  linewidth::Union{Int,Nothing}=nothing,
                  leftwidth::Union{Int,Nothing}=nothing, indentation::Int=2,
                  filler::Char='.', separator::AbstractString="",
                  delimiter::AbstractString=" ",
                  continuation_label::AbstractString="(continued)")
    descr, path = strip(descr_raw), strip(path_raw)

    # Determine effective linewidth
    effective_linewidth = something(linewidth, leftwidth !== nothing ? default_linewidth() : nothing, default_linewidth())

    # Try to fit on one line
    len_filling = get_len_filling(linewidth, leftwidth, indentation, descr, separator, path)
    required_length = indentation + length(descr) + length(separator) + max(len_filling, 1) + length(path)

    if required_length <= effective_linewidth
        # Fits on one line - use msg_qty
        return msg_qty(descr, path; linewidth=linewidth, leftwidth=leftwidth,
                      indentation=indentation, filler=filler, separator=separator,
                      delimiter=delimiter, newline=true)
    end

    # Path is too long - split into multiple lines
    return _msg_path_multiline(descr, path, effective_linewidth, indentation, filler,
                               separator, delimiter, continuation_label)
end

function _msg_path_multiline(descr, path, linewidth, indentation, filler, separator, delimiter, continuation_label)
    min_filler = 2 * length(delimiter)

    # First line: calculate budget for path
    first_line_budget = linewidth - indentation - length(descr) - length(separator) - min_filler
    break_idx = _find_path_break(path, first_line_budget)

    path_first = break_idx > 0 ? path[1:break_idx] : ""
    path_rest = break_idx > 0 ? path[break_idx+1:end] : path

    # Build first line
    first_line = if !isempty(path_first)
        len_fill = linewidth - indentation - length(descr) - length(separator) - length(path_first)
        _msg_qty(descr, path_first, len_fill, indentation, filler, separator, delimiter, true)
    else
        len_fill = linewidth - indentation - length(descr) - length(separator)
        _msg_qty(descr, "", len_fill, indentation, filler, separator, delimiter, true)
    end

    # Build continuation lines
    continuation_lines = _build_continuation_lines(path_rest, continuation_label, linewidth,
                                                   indentation, filler, delimiter)

    return first_line * continuation_lines
end

function _build_continuation_lines(remaining_path, label, linewidth, indentation, filler, delimiter)
    isempty(remaining_path) && return ""

    result = ""
    min_filler = 2 * length(delimiter)
    prefix_len = indentation + length(label)

    while !isempty(remaining_path)
        budget = linewidth - prefix_len - min_filler

        if length(remaining_path) <= budget
            # Last piece fits
            len_fill = linewidth - prefix_len - length(remaining_path)
            result *= _msg_qty(label, remaining_path, len_fill, indentation, filler, "", delimiter, true)
            break
        else
            # Need another split
            break_idx = _find_path_break(remaining_path, budget)
            piece = break_idx > 0 ? remaining_path[1:break_idx] : remaining_path[1:min(budget, end)]
            remaining_path = length(piece) < length(remaining_path) ? remaining_path[length(piece)+1:end] : ""

            len_fill = linewidth - prefix_len - length(piece)
            result *= _msg_qty(label, piece, len_fill, indentation, filler, "", delimiter, true)
        end
    end

    return result
end

function _find_path_break(path::AbstractString, max_length::Int)
    length(path) <= max_length && return length(path)
    max_length <= 0 && return 0

    # Find last directory separator within budget
    last_sep = findlast(c -> c == '/' || c == '\\', @view path[1:min(max_length, end)])

    return something(last_sep, min(max_length, length(path)))
end

function msg_vec(descr_raw::AbstractString, vec::AbstractVector;
                 linewidth::Union{Int,Nothing}=nothing,
                 leftwidth::Union{Int,Nothing}=nothing, indentation::Int=2,
                 filler::Char='.', separator::AbstractString="",
                 delimiter::AbstractString=" ",
                 continuation_label::AbstractString="(continued)",
                 vec_delimiter::AbstractString=", ",
                 vec_brackets::Tuple{AbstractString,AbstractString}=("[", "]"))
    descr = strip(descr_raw)

    # Format vector as string
    vec_str = _format_vector(vec, vec_delimiter, vec_brackets)

    # Determine effective linewidth
    effective_linewidth = something(linewidth, leftwidth !== nothing ? default_linewidth() : nothing, default_linewidth())

    # Try to fit on one line
    len_filling = get_len_filling(linewidth, leftwidth, indentation, descr, separator, vec_str)
    required_length = indentation + length(descr) + length(separator) + max(len_filling, 1) + length(vec_str)

    if required_length <= effective_linewidth
        # Fits on one line - use msg_qty
        return msg_qty(descr, vec_str; linewidth=linewidth, leftwidth=leftwidth,
                      indentation=indentation, filler=filler, separator=separator,
                      delimiter=delimiter, newline=true)
    end

    # Vector string is too long - split into multiple lines
    return _msg_vec_multiline(descr, vec_str, effective_linewidth, indentation, filler,
                              separator, delimiter, continuation_label)
end

function _format_vector(vec::AbstractVector, vec_delimiter::AbstractString,
                       vec_brackets::Tuple{AbstractString,AbstractString})
    isempty(vec) && return vec_brackets[1] * vec_brackets[2]
    elements = _format_elements(vec)
    vec_str = vec_brackets[1] * join(elements, vec_delimiter) * vec_brackets[2]
    return vec_str
end

function _format_elements(vec::AbstractVector{<:Real})
    return [@sprintf("%.7g", x) for x in vec]
end
function _format_elements(vec::AbstractVector{T}) where {T}
    return [string(x) for x in vec]
end

function _msg_vec_multiline(descr, vec_str, linewidth, indentation, filler, separator,
                            delimiter, continuation_label)
    min_filler = 2 * length(delimiter)

    # First line: calculate budget for vector string
    first_line_budget = linewidth - indentation - length(descr) - length(separator) - min_filler
    break_idx = _find_vec_break(vec_str, first_line_budget)

    vec_first = break_idx > 0 ? vec_str[1:break_idx] : ""
    vec_rest = break_idx > 0 ? vec_str[break_idx+1:end] : vec_str

    # Build first line
    first_line = if !isempty(vec_first)
        len_fill = linewidth - indentation - length(descr) - length(separator) - length(vec_first)
        _msg_qty(descr, vec_first, len_fill, indentation, filler, separator, delimiter, true)
    else
        len_fill = linewidth - indentation - length(descr) - length(separator)
        _msg_qty(descr, "", len_fill, indentation, filler, separator, delimiter, true)
    end

    # Build continuation lines
    continuation_lines = _build_vec_continuation_lines(vec_rest, continuation_label, linewidth,
                                                       indentation, filler, delimiter)

    return first_line * continuation_lines
end

function _build_vec_continuation_lines(remaining_vec, label, linewidth, indentation, filler, delimiter)
    isempty(remaining_vec) && return ""

    result = ""
    min_filler = 2 * length(delimiter)
    prefix_len = indentation + length(label)

    while !isempty(remaining_vec)
        budget = linewidth - prefix_len - min_filler

        if length(remaining_vec) <= budget
            # Last piece fits
            len_fill = linewidth - prefix_len - length(remaining_vec)
            result *= _msg_qty(label, remaining_vec, len_fill, indentation, filler, "", delimiter, true)
            break
        else
            # Need another split
            break_idx = _find_vec_break(remaining_vec, budget)
            piece = break_idx > 0 ? remaining_vec[1:break_idx] : remaining_vec[1:min(budget, end)]
            remaining_vec = length(piece) < length(remaining_vec) ? remaining_vec[length(piece)+1:end] : ""

            len_fill = linewidth - prefix_len - length(piece)
            result *= _msg_qty(label, piece, len_fill, indentation, filler, "", delimiter, true)
        end
    end

    return result
end

function _find_vec_break(vec_str::AbstractString, max_length::Int)
    length(vec_str) <= max_length && return length(vec_str)
    max_length <= 0 && return 0

    # Find last comma or space within budget (prefer comma)
    last_comma = findlast(==(','), @view vec_str[1:min(max_length, end)])

    if last_comma !== nothing
        # Include the comma and any following spaces
        break_point = last_comma
        while break_point < min(max_length, length(vec_str)) &&
              vec_str[break_point + 1] == ' '
            break_point += 1
        end
        return break_point
    end

    # If no comma, look for space
    last_space = findlast(==(' '), @view vec_str[1:min(max_length, end)])

    return something(last_space, min(max_length, length(vec_str)))
end

function msg_fields_inline(obj::T) where {T}
    fields = fieldnames(T)
    values = Tuple(getfield(obj, field) for field in fields)
    return _msg_fields_inline(fields, values)
end

function msg_fields_inline(obj, fields)
    values = Tuple(getfield(obj, field) for field in fields)
    return _msg_fields_inline(fields, values)
end

function _msg_fields_inline(fields, values)
    n_fields = length(fields)
    msg = ""
    for i in eachindex(fields, values)
        field, value = fields[i], values[i]
        msg *= msg_qty(field, value; linewidth=0, indentation=0, separator="=",
                       delimiter="", newline=false)
        if i != n_fields
            msg *= ", "
        end
    end
    return msg
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
