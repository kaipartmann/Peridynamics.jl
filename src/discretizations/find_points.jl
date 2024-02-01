function get_function_argument(f::F) where {F<:Function}
    func_method = get_method_of_function(f)
    argnames = get_argument_names_of_function(func_method)
    if length(argnames) > 1
        msg =  "too many arguments in specified function $(nameof(f))!\n\n"
        msg *= "Should only contain one argument with name x, y, z or p!\n"
        throw(ArgumentError(msg))
    end
    arg = first(argnames)
    if !in(arg, FIND_POINTS_ALLOWED_SYMBOLS)
        msg =  "unknown argument name in specified function $(nameof(f))!\n\n"
        msg *= "Should only contain argument with name x, y, z or p!\n"
        throw(ArgumentError(msg))
    end
    return arg
end

function find_points(f::F, position::Matrix{Float64}) where {F<:Function}
    func_method = get_method_of_function(f)
    argnames = get_argument_names_of_function(func_method)
    if isempty(argnames)
        msg =  "no argument in specified function $(nameof(f))!\n\n"
        msg *= "The filter function should contain one argument with name x, y, z or p!\n"
        throw(ArgumentError(msg))
    end
    if length(argnames) > 1
        msg =  "too many arguments in specified function $(nameof(f))!\n\n"
        msg *= "The filter function should contain one argument with name x, y, z or p!\n"
        throw(ArgumentError(msg))
    end
    arg = first(argnames)
    if !in(arg, FIND_POINTS_ALLOWED_SYMBOLS)
        msg =  "unknown argument name in specified function $(nameof(f))!\n\n"
        msg *= "The filter function should contain one argument with name x, y, z or p!\n"
        throw(ArgumentError(msg))
    end
    arg === :p && return findall(f, eachcol(position))
    dim = SYMBOL_TO_DIM[arg]
    return findall(f, @view(position[dim, :]))
end
