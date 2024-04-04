const FIND_POINTS_ALLOWED_SYMBOLS = (:x, :y, :z, :p)
const SYMBOL_TO_DIM = Dict(:x => 0x1, :y => 0x2, :z => 0x3)


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
