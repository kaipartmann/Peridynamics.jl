function get_method_of_function(f::F) where {F<:Function}
    func_method = collect(methods(f))
    if length(func_method) > 1
        msg =  "multiple methods defined for specified function $(nameof(f))!\n\n"
        msg *= "Better use anonymous functions as an argument!\n"
        throw(ArgumentError(msg))
    end
    return first(func_method)
end

function get_argument_names_of_function(func_method::Method)
    argnames = Base.method_argnames(func_method)
    isempty(argnames) && return argnames
    first(argnames) === Symbol("#self#") && popfirst!(argnames)
    return argnames
end

function check_kwargs(p::Dict{Symbol,Any}, allowed_kwargs::NTuple{N,Symbol}) where {N}
    for key in keys(p)
        in(key, allowed_kwargs) && continue
        throw(ArgumentError(unknown_kwarg_msg(key, allowed_kwargs)))
    end
    return nothing
end

# a misspelled keyword is the first error a new user hits, so the message says what would
# have worked instead of only that this did not
function unknown_kwarg_msg(key::Symbol, allowed_kwargs::NTuple{N,Symbol}) where {N}
    msg = "keyword `$(key)` not allowed!\n"
    closest = closest_kwarg(key, allowed_kwargs)
    isnothing(closest) || (msg *= "  Did you mean `$(closest)`?\n")
    msg *= "  allowed keywords: $(join(allowed_kwargs, ", "))\n"
    return msg
end

function closest_kwarg(key::Symbol, allowed_kwargs::NTuple{N,Symbol}) where {N}
    isempty(allowed_kwargs) && return nothing
    distances = map(kwarg -> kwarg_distance(string(key), string(kwarg)), allowed_kwargs)
    distance, idx = findmin(distances)
    # a suggestion only helps if it is close; two edits on a short name is already a guess
    distance > max(1, length(string(key)) ÷ 3) && return nothing
    return allowed_kwargs[idx]
end

# Levenshtein distance, so that `epsilon_C`, `Ec` and `bta` find their keyword
function kwarg_distance(a::String, b::String)
    previous = collect(0:length(b))
    current = similar(previous)
    for (i, ca) in enumerate(a)
        current[1] = i
        for (j, cb) in enumerate(b)
            current[j + 1] = min(previous[j + 1] + 1, current[j] + 1,
                                 previous[j] + (ca == cb ? 0 : 1))
        end
        copyto!(previous, current)
    end
    return last(previous)
end
