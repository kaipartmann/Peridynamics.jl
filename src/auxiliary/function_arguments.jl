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
        if !in(key, allowed_kwargs)
            throw(ArgumentError("keyword $key not allowed!\n"))
        end
    end
    return nothing
end
