function is_duplicate(condition::C, conditions::Vector{C}) where {C<:AbstractCondition}
    for existing_condition in conditions
        override_eachother(existing_condition, condition) && return true
    end
    return false
end

function check_condition_function(f::F) where {F<:Function}
    func_method = get_method_of_function(f)
    args = get_argument_names_of_function(func_method)
    if length(args) == 1
        return :sdbc
    elseif length(args) == 2 && args[1] === :p && args[2] === :t
        return :pdsdbc
    else
        msg = "wrong arguments for position dependent condition function!\n"
        throw(ArgumentError(msg))
    end
    return nothing
end
