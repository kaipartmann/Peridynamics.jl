function is_duplicate(condition::C, conditions::Vector{C}) where {C<:AbstractCondition}
    for existing_condition in conditions
        override_eachother(existing_condition, condition) && return true
    end
    return false
end

function check_condition_function(f::F) where {F<:Function}
    func_method = get_method_of_function(f)
    args = get_argument_names_of_function(func_method)
    if length(args) != 1
        error("too many arguments for condition! Only one (time) is allowed!\n")
    end
    return nothing
end
