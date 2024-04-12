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

required_fields_conditions() = (:velocity, :velocity_half, :b_ext)

function req_storage_fields_conditions(::Type{Storage}) where {Storage}
    parameters = fieldnames(Storage)
    for req_field in required_fields_conditions()
        if !in(req_field, parameters)
            msg = "required field $req_field not found in $(Storage)!\n"
            msg *= "The field is required for the boundary conditions!\n"
            error(msg)
        end
    end
    return nothing
end
