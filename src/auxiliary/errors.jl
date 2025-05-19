"""
    InterfaceError

$(internal_api_warning())

A type for a customized error that is thrown when a material model is not implemented
correctly.

# Fields

- `type::DataType`: Type that is used.
- `func::String`: Function that is used.
"""
struct InterfaceError <: Exception
    type::DataType
    func::String
    function InterfaceError(_type::T, _func::F) where {T,F}
        func = string(_func)
        type = isa(_type, DataType) ? _type : T
        return new(type, func)
    end
end

function Base.showerror(io::IO, e::InterfaceError)
    print(io, "interface method not correctly defined!")
    print(io, "\n  type:    ")
    printstyled(io, string(e.type); bold=true, color=:red)
    print(io, "\n  method:  ")
    printstyled(io, string(e.func); bold=true, color=:red)
    println(io)
    return nothing
end
