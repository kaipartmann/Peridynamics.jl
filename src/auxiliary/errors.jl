struct InterfaceError <: Exception
    type::DataType
    func::String
end

function InterfaceError(::T, f::F) where {T,F}
    return InterfaceError(T, string(f))
end

function Base.showerror(io::IO, e::InterfaceError)
    print(io, "interface method not correctly defined!")
    print(io, "\n  type:      ")
    printstyled(io, string(e.type); bold=true, color=:red)
    print(io, "\n  function:  ")
    printstyled(io, string(e.func); bold=true, color=:red)
    println(io)
    return nothing
end
