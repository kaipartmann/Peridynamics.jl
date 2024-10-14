struct InterfaceError{F} <: Exception
    type::DataType
    fun::F
end

const InterfaceSupertypes = Union{AbstractSystem,AbstractMaterial,AbstractPointParameters,
                                  AbstractStorage,AbstractTimeSolver}

function InterfaceError(::T, f::F) where {T<:InterfaceSupertypes,F<:Function}
    return InterfaceError(T, f)
end

function Base.showerror(io::IO, e::InterfaceError)
    print(io, "interface method not correctly defined!")
    print(io, "\n  type:      ")
    printstyled(io, string(e.type); bold=true, color=:red)
    print(io, "\n  function:  ")
    printstyled(io, string(e.fun); bold=true, color=:red)
    println(io)
    return nothing
end
