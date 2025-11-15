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

"""
    NaNError

$(internal_api_warning())

A type for a customized error that is thrown when NaN values are detected in the internal
force density field after the force density evaluation.

# Fields
- `time::Float64`: Simulation time when NaNs were detected.
- `step::Int`: Simulation step when NaNs were detected.
"""
struct NaNError <: Exception
    time::Float64
    step::Int
    function NaNError(t::Float64, n::Int)
        return new(t, n)
    end
end

function Base.showerror(io::IO, e::NaNError)
    print(io, "NaN values detected in simulation data!")
    print(io, "\n  time:    ")
    printstyled(io, string(e.time); bold=true, color=:red)
    print(io, "\n  step:    ")
    printstyled(io, string(e.step); bold=true, color=:red)
    println(io)
    return nothing
end
