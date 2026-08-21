"""
    InterfaceError

$(extension_api_note())

A type for a customized error that is thrown when a material model is not implemented
correctly.

# Fields

- `type::DataType`: Type that is used.
- `func::String`: Function that is used.
- `hint::String`: Optional explanation of how the interface method can be defined.
"""
struct InterfaceError <: Exception
    type::DataType
    func::String
    hint::String
    function InterfaceError(_type::T, _func::F, _hint="") where {T,F}
        func = string(_func)
        type = isa(_type, DataType) ? _type : T
        return new(type, func, string(_hint))
    end
end

function Base.showerror(io::IO, e::InterfaceError)
    print(io, "interface method not correctly defined!")
    print(io, "\n  type:    ")
    printstyled(io, string(e.type); bold=true, color=:red)
    print(io, "\n  method:  ")
    printstyled(io, string(e.func); bold=true, color=:red)
    println(io)
    isempty(e.hint) || println(io, "  ", e.hint)
    return nothing
end

"""
    StorageContractError

$(extension_api_note())

A type for a customized error that is thrown when the storage of a material does not
contain all fields that are required by the material, its damage model and the time solver
of the submitted job. See [`check_storage_contract`](@ref).

# Fields

- `storage::DataType`: Storage type that misses fields.
- `material::DataType`: Material type the storage belongs to.
- `solver::DataType`: Time solver of the job.
- `missing_fields::Vector{Pair{Symbol,String}}`: Missing fields and the reason why they are
    required.
"""
struct StorageContractError <: Exception
    storage::DataType
    material::DataType
    solver::DataType
    missing_fields::Vector{Pair{Symbol,String}}
end

function Base.showerror(io::IO, e::StorageContractError)
    print(io, "storage type does not contain all required fields!")
    print(io, "\n  storage:   ")
    printstyled(io, string(nameof(e.storage)); bold=true, color=:red)
    print(io, "\n  material:  ")
    printstyled(io, string(e.material); bold=true, color=:red)
    print(io, "\n  solver:    ")
    printstyled(io, string(e.solver); bold=true, color=:red)
    println(io, "\n\n  missing fields:")
    pad = maximum(length(string(first(x))) for x in e.missing_fields)
    for (field, reason) in e.missing_fields
        print(io, "    ")
        printstyled(io, rpad(string(field), pad); bold=true, color=:red)
        println(io, "   required by ", reason)
    end
    println(io, "\n  To fix this, either")
    println(io, "    - add the missing fields to the `Peridynamics.@storage` definition")
    println(io, "      of `", nameof(e.storage), "` and define a `Peridynamics.init_field`")
    println(io, "      method for every added field that is not already covered by the")
    println(io, "      system or the time solver, or")
    println(io, "    - use a damage model and a time solver that `", nameof(e.storage),
            "` supports.")
    return nothing
end

"""
    NaNError

$(internal_api_warning())

A type for a customized error that is thrown when `NaN` values are detected in the internal
force density field after the force density evaluation.

# Fields
- `time::Float64`: Simulation time when `NaN`s were detected.
- `step::Int`: Simulation step when `NaN`s were detected.
"""
struct NaNError <: Exception
    time::Float64
    step::Int
    function NaNError(t::Float64, n::Int)
        return new(t, n)
    end
end

function Base.showerror(io::IO, e::NaNError)
    print(io, "NaN values detected in force density field!")
    print(io, "\n  time:    ")
    printstyled(io, string(e.time); bold=true, color=:red)
    print(io, "\n  step:    ")
    printstyled(io, string(e.step); bold=true, color=:red)
    println(io)
    return nothing
end
