"""
    system_type(mat)
    system_type(mat, ::Type{FT})
    system_type(mat, ::Type{FT}, ::Val{N})

$(extension_api_note())

Return the system type of a material, instantiated for the float type `FT` and the number of
spatial dimensions `N` of the simulation. The returned instantiation holds the arrays of the
CPU. A material inherits the method from the material family it subtypes, so only a new
system needs one.

# Default

An [`InterfaceError`](@ref) that names the material family a material can subtype instead.

# Example

```julia
function Peridynamics.system_type(mat::AbstractMySystemMaterial,
                                  ::Type{FT}=Peridynamics.default_float_type(),
                                  ::Val{N}=Val(3)) where {FT,N}
    return Peridynamics.host_system_type(MySystem, Val(N), FT)
end
```

See also [`storage_type`](@ref), [`check_system_compat`](@ref),
[`host_system_type`](@ref).
"""
function system_type(mat::AbstractMaterial, ::Type=default_float_type(), ::Val=Val(3))
    return throw(InterfaceError(mat, "system_type", system_type_hint()))
end

function system_type_hint()
    msg = "A material inherits its system from the material family it subtypes, e.g.\n"
    msg *= "        struct MyMaterial <: Peridynamics.AbstractBondSystemMaterial{NoCorrection}\n"
    msg *= "    gives a `BondSystem`. Define the method directly only for a new system:\n"
    msg *= "        Peridynamics.system_type(::MyMaterial, ::Type{FT}, ::Val{N}) where {FT,N}"
    return msg
end

"""
    check_system_compat(::Type{System}, mat)

$(extension_api_note())

Check that `mat` can be discretized with `System` and throw an `ArgumentError` that names
both if it cannot. This is what the constructor of a system calls first, so that a body with
the wrong material family fails at setup instead of somewhere inside a kernel.

A system restricts its materials by defining the pair of methods that [`BondSystem`](@ref)
does, one that accepts its material family and one that rejects everything else.

# Default

Accepts every material, so a system without a restriction needs no method.

See also [`system_type`](@ref), [`@system`](@ref).
"""
@inline function check_system_compat(::Type{<:AbstractSystem}, ::AbstractMaterial)
    return nothing
end

"""
    max_n_chunks(mat)

$(extension_api_note())

Return the largest number of chunks a body of the material `mat` may be decomposed into. A
system that cannot be split defines the method and returns `1`, e.g. a system that
transforms the whole body at once. A threaded run clamps the number of chunks to this value,
an MPI run with more ranks than the material allows is an error, because the ranks are given
from the outside and cannot be clamped away.

# Default

`typemax(Int)`, so a material places no limit and the number of chunks follows the number of
threads or MPI ranks.

# Example

```julia
Peridynamics.max_n_chunks(::MyMaterial) = 1
```

See also [`system_type`](@ref), [`first_chunk`](@ref).
"""
@inline max_n_chunks(::AbstractMaterial) = typemax(Int)

"""
    first_chunk(dh)

$(extension_api_note())

Return the first body chunk of a data handler. Under multithreading this is the first of its
chunks and under MPI it is the one chunk of the rank, so a system that exists only once per
body, see [`max_n_chunks`](@ref), reads its chunk with it in `log_system` without knowing
which of the two data handlers it was given.

See also [`max_n_chunks`](@ref).
"""
function first_chunk end

#=
The number of chunks an MPI run would create is the number of ranks, which is given from the
outside and cannot be clamped, so a material that allows fewer says so here.
=#
function check_max_n_chunks(mat::AbstractMaterial, n_chunks::Int)
    n_max = max_n_chunks(mat)
    n_chunks ≤ n_max && return nothing
    msg = "the material $(nameof(typeof(mat))) allows at most $(n_max) "
    msg *= "chunk$(n_max == 1 ? "" : "s") per body, but this run has $(n_chunks) MPI "
    msg *= "ranks!\n"
    msg *= "  Start the simulation with at most $(n_max) rank$(n_max == 1 ? "" : "s"), or "
    msg *= "run it with multithreading, where the number of chunks is clamped instead.\n"
    return throw(ArgumentError(msg))
end

function get_system(::AbstractBody{M}, ::PointDecomposition, ::Int) where {M}
    msg = "system for material $M not specified!\n"
    return error(msg)
end

function log_system(options::AbstractJobOptions, dh::AbstractDataHandler)
    log_system(system_type(dh), options, dh)
    return nothing
end

"""
    get_n_loc_points(system)

$(extension_api_note())

Return the number of *local* points of a chunk, i.e. the points this chunk is responsible for
and whose equation of motion it integrates. Always `≤ get_n_points(system)`.

A field declared with the default extent [`LocalPoints`](@ref) has this many entries.

See also [`get_n_points`](@ref), [`each_point_idx`](@ref).
"""
@inline function get_n_loc_points(system::AbstractSystem)
    return get_n_loc_points(system.chunk_handler)
end

"""
    get_n_points(system)

$(extension_api_note())

Return the number of points of a chunk, local *and* halo. Halo points are the points of other
chunks that this chunk needs in order to evaluate the force density of its own points. They
are read from the neighboring chunks and never integrated here.

A field annotated with [`@lth`](@ref) or [`@htl`](@ref) has this many entries, i.e. the extent
[`HaloPoints`](@ref).

See also [`get_n_loc_points`](@ref).
"""
@inline function get_n_points(system::AbstractSystem)
    return get_n_points(system.chunk_handler)
end

@inline function get_point_ids(system::AbstractSystem)
    return get_point_ids(system.chunk_handler)
end

@inline function get_loc_points(system::AbstractSystem)
    return get_loc_points(system.chunk_handler)
end

@inline function get_halo_points(system::AbstractSystem)
    return get_halo_points(system.chunk_handler)
end

@inline function get_hidxs_by_src(system::AbstractSystem)
    return get_hidxs_by_src(system.chunk_handler)
end

@inline function get_localizer(system::AbstractSystem)
    return get_localizer(system.chunk_handler)
end

"""
    each_point_idx(system)

$(extension_api_note())

Return an iterator over the indices of the *local* points of a chunk. This is the loop the
force density calculation of a body chunk runs over, and it is what makes a material work
under multithreading and MPI without a change: every chunk only ever iterates its own points.

# Example

```julia
for i in Peridynamics.each_point_idx(system)
    force_density_point!(storage, system, mat, paramsetup, t, Δt, i)
end
```

See also [`each_bond_idx`](@ref), [`get_n_loc_points`](@ref).
"""
@inline function each_point_idx(system::AbstractSystem)
    return each_point_idx(system.chunk_handler)
end

@inline function each_point_idx_pair(system::AbstractSystem)
    return each_point_idx_pair(system.chunk_handler)
end

@inline function get_loc_view(a::AbstractArray, system::AbstractSystem)
    return get_loc_view(a, system.chunk_handler)
end

"""
    get_n_dim(x)

$(extension_api_note())

Return the number of spatial dimensions `N` of a system, of a storage or of a nested
constitutive or damage state, read from its type parameter. A storage is built with the `N`
of its system, so the two can never disagree. Every system of this package is 3D, so this
always returns `3` today.

Use [`dims`](@ref) to turn it into the `Val{N}` that the vector, tensor and symmetric tensor
accessors of `static_arrays.jl` take, see [`get_vector`](@ref) for example.
"""
function get_n_dim end

"""
    dims(x)

$(extension_api_note())

Return the number of spatial dimensions of `x` as a `Val{N}`, which is the last argument of
every accessor of `static_arrays.jl`. This is how a kernel or a hook names its dimension:
`dims(system)` where a system is in scope, `dims(storage)` or `dims(state)` where only a
storage or a nested model state is.

`get_n_dim` reads a type parameter, so the `Val` is a compile time constant and the accessor
call folds into plain indexing. `dims(body)` is *not* a constant, because a body reads its
dimension off the size of its position matrix, so it must never appear in a kernel.

# Example

```julia
Δxij = Peridynamics.get_vector_diff(storage.position, i, j, Peridynamics.dims(system))
```

See also [`get_n_dim`](@ref), [`get_vector`](@ref).
"""
@inline dims(x) = Val(get_n_dim(x))

"""
    each_bond_idx(system, i)

$(extension_api_note())

Return an iterator over the bond indices of point `i`. A bond index addresses every bond
field of the system and of the storage, that is every field declared with a `Bond...` field
shape.

# Example

```julia
for bond_id in Peridynamics.each_bond_idx(system, i)
    j = Peridynamics.get_neighbor(system, bond_id)
    Peridynamics.bond_is_active(storage, system, bond_id) || continue
end
```

See also [`each_point_idx`](@ref), [`get_n_bonds`](@ref).
"""
@inline Base.@propagate_inbounds each_bond_idx(system::AbstractSystem, i::Int) = system.bond_ids[i]

"""
    get_n_bonds(system)

$(extension_api_note())

Return the number of bonds of a chunk. Every storage field declared with a `Bond...` field
shape has this many entries, and so does the state of a constitutive model that indexes per
bond, see [`@cm_storage`](@ref).
"""
@inline get_n_bonds(system::AbstractSystem) = length(system.neighbor)

# The reads below are `@inbounds`: a bond id always comes from `each_bond_idx` or one of its
# relatives, and the check would otherwise keep the length of the array alive inside every
# force loop, which costs the cheap kernels measurably.
"""
    get_neighbor(system, bond_id)

$(extension_api_note())

Return the point index of the neighbor of bond `bond_id`, i.e. the point at the other end of
the bond. This is the twin of `get_point` for a bond instead of a degree of freedom.

The read is unchecked, so `bond_id` must come from [`each_bond_idx`](@ref) or one of its
relatives.

See also [`reference_bond_length`](@ref), [`bond_may_fail`](@ref), [`each_bond_idx`](@ref).
"""
@inline get_neighbor(system::AbstractSystem, bond_id::Int) = @inbounds system.neighbor[bond_id]

"""
    reference_bond_length(system, bond_id)

$(extension_api_note())

Return the length of bond `bond_id` in the reference configuration. The current length in the
deformed configuration is [`current_bond_length`](@ref), and the ratio of the two minus one is
[`bond_stretch`](@ref).

The read is unchecked, so `bond_id` must come from [`each_bond_idx`](@ref) or one of its
relatives.

See also [`get_neighbor`](@ref), [`bond_may_fail`](@ref), [`each_bond_idx`](@ref).
"""
@inline function reference_bond_length(system::AbstractSystem, bond_id::Int)
    return @inbounds system.bond_length[bond_id]
end

"""
    bond_may_fail(system, bond_id)

$(extension_api_note())

Return whether bond `bond_id` is allowed to fail. It is `false` for the bonds of a point that
[`no_failure!`](@ref) protects and for a body without fracture parameters, and it is what a
failure criterion checks alongside the stretch, e.g.
`ε > εc && Peridynamics.bond_may_fail(system, bond_id)`.

The read is unchecked, so `bond_id` must come from [`each_bond_idx`](@ref) or one of its
relatives.

See also [`get_neighbor`](@ref), [`reference_bond_length`](@ref), [`each_bond_idx`](@ref).
"""
@inline bond_may_fail(system::AbstractSystem, bond_id::Int) = @inbounds system.fail_permit[bond_id]

"""
    surface_correction_factor(system, bond_id)

$(extension_api_note())

Return the surface correction factor of bond `bond_id`. A material multiplies the force
density of that bond by it. The correction of a body chunk is `system.correction`, and which
one it is follows from the `Correction` type parameter of the material, see
`AbstractBondSystemMaterial`.

# Default

`1`, both with [`NoCorrection`](@ref) and on a system without a `correction` field, e.g. a
[`BondAssociatedSystem`](@ref). A material can therefore always write the multiplication and
pays nothing when no correction is used.

# Example

```julia
ω = Peridynamics.bond_is_active(storage, system, bond_id) *
    Peridynamics.surface_correction_factor(system, bond_id)
b = ω * params.bc * ε / l .* Δxij
```

See also [`bond_is_active`](@ref), [`each_bond_idx`](@ref).
"""
@inline function surface_correction_factor(system::AbstractSystem, bond_id::Int)
    if hasfield(typeof(system), :correction)
        return surface_correction_factor(system.correction, bond_id)
    end
    return 1
end

"""
    default_float_type()

$(internal_api_warning())

Return the floating point type of all simulation data. This is the single place that decides
it, so nothing below hardcodes `Float64` any more. It will become a property of [`Body`](@ref)
and flow into the system, the point parameters and the storage. Until then every
[`float_type`](@ref) method returns this value.
"""
@inline default_float_type() = Float64

"""
    float_type(x)

$(extension_api_note())

Return the floating point type of the simulation data of `x`. Every storage field declared
with a field shape that does not pin its element type, e.g. `PointVector` instead of
`PointVector{Float64}`, is allocated with this type, see [`@storage`](@ref).
"""
function float_type end

@inline float_type(::AbstractSystem) = default_float_type()

@inline function get_n_dof(system::AbstractSystem)
    return get_n_dim(system) * get_n_points(system)
end

@inline function get_n_loc_dof(system::AbstractSystem)
    return get_n_dim(system) * get_n_loc_points(system)
end

@inline function each_dim(system::AbstractSystem)
    return 1:get_n_dim(system)
end

@inline function get_dof(system::AbstractSystem, dim::Int, point_id::Int)
    return get_dof(get_n_dim(system), dim, point_id)
end

@inline function get_dof(n_dim::Int, dim::Int, point_id::Int)
    return (point_id - 1) * n_dim + dim
end

@inline function each_dof_idx(system::AbstractSystem)
    return each_dof_idx(get_n_dim(system), 1:get_n_points(system))
end

@inline function each_loc_dof_idx(system::AbstractSystem)
    return each_dof_idx(get_n_dim(system), 1:get_n_loc_points(system))
end

@inline function each_dof_idx(system::AbstractSystem, idxs::AbstractVector{<:Integer})
    return each_dof_idx(get_n_dim(system), idxs)
end

#=
This function generates a cartesian product of dof indices, dimensions, and point indices.
In Julia it is very convenient, because elements in a 2-dimensional array can be
accessed via:
    A[dof]
or via
    A[dim, i]
=#
@inline function each_dof_idx(n_dim::Int, idxs::AbstractVector{<:Integer})
    return ((get_dof(n_dim, dim, i), dim, i) for i in idxs, dim in 1:n_dim)
end

@inline function each_dof(system::AbstractSystem)
    return each_dof(get_n_dim(system), 1:get_n_points(system))
end

@inline function each_loc_dof(system::AbstractSystem)
    return each_dof(get_n_dim(system), 1:get_n_loc_points(system))
end

@inline function each_dof(system::AbstractSystem, idxs::AbstractVector{<:Integer})
    return each_dof(get_n_dim(system), idxs)
end

@inline function each_dof(n_dim::Int, idxs::AbstractVector{<:Integer})
    return (get_dof(n_dim, dim, i) for i in idxs, dim in 1:n_dim)
end

# Get the point index from a dof index.
get_point(system::AbstractSystem, idx::Int) = get_point(get_n_dim(system), idx)
get_point(n_dim::Int, idx::Int) = div(idx - 1, n_dim) + 1

# Get the dimension index from a dof index.
get_dim(system::AbstractSystem, idx::Int) = get_dim(get_n_dim(system), idx)
get_dim(n_dim::Int, idx::Int) = mod(idx - 1, n_dim) + 1

@inline function init_field_system(system, field)
    return nothing
end

# the current position of every local and halo point starts at the reference position
@inline function init_field_system(system::AbstractSystem, ::Val{:position})
    return copy(system.position)
end

function log_material(mat::M; indentation::Int=2) where {M}
    msg = msg_qty("material type", nameof(M); indentation)
    for prop in fieldnames(M)
        msg *= log_material_property(Val(prop), mat; indentation)
    end
    return msg
end

function log_material_property(::Val{S}, mat; indentation) where {S}
    return ""
end

# --------------------------------------------------------------------------------------
# the @system macro
# --------------------------------------------------------------------------------------

#=
The system twin of `core/storages.jl`.

A system is the discretization of a body chunk. It owns the same kind of shaped fields a
storage owns, so it is declared the same way: `@system` reads the field declarations of
`core/storage_fields.jl`, derives the type parameters from them and generates the struct,
the positional constructor, `Adapt.adapt_structure`, `host_system_type` and the accessors
that answer for the shape layer.

The one difference to `@storage` is that a system may declare type parameters of its own.
They come first in the generated header, so that a pattern like
`BondSystem{<:EnergySurfaceCorrection}` keeps dispatching after `N`, `FT` and the array
parameters were appended behind them.
=#

"""
    host_type(T, ::Val{N}, ::Type{FT})

$(internal_api_warning())

Return the concrete instantiation of `T` whose arrays live on the CPU. This is what
[`host_system_type`](@ref) fills a user type parameter of a system with, e.g. the
correction of a [`BondSystem`](@ref), naming the dimension and the float type of the
system along the way.

Anything that holds no arrays is its own host type, which is the default.
"""
function host_type end

host_type(::Type{T}, ::Val, ::Type) where {T} = T

"""
    host_system_type(::Type{System}, user_params..., ::Val{N}, ::Type{FT})

$(extension_api_note())

Return the instantiation of a system generated by [`@system`](@ref) whose arrays live on the
CPU, for `N` spatial dimensions and the float type `FT`. This is what
[`system_type`](@ref) returns, so that the body chunks of a data handler are concrete.

The `user_params` are the type parameters the system declares itself, in the order it
declares them, each of them already a host type.

# Default

The method is generated by [`@system`](@ref), so a system declared with that macro needs
nothing here.

See also [`@system`](@ref), [`system_type`](@ref).
"""
function host_system_type end

"""
    SystemSizes{N,FT}(chunk_handler, n_bonds)

$(extension_api_note())

The sizes a system has before the system exists. A constructor allocates the fields of its
system against this instead of writing `zeros(3, n_points)` by hand, so a shaped field of a
system is allocated by exactly the same [`alloc_field`](@ref) the storage side uses. It
answers [`get_n_dim`](@ref), [`float_type`](@ref), the point counts through the chunk
handler and [`get_n_bonds`](@ref), and nothing else.

# Example

```julia
sizes = SystemSizes{N,FT}(chunk_handler, length(neighbor))
kernels = alloc_field(BondScalar(), sizes, LocalPoints())
```

See also [`@system`](@ref), [`get_n_dim`](@ref), [`float_type`](@ref).
"""
struct SystemSizes{N,FT,CH<:AbstractChunkHandler} <: AbstractSystem
    chunk_handler::CH
    n_bonds::Int
end

function SystemSizes{N,FT}(chunk_handler::CH, n_bonds::Int) where {N,FT,CH}
    return SystemSizes{N,FT,CH}(chunk_handler, n_bonds)
end

@inline get_n_dim(::SystemSizes{N}) where {N} = N
@inline float_type(::SystemSizes{N,FT}) where {N,FT} = FT
@inline get_n_bonds(s::SystemSizes) = getfield(s, :n_bonds)

# --------------------------------------------------------------------------------------
# type parameters of the generated system
# --------------------------------------------------------------------------------------

"""
    SystemTypeParam

$(internal_api_warning())

One derived type parameter of a system generated by [`@system`](@ref): its `name`, the
expression it is `bound` by in the struct header and the `default` it has in the CPU
instantiation that [`host_system_type`](@ref) returns.
"""
struct SystemTypeParam
    name::Symbol
    bound::Any
    default::Any
end

#=
The chunk handler is not declared by the author of a system, it is injected by the macro
itself as the last field of the struct, so its type parameter has a fixed name and a fixed
position in the header, right after `FT` and before the derived array parameters.
=#
const CH_TYPE_PARAM = :CH

"""
    derive_system_type_params(decls, user_params, name)

$(internal_api_warning())

Derive the type parameters of a system from its field declarations. A field declared with a
field shape or with a concrete `Array` type contributes the same parameter [`@storage`](@ref)
derives for it, one per distinct pair of element type and number of array dimensions.
Everything else keeps the type it is declared with, except an abstract type, which is an
error: a system field needs a concrete type, or a user type parameter of the system bounded
by that abstract type. The parameters are returned in the order the fields first ask for
them.
"""
function derive_system_type_params(decls, user_params, name::AbstractString)
    params = Vector{SystemTypeParam}()
    array_params = Vector{StorageTypeParam}()
    field_params = Dict{Symbol,Symbol}()
    taken = Set{Symbol}(user_params)
    push!(taken, DIM_TYPE_PARAM)
    push!(taken, FLOAT_TYPE_PARAM)
    push!(taken, CH_TYPE_PARAM)
    for decl in decls
        key = storage_field_param_key(decl)
        if !isnothing(key)
            T, n_dims = key
            idx = findfirst(p -> p.eltype === T && p.n_dims == n_dims, array_params)
            if isnothing(idx)
                pname = unique_param_name(storage_param_name(array_params, T, n_dims), taken)
                push!(taken, pname)
                array_param = StorageTypeParam(pname, T, n_dims)
                push!(array_params, array_param)
                push!(params,
                      SystemTypeParam(pname, storage_param_bound(array_param),
                                      storage_param_default(array_param)))
                idx = lastindex(array_params)
            end
            field_params[decl.name] = array_params[idx].name
            continue
        end
        T = decl.type
        (isa(T, Base.Type) && isabstracttype(T)) || continue
        msg = "the system field `$(decl.name)` of `$(name)` is declared with the abstract "
        msg *= "type `$(T)`!\n"
        msg *= "  A system field needs a concrete type. Either give `$(decl.name)` a "
        msg *= "concrete type, or declare a type parameter of the system bounded by "
        msg *= "`$(T)` and use that for the field, e.g. `P<:$(T)` and `$(decl.name)::P`.\n"
        throw(ArgumentError(msg))
    end
    return (; params, field_params)
end

function unique_param_name(name::Symbol, taken)
    in(name, taken) || return name
    n = 2
    while in(Symbol(name, "_", n), taken)
        n += 1
    end
    return Symbol(name, "_", n)
end

# --------------------------------------------------------------------------------------
# the macros
# --------------------------------------------------------------------------------------

"""
    @system struct MySystem ... end
    @system struct MySystem{P} <: AbstractSystem ... end

$(extension_api_note())

Define a system, the discretization of a body chunk. This is [`@storage`](@ref) for a
system, and the body accepts the same field declarations. Unlike a storage, a system has no
reusable field blocks. `@storage_fields` and `@inherit` are for a storage, and writing
`@inherit` inside `@system` is an error, so a system lists every one of its fields directly.

# Declaring a field

A field is declared with a **field shape**, which says what the field *is* and from which the
container type, the element type and the number of entries follow:

| shape | container | rows | entries |
|:---|:---|:---|:---|
| [`PointScalar`](@ref) | `Vector` | – | points |
| [`PointVector`](@ref) | `Matrix` | `N` | points |
| [`BondScalar`](@ref) | `Vector` | – | bonds |
| [`BondVector`](@ref) | `Matrix` | `N` | bonds |

and so on for every shape of [`@storage`](@ref). A field can also be declared with a plain
concrete type, e.g. `two_nis::Vector{TwoNeighborInteraction}`. A field declared with an
abstract type is an error, because a system field needs a concrete type: give it one, or add
a type parameter of the system bounded by that abstract type and declare the field with it.

Unlike a storage field, a system field never specifies an initial value with `= value`: the
constructor of the system fills every field, so an initial value in the declaration would be
written and never read.

A system never declares a `chunk_handler` field itself, doing so is an error. The macro
injects it as the last field of the struct: it is what answers the point counts and the halo
bookkeeping, see [`get_n_loc_points`](@ref) and [`each_point_idx`](@ref).

# The generated type

The type parameters of the system are the ones it declares itself, then `N` for the number of
spatial dimensions, then `FT` for the float type of the simulation, then `CH` for the type of
the injected chunk handler, then one parameter per distinct array type of its fields. So

```julia
Peridynamics.@system struct BondSystem{Correction<:AbstractCorrection} <: AbstractBondSystem
    position::PointVector{Float64}
    volume::PointScalar
    neighbor::BondScalar{Int}
    bond_length::BondScalar
    fail_permit::BondScalar{Bool}
    n_neighbors::PointScalar{Int}
    bond_ids::PointScalar{UnitRange{Int}}
    kernels::BondScalar
    correction::Correction
end
```

becomes a struct whose header starts with `BondSystem{Correction,N,FT,CH,...}`. The declared
parameters come first, so `BondSystem{<:EnergySurfaceCorrection}` still selects a method,
and every dispatch on a system has to be written with `<:` for that reason.

`N` and `FT` are what [`get_n_dim`](@ref) and [`float_type`](@ref) read off the system, and
a storage is built with the `N` of its system, so the two can never disagree. The generated
constructor checks that `position` really has `N` rows. A user type parameter named `N`,
`FT` or `CH` is an error, because the macro fills those in itself.

The macro generates the positional constructor `MySystem{user...,N,FT}(fields...,
chunk_handler)`, which takes the declared fields in the order they were declared and the
chunk handler last, and infers `CH` and the array parameters from the values given.
It also generates `Adapt.adapt_structure`, [`host_system_type`](@ref),
[`get_n_dim`](@ref), [`float_type`](@ref) and, for a system that declares a bond field,
[`get_n_bonds`](@ref). What a system defines itself is [`system_type`](@ref), its constructor
`MySystem(body, pd, chunk_id)` and `get_system`.

# Example

```julia
Peridynamics.@system struct MySystem <: Peridynamics.AbstractSystem
    position::PointVector{Float64}
    volume::PointScalar
end
```

See also [`@storage`](@ref), [`SystemSizes`](@ref), [`system_type`](@ref),
[`check_system_compat`](@ref), [`host_system_type`](@ref).
"""
macro system(system)
    macrocheck_input_system_struct(system)
    return __system(system, __module__)
end

function __system(system, mod::Module)
    local _data = get_system_structdef(system, mod)
    local _system_type = _data.system_type
    local _user_params = _data.user_params
    local _user_param_names = _data.user_param_names
    local _decls = _data.decls
    local _params = _data.params
    local _field_params = _data.field_params

    local _struct = _data.system_struct
    local _constructor = system_constructor_expr(_system_type, _user_params,
                                                 _user_param_names, _decls, _params,
                                                 _field_params)
    local _adapt_structure = system_adapt_expr(_system_type, _user_params,
                                               _user_param_names, _decls)
    local _host_system_type = host_system_type_expr(_system_type, _user_params,
                                                    _user_param_names, _params)
    local _accessors = system_accessor_exprs(_system_type, _user_param_names, _decls)

    # `chunk_handler` is not one of `_decls`, the macro injects it itself, but the field
    # table has to show it, so a synthetic declaration is appended for display only
    local _display_decls = vcat(_decls,
                                StorageFieldDecl(:chunk_handler, :none, CH_TYPE_PARAM,
                                                 nothing, nothing))
    local _storage_fields_expr = quote
        function Peridynamics.storage_fields_expr(::Base.Type{<:$(esc(_system_type))})
            return $(QuoteNode(_display_decls))
        end
    end

    # the docstring is attached last, so that a `$(block_table(Name))` in it can already read
    # the declarations that were just registered
    local _doc = quote
        Base.@__doc__ $(esc(_system_type))
    end

    return Expr(:block, _struct, _constructor, _adapt_structure, _host_system_type,
                _accessors..., _storage_fields_expr, _doc)
end

function get_system_structdef(system_expr, mod::Module)
    system_type, user_params, supertype = get_system_header(system_expr)
    user_param_names = [user_param_name(p) for p in user_params]
    name = string(system_type)
    check_user_param_names(user_param_names, name)
    check_no_system_inherit(system_expr.args[3].args, name)
    decls = get_field_decls(system_expr.args[3].args, mod)
    check_system_decls(decls, name)
    (; params, field_params) = derive_system_type_params(decls, user_param_names, name)
    fields_exprs = [Expr(:(::), d.name, storage_field_decl_type(d, field_params))
                    for d in decls]
    push!(fields_exprs, Expr(:(::), :chunk_handler, CH_TYPE_PARAM))
    fields_block = Expr(:block, fields_exprs...)
    header = system_header_expr(system_type, supertype, user_params, params)
    system_struct = Expr(:struct, system_expr.args[1], header, fields_block)
    return (; system_struct, system_type, user_params, user_param_names, decls, params,
            field_params)
end

"""
    get_system_header(system_expr)

$(internal_api_warning())

Return the name, the declared type parameters and the supertype of the struct of a
[`@system`](@ref) definition. Unlike a storage, a system may declare type parameters of its
own, and they stay the leading parameters of the generated type.
"""
function get_system_header(system_expr)
    header = system_expr.args[2]
    supertype = :(Peridynamics.AbstractSystem)
    if header isa Expr && header.head === :(<:)
        supertype = header.args[2]
        header = header.args[1]
    end
    if header isa Symbol
        return header, Any[], supertype
    elseif header isa Expr && header.head === :curly && header.args[1] isa Symbol
        return header.args[1], Any[header.args[2:end]...], supertype
    end
    msg = "system header `$(header)` not supported!\n"
    msg *= "  Write a system as `struct MySystem <: AbstractSystem`, with its own type "
    msg *= "parameters if it needs them, e.g. "
    msg *= "`struct MySystem{P<:AbstractCorrection} <: AbstractSystem`.\n"
    return throw(ArgumentError(msg))
end

function user_param_name(param)
    param isa Symbol && return param
    if param isa Expr && param.head === :(<:) && param.args[1] isa Symbol
        return param.args[1]
    end
    msg = "type parameter `$(param)` of a `@system` definition not supported!\n"
    msg *= "  A system declares its type parameters as `P` or as `P<:Bound`.\n"
    return throw(ArgumentError(msg))
end

function system_header_expr(system_type, supertype, user_params, params)
    param_exprs = Any[user_params...]
    push!(param_exprs, DIM_TYPE_PARAM)
    push!(param_exprs, Expr(:(<:), FLOAT_TYPE_PARAM, Real))
    push!(param_exprs, Expr(:(<:), CH_TYPE_PARAM, :(Peridynamics.AbstractChunkHandler)))
    for param in params
        push!(param_exprs, Expr(:(<:), param.name, param.bound))
    end
    return Expr(:(<:), Expr(:curly, system_type, param_exprs...), supertype)
end

#=
The where clause of every generated method: the declared parameters with their bounds, then
`N`, `FT` and `CH`, then the derived parameters. The order matters, because a bound may name
a parameter that was declared before it, e.g. `V_FT<:AbstractVector{FT}`.
=#
function system_where_params(user_params, params)
    where_params = Any[user_params...]
    push!(where_params, DIM_TYPE_PARAM)
    push!(where_params, Expr(:(<:), FLOAT_TYPE_PARAM, Real))
    push!(where_params, Expr(:(<:), CH_TYPE_PARAM, :(Peridynamics.AbstractChunkHandler)))
    for param in params
        push!(where_params, Expr(:(<:), param.name, param.bound))
    end
    return where_params
end

function system_curly_params(user_param_names, params)
    curly_params = Any[user_param_names...]
    push!(curly_params, DIM_TYPE_PARAM)
    push!(curly_params, FLOAT_TYPE_PARAM)
    push!(curly_params, CH_TYPE_PARAM)
    for param in params
        push!(curly_params, param.name)
    end
    return curly_params
end

"""
    system_constructor_expr(system_type, user_params, user_param_names, decls, params,
                            field_params)

$(internal_api_warning())

Assemble the positional constructor `MySystem{user...,N,FT}(fields...,chunk_handler)` of a
generated system: the declared fields in the order they were declared, then the injected
chunk handler last. Neither `N` nor `FT` follows from a field type, so the constructor
names both and lets the field types and the chunk handler answer every other parameter. It
is what the constructor of the system and `Adapt.adapt_structure` build through.
"""
function system_constructor_expr(system_type, user_params, user_param_names, decls, params,
                                 field_params)
    where_params = system_where_params(user_params, params)
    curly_params = system_curly_params(user_param_names, params)
    head_params = Any[user_param_names...]
    push!(head_params, DIM_TYPE_PARAM)
    push!(head_params, FLOAT_TYPE_PARAM)
    args = [Expr(:(::), d.name, storage_field_decl_type(d, field_params)) for d in decls]
    push!(args, Expr(:(::), :chunk_handler, CH_TYPE_PARAM))
    names = [d.name for d in decls]
    push!(names, :chunk_handler)
    call = Expr(:call, Expr(:curly, esc(system_type), head_params...), args...)
    signature = Expr(:where, call, where_params...)
    body = Any[]
    if any(d -> d.name === :position && field_n_dims_of(d) == 2, decls)
        check = :(Peridynamics.check_system_n_dim(position,
                                                  Base.Val($(DIM_TYPE_PARAM))))
        push!(body, check)
    end
    push!(body,
          Expr(:return,
               Expr(:call, Expr(:curly, esc(system_type), curly_params...), names...)))
    return Expr(:function, signature, Expr(:block, body...))
end

@inline field_n_dims_of(decl::StorageFieldDecl) = _field_n_dims_of(decl, decl.shape)
@inline _field_n_dims_of(decl, shape::AbstractFieldShape) = field_n_dims(shape)
function _field_n_dims_of(decl, ::Nothing)
    T = decl.type
    (isa(T, Base.Type) && T <: Base.Array && isconcretetype(T)) || return 0
    return ndims(T)
end

"""
    check_system_n_dim(position, ::Val{N})

$(internal_api_warning())

Check that the `position` matrix of a system really has `N` rows, so that the number of
spatial dimensions of a system and its positions can never disagree.
"""
@inline function check_system_n_dim(position::AbstractMatrix, ::Val{N}) where {N}
    Base.size(position, 1) === N && return nothing
    msg = "the system is built for $(N) spatial dimensions, but its `position` has "
    msg *= "$(Base.size(position, 1)) rows!\n"
    return throw(DimensionMismatch(msg))
end

#=
Moving a system to another array backend is a matter of moving every field. A field whose
type is a declared parameter of the system, e.g. `correction::Correction`, may change its
type on the way, so the parameter is re-read from the moved value; `N` and `FT` are carried
over unchanged, because no field answers them.
=#
function system_adapt_expr(system_type, user_params, user_param_names, decls)
    head_params = Any[user_param_names...]
    push!(head_params, DIM_TYPE_PARAM)
    push!(head_params, FLOAT_TYPE_PARAM)
    where_params = Any[user_params...]
    push!(where_params, DIM_TYPE_PARAM)
    push!(where_params, FLOAT_TYPE_PARAM)
    assignments = [:($(d.name) = Adapt.adapt(to, Base.getfield(s, $(QuoteNode(d.name)))))
                   for d in decls]
    push!(assignments, :(chunk_handler = Adapt.adapt(to, Base.getfield(s, :chunk_handler))))
    # a declared parameter that is the type of exactly one field follows that field
    curly_params = Any[]
    for name in user_param_names
        idx = findfirst(d -> d.type === name, decls)
        push!(curly_params, isnothing(idx) ? name : :(Base.typeof($(decls[idx].name))))
    end
    push!(curly_params, DIM_TYPE_PARAM)
    push!(curly_params, FLOAT_TYPE_PARAM)
    call = Expr(:call, :(Adapt.adapt_structure), :to,
                Expr(:(::), :s, Expr(:curly, esc(system_type), head_params...)))
    signature = Expr(:where, call, where_params...)
    body = Expr(:block, assignments...,
                Expr(:return,
                     Expr(:call, Expr(:curly, esc(system_type), curly_params...),
                          (d.name for d in decls)..., :chunk_handler)))
    return Expr(:function, signature, body)
end

function host_system_type_expr(system_type, user_params, user_param_names, params)
    args = Any[Expr(:(::), Expr(:curly, :(Base.Type), esc(system_type)))]
    for name in user_param_names
        push!(args, Expr(:(::), Expr(:curly, :(Base.Type), name)))
    end
    push!(args, Expr(:(::), Expr(:curly, :(Base.Val), DIM_TYPE_PARAM)))
    push!(args,
          Expr(:kw, Expr(:(::), Expr(:curly, :(Base.Type), FLOAT_TYPE_PARAM)),
               :(Peridynamics.default_float_type())))
    where_params = Any[user_params...]
    push!(where_params, DIM_TYPE_PARAM)
    push!(where_params, FLOAT_TYPE_PARAM)
    curly_params = Any[user_param_names...]
    push!(curly_params, DIM_TYPE_PARAM)
    push!(curly_params, FLOAT_TYPE_PARAM)
    push!(curly_params, :(Peridynamics.ChunkHandler))
    for param in params
        push!(curly_params, param.default)
    end
    call = Expr(:call, :(Peridynamics.host_system_type), args...)
    signature = Expr(:where, call, where_params...)
    body = Expr(:block,
                Expr(:return, Expr(:curly, esc(system_type), curly_params...)))
    return Expr(:function, signature, body)
end

#=
`get_n_dim` and `float_type` read the parameters of the system, so both fold to a constant in
a kernel. `get_n_bonds` exists only for a system that really has bonds, and reads the first
bond-shaped field it declares.
=#
function system_accessor_exprs(system_type, user_param_names, decls)
    dim_pattern = Expr(:curly, esc(system_type), user_param_names..., DIM_TYPE_PARAM)
    dim_where = Any[user_param_names..., DIM_TYPE_PARAM]
    float_pattern = Expr(:curly, esc(system_type), user_param_names..., DIM_TYPE_PARAM,
                         FLOAT_TYPE_PARAM)
    float_where = Any[user_param_names..., DIM_TYPE_PARAM, FLOAT_TYPE_PARAM]
    exprs = Any[]
    dim_body = Expr(:block, Expr(:return, DIM_TYPE_PARAM))
    float_body = Expr(:block, Expr(:return, FLOAT_TYPE_PARAM))
    for (fn, arg, where_params, fn_body) in
        ((:(Peridynamics.get_n_dim), dim_pattern, dim_where, dim_body),
         (:(Peridynamics.get_n_dim), Expr(:curly, :(Base.Type), Expr(:(<:), dim_pattern)),
          dim_where, dim_body),
         (:(Peridynamics.float_type), float_pattern, float_where, float_body))
        call = Expr(:call, fn, Expr(:(::), arg))
        method = Expr(:function, Expr(:where, call, where_params...), fn_body)
        push!(exprs, Expr(:macrocall, Symbol("@inline"), nothing, method))
    end
    bond_decl = findfirst(d -> isa(d.shape, AbstractBondFieldShape), decls)
    if !isnothing(bond_decl)
        decl = decls[bond_decl]
        n_bonds = if field_n_dims(decl.shape) == 1
            :(Base.length(Base.getfield(s, $(QuoteNode(decl.name)))))
        else
            :(Base.size(Base.getfield(s, $(QuoteNode(decl.name))), 2))
        end
        push!(exprs, quote
                  @inline function Peridynamics.get_n_bonds(s::$(esc(system_type)))
                      return $(n_bonds)
                  end
              end)
    end
    return exprs
end

# --------------------------------------------------------------------------------------
# what a system may not declare
# --------------------------------------------------------------------------------------

function check_system_decls(decls, name::AbstractString)
    for decl in decls
        if decl.name === :chunk_handler
            msg = "the system field `chunk_handler` of `$(name)` is declared explicitly!\n"
            msg *= "  The `@system` macro provides this field itself, as the last field of "
            msg *= "the struct and the type parameter `CH<:AbstractChunkHandler`. Remove "
            msg *= "the declaration.\n"
            throw(ArgumentError(msg))
        end
        if !isnothing(decl.init)
            msg = "the system field `$(decl.name)` of `$(name)` specifies the initial "
            msg *= "value `$(decl.init)`!\n"
            msg *= "  A system field never has one: the constructor of the system fills "
            msg *= "every field, so an initial value would be written and never read.\n"
            throw(ArgumentError(msg))
        end
        if is_halo_decl(decl)
            msg = "the system field `$(decl.name)` of `$(name)` is annotated with "
            msg *= "`@$(decl.annotation)`!\n"
            msg *= "  A system is built once during setup and never changes, so nothing of "
            msg *= "it is exchanged between chunks. Only storage fields are annotated.\n"
            throw(ArgumentError(msg))
        end
        if is_nested_state_decl(decl)
            marker = nested_state_marker(is_cm_state_decl(decl) ? :cm : :dmg)
            msg = "the system field `$(decl.name)` of `$(name)` is declared with "
            msg *= "`$(marker)`!\n"
            msg *= "  A system carries no model state, that is what the storage of the "
            msg *= "material is for.\n"
            throw(ArgumentError(msg))
        end
    end
    return nothing
end

#=
Reusable field blocks (`@storage_fields` and `@inherit`) exist only for a storage. A system
lists every one of its fields directly, so `@inherit` inside `@system` is rejected before the
fields are even parsed.
=#
function check_no_system_inherit(block_args, name::AbstractString)
    for field in block_args
        field isa LineNumberNode && continue
        macro_name = (field isa Expr && field.head === :macrocall) ?
                     get_macro_name(field.args[1]) : nothing
        macro_name === Symbol("@inherit") || continue
        msg = "the system `$(name)` uses `@inherit`!\n"
        msg *= "  Reusable field blocks exist only for storages, declared with "
        msg *= "`@storage_fields` and included with `@inherit`. A system lists every one "
        msg *= "of its fields directly.\n"
        throw(ArgumentError(msg))
    end
    return nothing
end

#=
`N`, `FT` and `CH` are filled in by the macro itself, so a user type parameter with one of
those names would collide with the parameter the macro generates for it.
=#
function check_user_param_names(user_param_names, name::AbstractString)
    for pname in user_param_names
        pname in (DIM_TYPE_PARAM, FLOAT_TYPE_PARAM, CH_TYPE_PARAM) || continue
        msg = "the type parameter `$(pname)` of `$(name)` collides with the parameter "
        msg *= "`@system` reserves for it!\n"
        msg *= "  `N`, `FT` and `CH` are filled in by the macro, rename the type "
        msg *= "parameter.\n"
        throw(ArgumentError(msg))
    end
    return nothing
end

function macrocheck_input_system_struct(system)
    (system isa Expr && system.head === :struct) && return nothing
    (system isa Expr && system.head === :escape) && return nothing
    msg = "specified input is not a valid system struct expression!\n"
    return throw(ArgumentError(msg))
end
