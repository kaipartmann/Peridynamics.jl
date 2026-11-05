function system_type(mat::AbstractMaterial)
    return throw(InterfaceError(mat, "system_type", system_type_hint()))
end

function system_type_hint()
    msg = "A material inherits its system from the material family it subtypes, e.g.\n"
    msg *= "        struct MyMaterial <: Peridynamics.AbstractBondSystemMaterial{NoCorrection}\n"
    msg *= "    gives a `BondSystem`. Define the method directly only for a new system:\n"
    msg *= "        Peridynamics.system_type(::MyMaterial) = MySystem"
    return msg
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
chunks that this chunk needs in order to evaluate the force density of its own points; they
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
    force_density_point!(storage, system, mat, params, t, Δt, i)
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

@inline function get_n_dim(::AbstractSystem)
    return 3 # 3D only for now
end

"""
    default_float_type()

$(internal_api_warning())

Return the floating point type of all simulation data. This is the single place that decides
it, so nothing below hardcodes `Float64` any more. It will become a property of [`Body`](@ref)
and flow into the system, the point parameters and the storage; until then every
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
