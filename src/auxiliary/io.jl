@inline default_export_fields() = [:displacement, :damage]

@inline function default_fields_spec(ms::AbstractMultibodySetup)
    fields_spec = Dict{Symbol,Vector{Symbol}}()
    for body_name in each_body_name(ms)
        fields_spec[body_name] = default_export_fields()
    end
    return fields_spec
end

function get_export_fields(body::AbstractBody, solver::AbstractTimeSolver,
                           o::Dict{Symbol,Any})
    if haskey(o, :fields)
        fields = extract_export_fields(o[:fields])
    else
        fields = default_export_fields()
    end

    check_export_fields(body, solver, fields)

    return fields
end

function get_export_fields(ms::AbstractMultibodySetup, solver::AbstractTimeSolver,
                           o::Dict{Symbol,Any})
    if haskey(o, :fields)
        fields_spec = extract_fields_spec(ms, o[:fields])
    else
        fields_spec = default_fields_spec(ms)
    end

    check_export_fields(ms, solver, fields_spec)

    return fields_spec
end

@inline function extract_export_fields(o)
    fields::Vector{Symbol} = _extract_export_fields(o)
    return fields
end

@inline _extract_export_fields(fields::Symbol) = [fields]
@inline _extract_export_fields(fields::NTuple{N,Symbol}) where {N} = [f for f in fields]
@inline _extract_export_fields(fields::Vector{Symbol}) = fields
@inline function _extract_export_fields(arbitrary_fields::T) where {T}
    fields::Vector{Symbol} = try
        convert(Vector{Symbol}, arbitrary_fields)
    catch
        msg = "wrong type of specification for keyword `fields`!\n"
        msg *= "Cannot convert input of type `$(T)` to `Vector{Symbol}`!\n"
        throw(ArgumentError(msg))
    end
    return fields
end

@inline function extract_fields_spec(ms::AbstractMultibodySetup, o)
    fields_spec::Dict{Symbol,Vector{Symbol}} = _extract_fields_spec(ms, o)
    for body_name in each_body_name(ms)
        if !haskey(fields_spec, body_name)
            fields_spec[body_name] = default_export_fields()
        end
    end
    return fields_spec
end

function _extract_fields_spec(ms::AbstractMultibodySetup, o::Dict{Symbol,T}) where {T}
    fields_spec = Dict{Symbol,Vector{Symbol}}()
    for (body_name, fields) in o
        check_if_bodyname_is_defined(ms, body_name)
        fields_spec[body_name] = extract_export_fields(fields)
    end
    return fields_spec
end

function _extract_fields_spec(ms::AbstractMultibodySetup, o)
    fields = extract_export_fields(o)
    fields_spec = Dict{Symbol,Vector{Symbol}}()
    for body_name in each_body_name(ms)
        fields_spec[body_name] = fields
    end
    return fields_spec
end

function check_export_fields(body::AbstractBody, solver::AbstractTimeSolver,
                             fields::Vector{Symbol})
    S = storage_type(body)
    check_export_fields(S, fields)
    return nothing
end

function check_export_fields(ms::AbstractMultibodySetup, solver::AbstractTimeSolver,
                             fields_spec::Dict{Symbol,Vector{Symbol}})
    for body_name in each_body_name(ms)
        body = get_body(ms, body_name)
        S = storage_type(body)
        fields = fields_spec[body_name]
        check_export_fields(S, fields)
    end
    return nothing
end

function check_export_fields(::Type{S}, fields::Vector{Symbol}) where {S}
    allowed_fields = point_data_fields(S)
    for f in fields
        if !in(f, allowed_fields) && !custom_field(S, f)
            msg = "unknown point data field `:$(f)` specified for export!\n"
            msg *= "If you intend to export a custom field, please define a method:"
            msg *= "\n   `Peridynamics.custom_field(::Type{<:$(nameof(S))}, ::Val{:$(f)}) = true`\n"
            msg *= "Otherwise, see here all available point data fields of $(nameof(S)):\n"
            for allowed_name in allowed_fields
                msg *= "  - $allowed_name\n"
            end
            throw(ArgumentError(msg))
        end
    end
    return nothing
end

"""
    custom_field(::Type{<:MyStorage}, ::Val{field})

$(extension_api_note())

Return `true` for every field name that a storage exports but that is not one of its own
field names, e.g. a quantity that [`export_field`](@ref) derives on the fly. Without this,
naming the field in `Job(...; fields=(...,))` is rejected as a typo.

Defaults to `false`, so a field that *is* a field of the storage needs nothing.

# Example

```julia
Peridynamics.custom_field(::Type{<:MyStorage}, ::Val{:von_mises_stress}) = true
```

See also [`export_field`](@ref).
"""
custom_field(S::Type{<:AbstractStorage}, field::Symbol) = custom_field(S, Val(field))

custom_field(::Type{<:AbstractStorage}, ::Val{field}) where {field} = false

function get_vtk_filebase(body::AbstractBody, root::AbstractString)
    body_name = replace(string(get_name(body)), " " => "_")
    filebase = isempty(body_name) ? "timestep" : body_name * "_timestep"
    vtk_filebase::String = joinpath(root, "vtk", filebase)
    return vtk_filebase
end

function get_vtk_filebase(ms::AbstractMultibodySetup, root::AbstractString)
    vtk_filebase = Dict{Symbol,String}()
    for body_name in each_body_name(ms)
        vtk_filebase[body_name] = get_vtk_filebase(get_body(ms, body_name), root)
    end
    return vtk_filebase
end

function _export_results(options::AbstractJobOptions, chunk::AbstractBodyChunk,
                         chunk_id::Int, n_chunks::Int, n::Int, t::Float64)
    filename = get_filename(options, chunk.body_name, n)
    position = get_loc_position(chunk)
    pvtk_grid(filename, position, chunk.cells; part=chunk_id, nparts=n_chunks) do vtk
        export_fields!(vtk, chunk, options.fields, t)
        vtk["time", VTKFieldData()] = t
    end
    return nothing
end

@inline function get_filename(options::AbstractJobOptions, body_name::Symbol, n::Int)
    return _get_filename(options.vtk_filebase, body_name, n)
end

@inline function _get_filename(vtk_filebase::String, ::Symbol, n::Int)
    return @sprintf("%s_%d", vtk_filebase, n)
end

@inline function _get_filename(vtk_filebase::Dict{Symbol,String}, body_name::Symbol, n::Int)
    return @sprintf("%s_%d", vtk_filebase[body_name], n)
end

function export_fields!(vtk, chunk, fields::Vector{Symbol}, t)
    (; mat, system, storage, paramsetup) = chunk
    for field in fields
        point_data = export_field(Val(field), mat, system, storage, paramsetup, t)
        vtk[string(field), VTKPointData()] = point_data
    end
    return nothing
end

function export_fields!(vtk, chunk, fields_spec::Dict{Symbol,Vector{Symbol}}, t)
    fields = fields_spec[chunk.body_name]
    export_fields!(vtk, chunk, fields, t)
    return nothing
end

"""
    export_field(::Val{field}, mat, system, storage, paramsetup, t)

$(extension_api_note())

Return the point data that is written to the VTK file for `field`. The default returns the
local points of the storage field of that name, so a field of the storage is exported without
any further work.

Specialize it to derive a quantity that is not a storage field, or to reduce a bond field to
a point field. A derived name also has to be announced with [`custom_field`](@ref).

The returned array must have one column per *local* point, i.e. `get_n_loc_points(system)`
of them; halo entries are owned by another chunk and must not be exported twice.

# Example

```julia
Peridynamics.custom_field(::Type{<:MyStorage}, ::Val{:bond_damage_avg}) = true

function Peridynamics.export_field(::Val{:bond_damage_avg}, mat, system, storage::MyStorage,
                                   paramsetup, t)
    n = Peridynamics.get_n_loc_points(system)
    out = zeros(n)
    for i in 1:n
        bond_ids = Peridynamics.each_bond_idx(system, i)
        out[i] = sum(@view storage.bond_damage[bond_ids]) / length(bond_ids)
    end
    return out
end
```
"""
function export_field(::Val{field}, mat, system, storage, paramsetup, t) where {field}
    return get_loc_point_data(storage, system, field)
end

@inline function get_loc_position(chunk::AbstractBodyChunk)
    return @views chunk.storage.position[:, 1:get_n_loc_points(chunk)]
end

function msg_export_fields(fields::Vector{Symbol}; indentation=2)
    return msg_vec("exported fields", fields; continuation_label=" "^15, indentation)
end

function msg_export_fields(fields_spec::Dict{Symbol,Vector{Symbol}}; indentation=2)
    msg = " "^indentation * "EXPORTED FIELDS PER BODY\n"
    for (body_name, fields) in fields_spec
        body_name_str = string(body_name)
        n = length(body_name_str)
        ind = indentation + 2
        msg *= msg_vec(body_name_str, fields; continuation_label=" "^n, indentation=ind)
    end
    return msg
end
