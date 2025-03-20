"""
    BodyChunk{System, Material, Params, Storage}

$(internal_api_warning())

A type that contains all data of a body chunk. For parallel simulations, the body is divided into multiple chunks. Each `BodyChunk` instance contains all necessary information for the simulation on this specific chunk. This type is used for multithreading and MPI.

# Type Parameters

- `System<:AbstractSystem`: Type of the system.
- `Material<:AbstractMaterial`: Type of the material model of the system.
- `Params<:AbstractParameterSetup`: Material parameters of the points in the body chunk.
- `Storage<:AbstractStorage`: Storage of all information that changes during the simulation.

# Fields

- `body_name::Symbol`: Name of the body in multibody simulations.
- `system::System`: System with all information that is known before the simulation.
- `mat::Material`: Material model of the system.
- `paramsetup::Params`: Material parameters of the points in the body chunk.
- `storage::Storage`: Storage of all information that changes during the simulation.
- `psets::Dict{Symbol,Vector{Int}}`: Point sets of the chunk with local indices.
- `sdbcs::Vector{SingleDimBC}`: Single dimension boundary conditions.
- `pdsdbcs::Vector{PosDepSingleDimBC}`: Position dependent single dimension boundary
    conditions.
- `cells::Vector{MeshCell{VTKCellType,Tuple{Int64}}}`: Cells for vtk export.
"""
struct BodyChunk{System<:AbstractSystem,
                 Material<:AbstractMaterial,
                 Params<:AbstractParameterSetup,
                 Storage<:AbstractStorage} <: AbstractBodyChunk{System,Material}
    body_name::Symbol
    system::System
    mat::Material
    paramsetup::Params
    storage::Storage
    psets::Dict{Symbol,Vector{Int}}
    sdbcs::Vector{SingleDimBC}
    pdsdbcs::Vector{PosDepSingleDimBC}
    cells::Vector{MeshCell{VTKCellType,Tuple{Int64}}}
end

function BodyChunk(body::AbstractBody, solver::AbstractTimeSolver, pd::PointDecomposition,
                   chunk_id::Int, param_spec::AbstractParamSpec)
    body_name = get_name(body)
    mat = body.mat
    system = get_system(body, pd, chunk_id)
    paramsetup = get_paramsetup(body, system.chunk_handler, param_spec)
    storage = get_storage(mat, solver, system)
    psets = localized_point_sets(body.point_sets, system.chunk_handler)
    sdbcs = body.single_dim_bcs
    pdsdbcs = body.posdep_single_dim_bcs
    cells = get_cells(get_n_loc_points(system))
    chunk = BodyChunk(body_name, system, mat, paramsetup, storage, psets, sdbcs,
                      pdsdbcs, cells)
    return chunk
end

@inline function get_params(chunk::BodyChunk, point_id::Int)
    return get_params(chunk.paramsetup, point_id)
end

@inline function body_chunk_type(body::AbstractBody{Material}, solver::AbstractTimeSolver,
                                 param_spec::AbstractParamSpec) where {Material}
    System = system_type(body.mat)
    Params = parameter_setup_type(body, param_spec)
    Storage = storage_type(body)
    return BodyChunk{System,Material,Params,Storage}
end

@inline each_point_idx(chunk::AbstractBodyChunk) = each_point_idx(chunk.system)
@inline each_point_idx_pair(chunk::AbstractBodyChunk) = each_point_idx_pair(chunk.system)

@inline n_loc_points(chunk::AbstractBodyChunk) = get_n_loc_points(chunk.system)
@inline n_points(chunk::AbstractBodyChunk) = get_n_points(chunk.system)

function initialize!(::AbstractBodyChunk)
    return nothing
end
