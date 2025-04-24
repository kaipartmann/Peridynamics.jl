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
    databcs::Vector{DataBC}
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
    databcs = body.data_bcs
    cells = get_cells(get_n_loc_points(system))
    chunk = BodyChunk(body_name, system, mat, paramsetup, storage, psets, sdbcs,
                      pdsdbcs, databcs, cells)
    return chunk
end

@inline function get_params(b::BodyChunk, point_id::Int)
    return get_params(b.paramsetup, point_id)
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
