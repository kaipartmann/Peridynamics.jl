struct MultiParamBodyChunk{M<:AbstractMaterial,P<:AbstractPointParameters,
                           D<:AbstractSystem,
                           S<:AbstractStorage} <: AbstractBodyChunk{M}
    mat::M
    system::D
    storage::S
    param::Vector{P}
    paramap::Vector{Int}
    psets::Dict{Symbol,Vector{Int}}
    sdbcs::Vector{SingleDimBC}
    pdsdbcs::Vector{PosDepSingleDimBC}
    ch::ChunkHandler
    cells::Vector{MeshCell{VTKCellType, Tuple{Int64}}}
end

function MultiParamBodyChunk(body::Body{M,P}, ts::T, pd::PointDecomposition,
                             chunk_id::Int) where {M,P,T}
    mat = body.mat
    system, ch = get_system(body, pd, chunk_id)
    storage = init_storage(mat, ts, system, ch)
    paramap = body.params_map[ch.loc_points]
    param = body.point_params
    psets = localized_point_sets(body.point_sets, ch)
    sdbcs = body.single_dim_bcs
    pdsdbcs = body.posdep_single_dim_bcs
    cells = get_cells(ch.n_loc_points)
    return MultiParamBodyChunk(mat, system, storage, param, paramap, psets, sdbcs, pdsdbcs, ch,
                               cells)
end

@inline function get_param(b::MultiParamBodyChunk, point_id::Int)
    return b.param[b.paramap[point_id]]
end

@inline get_material_type(::MultiParamBodyChunk{M,P,D,S}) where {M,P,D,S} = M
@inline get_material_type(::Vector{MultiParamBodyChunk{M,P,D,S}}) where {M,P,D,S} = M
