struct BodyChunk{M<:AbstractMaterial,P<:AbstractPointParameters,D<:AbstractSystem,
                 S<:AbstractStorage} <: AbstractBodyChunk{M}
    mat::M
    system::D
    storage::S
    param::P
    psets::Dict{Symbol,Vector{Int}}
    sdbcs::Vector{SingleDimBC}
    pdsdbcs::Vector{PosDepSingleDimBC}
    ch::ChunkHandler
    cells::Vector{MeshCell{VTKCellType, Tuple{Int64}}}
end

function BodyChunk(body::Body{M,P}, ts::T, pd::PointDecomposition,
                   chunk_id::Int) where {M,P,T}
    mat = body.mat
    system, ch = init_system(body, pd, chunk_id)
    storage = init_storage(body, ts, system, ch)
    @assert length(body.point_params) == 1
    param = first(body.point_params)
    psets = localized_point_sets(body.point_sets, ch)
    sdbcs = body.single_dim_bcs
    pdsdbcs = body.posdep_single_dim_bcs
    cells = get_cells(ch.n_loc_points)
    return BodyChunk(mat, system, storage, param, psets, sdbcs, pdsdbcs, ch, cells)
end

@inline function get_param(b::BodyChunk, ::Int)
    return b.param
end

@inline get_material_type(::BodyChunk{M,P,D,S}) where {M,P,D,S} = M
@inline get_material_type(::Vector{BodyChunk{M,P,D,S}}) where {M,P,D,S} = M

@inline each_point_idx(b::AbstractBodyChunk) = each_point_idx(b.ch)
