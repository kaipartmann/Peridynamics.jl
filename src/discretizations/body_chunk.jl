struct MultiParamBodyChunk{M<:AbstractMaterial,P<:AbstractPointParameters,
                           D<:AbstractDiscretization,
                           S<:AbstractStorage} <: AbstractBodyChunk
    mat::M
    dscr::D
    store::S
    param::Vector{P}
    paramap::Vector{Int}
    psets::Dict{Symbol,Vector{Int}}
    sdbcs::Vector{SingleDimBC}
    ch::ChunkHandler
    cells::Vector{MeshCell{VTKCellType, Tuple{Int64}}}
end

function MultiParamBodyChunk(body::Body{M,P}, ts::T, pd::PointDecomposition,
                             chunk_id::Int) where {M,P,T}
    mat = body.mat
    dscr, ch = init_discretization(body, pd, chunk_id)
    store = init_storage(body, ts, dscr, ch)
    paramap = body.params_map[ch.loc_points]
    param = body.point_params
    psets = localized_point_sets(body.point_sets, ch)
    sdbcs = body.single_dim_bcs
    cells = get_cells(ch.n_loc_points)
    return MultiParamBodyChunk(mat, dscr, store, param, paramap, psets, sdbcs, ch, cells)
end

struct BodyChunk{M<:AbstractMaterial,P<:AbstractPointParameters,D<:AbstractDiscretization,
                 S<:AbstractStorage} <: AbstractBodyChunk
    mat::M
    dscr::D
    store::S
    param::P
    psets::Dict{Symbol,Vector{Int}}
    sdbcs::Vector{SingleDimBC}
    ch::ChunkHandler
    cells::Vector{MeshCell{VTKCellType, Tuple{Int64}}}
end

function BodyChunk(body::Body{M,P}, ts::T, pd::PointDecomposition,
                   chunk_id::Int) where {M,P,T}
    mat = body.mat
    dscr, ch = init_discretization(body, pd, chunk_id)
    store = init_storage(body, ts, dscr, ch)
    @assert length(body.point_params) == 1
    param = first(body.point_params)
    psets = localized_point_sets(body.point_sets, ch)
    sdbcs = body.single_dim_bcs
    cells = get_cells(ch.n_loc_points)
    return BodyChunk(mat, dscr, store, param, psets, sdbcs, ch, cells)
end

@inline function get_param(b::MultiParamBodyChunk, point_id::Int)
    return b.param[b.paramap[point_id]]
end

@inline function get_param(b::BodyChunk, ::Int)
    return b.param
end

function chop_body_threads(body::Body{M,P}, ts::T, pd::PointDecomposition,
                           v::Val{N}) where {M,P,T,N}
    D = discretization_type(body.mat)
    S = storage_type(body.mat, ts)
    body_chunks = _chop_body_threads(body, ts, pd, D, S, v)
    return body_chunks
end

function _chop_body_threads(body::Body{M,P}, ts::T, pd::PointDecomposition, ::Type{D},
                            ::Type{S}, ::Val{1}) where {M,P,D,S,T}
    body_chunks = Vector{BodyChunk{M,P,D,S}}(undef, pd.n_chunks)

    @threads :static for chunk_id in eachindex(pd.decomp)
        body_chunk = BodyChunk(body, ts, pd, chunk_id)
        apply_precracks!(body_chunk, body)
        apply_initial_conditions!(body_chunk, body)
        body_chunks[chunk_id] = body_chunk
    end

    return body_chunks
end

function _chop_body_threads(body::Body{M,P}, ts::T, pd::PointDecomposition, ::Type{D},
                            ::Type{S}, ::Val{N}) where {M,P,D,S,T,N}
    body_chunks = Vector{MultiParamBodyChunk{M,P,D,S}}(undef, pd.n_chunks)

    @threads :static for chunk_id in eachindex(pd.decomp)
        body_chunk = MultiParamBodyChunk(body, ts, pd, chunk_id)
        apply_precracks!(body_chunk, body)
        apply_initial_conditions!(body_chunk, body)
        body_chunks[chunk_id] = body_chunk
    end

    return body_chunks
end

@inline get_material_type(::MultiParamBodyChunk{M,P,D,S}) where {M,P,D,S} = M
@inline get_material_type(::BodyChunk{M,P,D,S}) where {M,P,D,S} = M
@inline get_material_type(::Vector{MultiParamBodyChunk{M,P,D,S}}) where {M,P,D,S} = M
@inline get_material_type(::Vector{BodyChunk{M,P,D,S}}) where {M,P,D,S} = M

@inline each_point_idx(b::AbstractBodyChunk) = each_point_idx(b.ch)

@inline function calc_damage!(b::AbstractBodyChunk)
    for point_id in each_point_idx(b)
        dmg = 1 - b.store.n_active_bonds[point_id] / b.dscr.n_neighbors[point_id]
        b.store.damage[point_id] = dmg
    end
    return nothing
end

function apply_initial_conditions!(b::AbstractBodyChunk, body::Body)
    apply_single_dim_ic!(b, body)
    return nothing
end

@inline function apply_single_dim_ic!(b::AbstractBodyChunk, body::Body)
    for ic in body.single_dim_ics
        apply_ic!(b, ic)
    end
    return nothing
end

function apply_precracks!(b::AbstractBodyChunk, body::Body)
    for precrack in body.point_sets_precracks
        apply_precrack!(b, body, precrack)
    end
    calc_damage!(b)
    return nothing
end

function apply_precrack!(b::AbstractBodyChunk, body::Body, precrack::PointSetsPreCrack)
    set_a = filter(x -> in(x, b.ch.point_ids), body.point_sets[precrack.set_a])
    set_b = filter(x -> in(x, b.ch.point_ids), body.point_sets[precrack.set_b])
    localize!(set_a, b.ch.localizer)
    localize!(set_b, b.ch.localizer)
    if isempty(set_a) || isempty(set_b)
        return nothing
    end
    _apply_precrack!(b.store, b.dscr, b.ch, set_a, set_b)
    return nothing
end

function _apply_precrack!(s::AbstractStorage, bd::BondDiscretization, ch::ChunkHandler,
                          set_a::Vector{Int}, set_b::Vector{Int})
    s.n_active_bonds .= 0
    for point_id in each_point_idx(ch)
        for bond_id in each_bond_idx(bd, point_id)
            bond = bd.bonds[bond_id]
            neighbor_id = bond.neighbor
            point_in_a = in(point_id, set_a)
            point_in_b = in(point_id, set_b)
            neigh_in_a = in(neighbor_id, set_a)
            neigh_in_b = in(neighbor_id, set_b)
            if (point_in_a && neigh_in_b) || (point_in_b && neigh_in_a)
                s.bond_active[bond_id] = false
            end
            s.n_active_bonds[point_id] += s.bond_active[bond_id]
        end
    end
    return nothing
end
