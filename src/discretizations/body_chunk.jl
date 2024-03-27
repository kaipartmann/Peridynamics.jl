struct MultiParamBodyChunk{M<:AbstractMaterial,P<:AbstractPointParameters,
                           D<:AbstractSystem,
                           S<:AbstractStorage} <: AbstractBodyChunk
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
    system, ch = init_discretization(body, pd, chunk_id)
    storage = init_storage(body, ts, system, ch)
    paramap = body.params_map[ch.loc_points]
    param = body.point_params
    psets = localized_point_sets(body.point_sets, ch)
    sdbcs = body.single_dim_bcs
    pdsdbcs = body.posdep_single_dim_bcs
    cells = get_cells(ch.n_loc_points)
    return MultiParamBodyChunk(mat, system, storage, param, paramap, psets, sdbcs, pdsdbcs, ch,
                               cells)
end

struct BodyChunk{M<:AbstractMaterial,P<:AbstractPointParameters,D<:AbstractSystem,
                 S<:AbstractStorage} <: AbstractBodyChunk
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
    system, ch = init_discretization(body, pd, chunk_id)
    storage = init_storage(body, ts, system, ch)
    @assert length(body.point_params) == 1
    param = first(body.point_params)
    psets = localized_point_sets(body.point_sets, ch)
    sdbcs = body.single_dim_bcs
    pdsdbcs = body.posdep_single_dim_bcs
    cells = get_cells(ch.n_loc_points)
    return BodyChunk(mat, system, storage, param, psets, sdbcs, pdsdbcs, ch, cells)
end

@inline function get_param(b::MultiParamBodyChunk, point_id::Int)
    return b.param[b.paramap[point_id]]
end

@inline function get_param(b::BodyChunk, ::Int)
    return b.param
end

@inline get_material_type(::MultiParamBodyChunk{M,P,D,S}) where {M,P,D,S} = M
@inline get_material_type(::BodyChunk{M,P,D,S}) where {M,P,D,S} = M
@inline get_material_type(::Vector{MultiParamBodyChunk{M,P,D,S}}) where {M,P,D,S} = M
@inline get_material_type(::Vector{BodyChunk{M,P,D,S}}) where {M,P,D,S} = M

@inline each_point_idx(b::AbstractBodyChunk) = each_point_idx(b.ch)

@inline function calc_damage!(b::AbstractBodyChunk)
    for point_id in each_point_idx(b)
        dmg = 1 - b.storage.n_active_bonds[point_id] / b.system.n_neighbors[point_id]
        b.storage.damage[point_id] = dmg
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
    _apply_precrack!(b.storage, b.system, b.ch, set_a, set_b)
    return nothing
end

function _apply_precrack!(s::AbstractStorage, bd::BondSystem, ch::ChunkHandler,
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
