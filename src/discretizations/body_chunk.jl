struct BodyChunk{M<:AbstractMaterial,P<:AbstractPointParameters,D<:AbstractDiscretization,
                 S<:AbstractStorage} <: AbstractBodyChunk
    mat::M
    discret::D
    storage::S
    point_params::Vector{P}
    params_map::Vector{Int} # TODO: BodyChunk type with only 1 point param for entire body
    point_sets::Dict{Symbol,Vector{Int}}
    single_dim_bcs::Vector{SingleDimBC}
    ch::ChunkHandler
end

function init_body_chunk(body::Body{M,P}, ts::T, pd::PointDecomposition,
                         chunk_id::Int) where {M,P,T}
    mat = body.mat
    discret, ch = init_discretization(body, pd, chunk_id)
    storage = init_storage(body, ts, discret, ch)
    params_map = body.params_map[ch.loc_points]
    point_params = body.point_params
    point_sets = localized_point_sets(body.point_sets, ch)
    single_dim_bcs = body.single_dim_bcs
    return BodyChunk(mat, discret, storage, point_params, params_map, point_sets,
                     single_dim_bcs, ch)
end

@inline function get_point_param(b::BodyChunk, key::Symbol, i::Int)
    return getfield(b.point_params[b.params_map[i]], key)
end

@inline function get_point_param(b::BodyChunk, point_id::Int)
    return b.point_params[b.params_map[point_id]]
end

function chop_body_threads(body::Body{M,P}, ts::T, pd::PointDecomposition) where {M,P,T}
    D = discretization_type(body.mat)
    S = storage_type(body.mat, ts)
    body_chunks = Vector{BodyChunk{M,P,D,S}}(undef, pd.n_chunks)

    @threads :static for chunk_id in eachindex(pd.decomp)
        body_chunk = init_body_chunk(body, ts, pd, chunk_id)
        apply_precracks!(body_chunk, body)
        apply_initial_conditions!(body_chunk, body)
        body_chunks[chunk_id] = body_chunk
    end

    return body_chunks
end

function chop_body_threads(body::MultibodySetup{M,P}, ts::T,
                           point_decomp::V) where {M,P,T,V}
    #TODO
    return nothing
end

@inline each_point_idx(b::BodyChunk) = each_point_idx(b.ch)

function apply_precracks!(b::BodyChunk, body::Body)
    for precrack in body.point_sets_precracks
        apply_precrack!(b, precrack)
    end
    return nothing
end

@inline function calc_damage!(b::BodyChunk)
    for point_id in each_point_idx(b)
        dmg = 1 - b.n_active_bonds[point_id] / b.discret.n_neighbors[point_id]
        b.damage[point_id] = dmg
    end
    return nothing
end

function apply_initial_conditions!(b::BodyChunk, body::Body)
    apply_single_dim_ic!(b, body)
    return nothing
end

@inline function apply_single_dim_ic!(b::BodyChunk, body::Body)
    for ic in body.single_dim_ics
        apply_ic!(b, ic)
    end
    return nothing
end

function apply_precrack!(b::AbstractBodyChunk, precrack::PointSetsPreCrack)
    set_a = b.point_sets[precrack.set_a]
    set_b = b.point_sets[precrack.set_b]
    if isempty(set_a) || isempty(set_b)
        return nothing
    end
    _apply_precrack!(b.storage, b.discret, b.ch, set_a, set_b)
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
