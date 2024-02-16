struct BodyChunk{M<:AbstractMaterial,P<:AbstractPointParameters,D<:AbstractDiscretization,
                 S<:AbstractStorage}
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
    loc_points = pd.decomp[chunk_id]
    mat = body.mat
    discret, ch = init_discretization(body, pd, chunk_id)
    storage = init_storage(body, ts, discret, ch)
    params_map = localize(body.params_map, ch)
    point_params = body.point_params[sort(unique(params_map))]
    point_sets = localized_point_sets(body.point_sets, ch)
    single_dim_bcs = body.single_dim_bcs
    return BodyChunk(mat, discret, storage, point_params, params_map, point_sets,
                     single_dim_bcs, ch)
end

@inline function get_point_param(b::BodyChunk, key::Symbol, i::Int)
    return getfield(b.point_params[b.params_map[i]], key)
end

function chop_body_threads(body::Body{M,P}, ts::T, pd::PointDecomposition) where {M,P,T}
    D = discretization_type(body.mat)
    S = storage_type(body.mat, ts)
    body_chunks = Vector{BodyChunk{M,P,D,S}}(undef, pd.n_chunks)

    @threads :static for chunk_id in eachindex(pd.decomp)
        body_chunk = init_body_chunk(body, ts, pd, chunk_id)
        #TODO: define precracks
        #TODO: apply initial conditions
        #TODO: get ExchangePulls
        body_chunks[chunk_id] = body_chunk
    end

    return body_chunks
end

function chop_body_threads(body::MultibodySetup{M,P}, ts::T,
                           point_decomp::V) where {M,P,T,V}
    #TODO
    return nothing
end
