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

function init_body_chunk(body::Body{M,P}, ts::T, loc_points::UnitRange{Int}) where {M,P,T}
    mat = body.mat
    discret, ch = init_discretization(body, loc_points)
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

function chop_body_threads(body::Body{M,P}, ts::T, point_decomp::V) where {M,P,T,V}
    D = discretization_type(body.mat)
    S = storage_type(body.mat, ts)
    body_chunks = Vector{BodyChunk{M,P,D,S}}(undef, length(point_decomp))

    @threads :static for cid in eachindex(point_decomp)
        loc_points = point_decomp[cid]
        body_chunks[cid] = init_body_chunk(body, ts, loc_points)
        #TODO: define precracks
        #TODO: apply initial conditions
    end

    return body_chunks
end

function chop_body_threads(body::MultibodySetup{M,P}, ts::T,
                           point_decomp::V) where {M,P,T,V}
    #TODO
    return nothing
end
