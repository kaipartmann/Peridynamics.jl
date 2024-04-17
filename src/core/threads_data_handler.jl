struct HaloExchange
    src_chunk_id::Int
    dest_chunk_id::Int
    src_idxs::Vector{Int}
    dest_idxs::Vector{Int}
end

struct ThreadsDataHandler{C<:AbstractBodyChunk} <: AbstractThreadsDataHandler
    n_chunks::Int
    chunks::Vector{C}
    loc_to_halo_exs::Vector{Vector{HaloExchange}}
    halo_to_loc_exs::Vector{Vector{HaloExchange}}
end

function ThreadsDataHandler(body::AbstractBody, time_solver::AbstractTimeSolver,
                            point_decomp::PointDecomposition)
    n_param = length(body.point_params)
    v = n_param == 1 ? Val{1}() : Val{2}()
    return ThreadsDataHandler(body, time_solver, point_decomp, v)
end

function ThreadsDataHandler(body::AbstractBody, time_solver::AbstractTimeSolver,
                            point_decomp::PointDecomposition, v::Val{N}) where {N}
    chunks = chop_body_threads(body, time_solver, point_decomp, v)
    n_chunks = length(chunks)
    loc_to_halo_exs, halo_to_loc_exs = find_halo_exchanges(chunks)
    return ThreadsDataHandler(n_chunks, chunks, loc_to_halo_exs, halo_to_loc_exs)
end

function ThreadsDataHandler(multibody::AbstractMultibodySetup,
                            time_solver::AbstractTimeSolver,
                            point_decomp::PointDecomposition)
    error("MultibodySetup not yet implemented!\n")
end

function chop_body_threads(body::Body{M,P}, ts::T, pd::PointDecomposition,
                           v::Val{N}) where {M,P,T,N}
    D = system_type(body.mat)
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

function find_halo_exchanges(chunks::Vector{B}) where {B<:AbstractBodyChunk}
    has_lth = has_loc_to_halo_exs(chunks)
    has_htl = has_halo_to_loc_exs(chunks)
    loc_to_halo_exs = Vector{Vector{HaloExchange}}(undef, length(chunks))
    halo_to_loc_exs = Vector{Vector{HaloExchange}}(undef, length(chunks))
    @threads :static for chunk_id in eachindex(chunks)
        lth_exs, htl_exs = find_exs(chunks, chunk_id, has_lth, has_htl)
        loc_to_halo_exs[chunk_id] = lth_exs
        halo_to_loc_exs[chunk_id] = htl_exs
    end
    reorder_htl_exs!(halo_to_loc_exs)
    return loc_to_halo_exs, halo_to_loc_exs
end

function has_loc_to_halo_exs(chunks::Vector{B}) where {B<:AbstractBodyChunk}
    return has_loc_to_halo_exs(first(chunks))
end

function has_halo_to_loc_exs(chunks::Vector{B}) where {B<:AbstractBodyChunk}
    return has_halo_to_loc_exs(first(chunks))
end

function find_exs(chunks::Vector{B}, chunk_id::Int, has_lth::Bool,
                  has_htl::Bool) where {B<:AbstractBodyChunk}
    lth_exs = Vector{HaloExchange}()
    htl_exs = Vector{HaloExchange}()
    chunk = chunks[chunk_id]
    for (halo_chunk_id, idxs) in chunk.ch.hidxs_by_src
        halo_chunk = chunks[halo_chunk_id]
        halo_idxs = chunk.ch.point_ids[idxs]
        localize!(halo_idxs, halo_chunk.ch.localizer)
        has_lth && push!(lth_exs, HaloExchange(halo_chunk_id, chunk_id, halo_idxs, idxs))
        has_htl && push!(htl_exs, HaloExchange(chunk_id, halo_chunk_id, idxs, halo_idxs))
    end
    return lth_exs, htl_exs
end

function reorder_htl_exs!(halo_to_loc_exs::Vector{Vector{HaloExchange}})
    all_write_exs = reduce(vcat, halo_to_loc_exs)
    @threads :static for chunk_id in eachindex(halo_to_loc_exs)
        halo_to_loc_exs[chunk_id] = filter(x -> x.dest_chunk_id == chunk_id, all_write_exs)
    end
    return nothing
end

function calc_stable_timestep(dh::ThreadsDataHandler, safety_factor::Float64)
    Δt = zeros(length(dh.chunks))
    @threads :static for chunk_id in eachindex(dh.chunks)
        Δt[chunk_id] = calc_timestep(dh.chunks[chunk_id])
    end
    return minimum(Δt) * safety_factor
end

get_cells(n::Int) = [MeshCell(VTKCellTypes.VTK_VERTEX, (i,)) for i in 1:n]

function exchange_loc_to_halo!(dh::ThreadsDataHandler, chunk_id::Int)
    for he in dh.loc_to_halo_exs[chunk_id]
        dest_storage = dh.chunks[he.dest_chunk_id].storage
        src_storage = dh.chunks[he.src_chunk_id].storage
        _exchange_loc_to_halo!(dest_storage, src_storage, he)
    end
    return nothing
end

function _exchange_loc_to_halo!(dest_storage::S, src_storage::S, he::HaloExchange) where {S}
    for field in halo_read_fields(dest_storage)
        dest_field = get_point_data(dest_storage, field)
        src_field = get_point_data(src_storage, field)
        exchange!(dest_field, src_field, he.dest_idxs, he.src_idxs)
    end
    return nothing
end

function exchange_halo_to_loc!(dh::ThreadsDataHandler, chunk_id::Int)
    for he in dh.halo_to_loc_exs[chunk_id]
        dest_storage = dh.chunks[he.dest_chunk_id].storage
        src_storage = dh.chunks[he.src_chunk_id].storage
        _exchange_halo_to_loc!(dest_storage, src_storage, he)
    end
    return nothing
end

function _exchange_halo_to_loc!(dest_storage::S, src_storage::S,
                                 he::HaloExchange) where {S}
    for field in halo_write_fields(dest_storage)
        dest_field = get_point_data(dest_storage, field)
        src_field = get_point_data(src_storage, field)
        exchange_add!(dest_field, src_field, he.dest_idxs, he.src_idxs)
    end
    return nothing
end

function export_results(dh::ThreadsDataHandler, options::AbstractOptions, chunk_id::Int,
                        timestep::Int, time::Float64)
    options.exportflag || return nothing
    if mod(timestep, options.freq) == 0
        _export_results(dh.chunks[chunk_id], chunk_id, dh.n_chunks, options, timestep, time)
    end
    return nothing
end

function export_reference_results(dh::ThreadsDataHandler, options::AbstractOptions)
    options.exportflag || return nothing
    @threads :static for chunk_id in eachindex(dh.chunks)
        _export_results(dh.chunks[chunk_id], chunk_id, dh.n_chunks, options, 0, 0.0)
    end
    return nothing
end
