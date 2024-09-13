struct ThreadsBodyDataHandler{Sys,M,P,S} <: AbstractThreadsBodyDataHandler{Sys,M,P,S}
    n_chunks::Int
    chunks::Vector{BodyChunk{Sys,M,P,S}}
    lth_exs::Vector{Vector{HaloExchange}}
    htl_exs::Vector{Vector{HaloExchange}}
end

function threads_data_handler(body::AbstractBody, solver::AbstractTimeSolver, n_chunks::Int)
    point_decomp = PointDecomposition(body, n_chunks)
    param_spec = get_param_spec(body)
    chunks = chop_body_threads(body, solver, point_decomp, param_spec)
    n_chunks = length(chunks)
    lth_exs, htl_exs = find_halo_exchanges(chunks)
    return ThreadsBodyDataHandler(n_chunks, chunks, lth_exs, htl_exs)
end

function chop_body_threads(body::AbstractBody, solver::AbstractTimeSolver,
                           point_decomp::PointDecomposition, param_spec::AbstractParamSpec)
    ChunkType = body_chunk_type(body, solver, param_spec)
    chunks = _chop_body_threads(ChunkType, body, solver, point_decomp, param_spec)
    return chunks
end

function _chop_body_threads(::Type{ChunkType}, body::AbstractBody,
                            solver::AbstractTimeSolver, point_decomp::PointDecomposition,
                            param_spec::AbstractParamSpec) where {ChunkType}
    chunks = Vector{ChunkType}(undef, point_decomp.n_chunks)
    @threads :static for chunk_id in eachindex(point_decomp.decomp)
        chunk = BodyChunk(body, solver, point_decomp, chunk_id, param_spec)
        apply_precracks!(chunk, body)
        apply_initial_conditions!(chunk, body)
        initialize!(chunk)
        chunks[chunk_id] = chunk
    end
    return chunks
end

function find_halo_exchanges(chunks::Vector{B}) where {B<:AbstractBodyChunk}
    lth_exs = Vector{Vector{HaloExchange}}(undef, length(chunks))
    htl_exs = Vector{Vector{HaloExchange}}(undef, length(chunks))
    @threads :static for chunk_id in eachindex(chunks)
        _lth_exs, _htl_exs = find_exs(chunks, chunk_id)
        lth_exs[chunk_id] = _lth_exs
        htl_exs[chunk_id] = _htl_exs
    end
    reorder_htl_exs!(htl_exs)
    return lth_exs, htl_exs
end

function find_exs(chunks::Vector{B}, chunk_id::Int) where {B<:AbstractBodyChunk}
    lth_exs = Vector{HaloExchange}()
    htl_exs = Vector{HaloExchange}()
    chunk = chunks[chunk_id]
    for (halo_chunk_id, idxs) in chunk.ch.hidxs_by_src
        halo_chunk = chunks[halo_chunk_id]
        halo_idxs = chunk.ch.point_ids[idxs]
        localize!(halo_idxs, halo_chunk.ch.localizer)
        push!(lth_exs, HaloExchange(0, halo_chunk_id, chunk_id, halo_idxs, idxs))
        push!(htl_exs, HaloExchange(0, chunk_id, halo_chunk_id, idxs, halo_idxs))
    end
    return lth_exs, htl_exs
end

function reorder_htl_exs!(htl_exs::Vector{Vector{HaloExchange}})
    all_htl_exs = reduce(vcat, htl_exs)
    @threads :static for chunk_id in eachindex(htl_exs)
        htl_exs[chunk_id] = filter(x -> x.dest_chunk_id == chunk_id, all_htl_exs)
    end
    return nothing
end

get_cells(n::Int) = [MeshCell(VTKCellTypes.VTK_VERTEX, (i,)) for i in 1:n]

function exchange_loc_to_halo!(dh::ThreadsBodyDataHandler, chunk_id::Int)
    fields = loc_to_halo_fields(dh.chunks[chunk_id].storage)
    isempty(fields) && return nothing
    exchange_loc_to_halo!(dh, chunk_id, fields)
    return nothing
end

function exchange_loc_to_halo!(dh::ThreadsBodyDataHandler, chunk_id::Int,
                               fields::NTuple{N,Symbol}) where {N}
    for field in fields
        exchange_loc_to_halo!(dh, chunk_id, field)
    end
    return nothing
end

function exchange_loc_to_halo!(dh::ThreadsBodyDataHandler, chunk_id::Int, field::Symbol)
    for ex in dh.lth_exs[chunk_id]
        dest_chunk = dh.chunks[ex.dest_chunk_id]
        src_chunk = dh.chunks[ex.src_chunk_id]
        _exchange_loc_to_halo!(dest_chunk, src_chunk, ex, field)
    end
    return nothing
end

function _exchange_loc_to_halo!(dest_chunk::C, src_chunk::C, ex::HaloExchange,
                                field::Symbol) where {C}
    dest_field = get_point_data(dest_chunk.storage, field)
    src_field = get_point_data(src_chunk.storage, field)
    exchange!(dest_field, src_field, ex.dest_idxs, ex.src_idxs)
    return nothing
end

function exchange_loc_to_halo!(get_field_function::F, dh::ThreadsBodyDataHandler,
                               chunk_id::Int) where {F<:Function}
    for ex in dh.lth_exs[chunk_id]
        dest_field = get_field_function(dh.chunks[ex.dest_chunk_id])
        src_field = get_field_function(dh.chunks[ex.src_chunk_id])
        exchange!(dest_field, src_field, ex.dest_idxs, ex.src_idxs)
    end
    return nothing
end

function exchange_halo_to_loc!(dh::ThreadsBodyDataHandler, chunk_id::Int)
    fields = halo_to_loc_fields(dh.chunks[chunk_id].storage)
    isempty(fields) && return nothing
    exchange_halo_to_loc!(dh, chunk_id, fields)
    return nothing
end

function exchange_halo_to_loc!(dh::ThreadsBodyDataHandler, chunk_id::Int,
                               fields::NTuple{N,Symbol}) where {N}
    for field in fields
        exchange_halo_to_loc!(dh, chunk_id, field)
    end
    return nothing
end

function exchange_halo_to_loc!(dh::ThreadsBodyDataHandler, chunk_id::Int, field::Symbol)
    for ex in dh.htl_exs[chunk_id]
        dest_chunk = dh.chunks[ex.dest_chunk_id]
        src_chunk = dh.chunks[ex.src_chunk_id]
        _exchange_halo_to_loc!(dest_chunk, src_chunk, ex, field)
    end
    return nothing
end

function _exchange_halo_to_loc!(dest_chunk::C, src_chunk::C, ex::HaloExchange,
                                field::Symbol) where {C}
    dest_field = get_point_data(dest_chunk.storage, field)
    src_field = get_point_data(src_chunk.storage, field)
    exchange_add!(dest_field, src_field, ex.dest_idxs, ex.src_idxs)
    return nothing
end

function exchange_halo_to_loc!(get_field_function::F, dh::ThreadsBodyDataHandler,
                               chunk_id::Int) where {F<:Function}
    for ex in dh.htl_exs[chunk_id]
        dest_field = get_field_function(dh.chunks[ex.dest_chunk_id])
        src_field = get_field_function(dh.chunks[ex.src_chunk_id])
        exchange_add!(dest_field, src_field, ex.dest_idxs, ex.src_idxs)
    end
    return nothing
end

function export_results(dh::ThreadsBodyDataHandler, options::AbstractJobOptions,
                        chunk_id::Int, timestep::Int, time::Float64)
    options.export_allowed || return nothing
    if mod(timestep, options.freq) == 0
        _export_results(options, dh.chunks[chunk_id], chunk_id, dh.n_chunks, timestep, time)
    end
    return nothing
end

function export_reference_results(dh::ThreadsBodyDataHandler, options::AbstractJobOptions)
    options.export_allowed || return nothing
    @threads :static for chunk_id in eachindex(dh.chunks)
        _export_results(options, dh.chunks[chunk_id], chunk_id, dh.n_chunks, 0, 0.0)
    end
    return nothing
end

function log_data_handler(options::AbstractJobOptions, dh::AbstractThreadsBodyDataHandler)
    log_system(options, dh)
    return nothing
end

@inline function system_type(dh::ThreadsBodyDataHandler{Sys}) where {Sys}
    return Sys
end

@inline function get_body_name(dh::AbstractThreadsBodyDataHandler)
    return dh.chunks[1].body_name
end
