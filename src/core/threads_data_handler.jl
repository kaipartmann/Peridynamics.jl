struct ThreadsDataHandler{Sys,M,P,S} <: AbstractThreadsDataHandler{Sys,M,P,S}
    n_chunks::Int
    chunks::Vector{BodyChunk{Sys,M,P,S}}
    lth_exs::Vector{Vector{HaloExchange}}
    htl_exs::Vector{Vector{HaloExchange}}
end

function ThreadsDataHandler(body::AbstractBody, solver::AbstractTimeSolver,
                            point_decomp::PointDecomposition)
    param_spec = get_param_spec(body)
    chunks = chop_body_threads(body, solver, point_decomp, param_spec)
    n_chunks = length(chunks)
    lth_exs, htl_exs = find_halo_exchanges(chunks)
    return ThreadsDataHandler(n_chunks, chunks, lth_exs, htl_exs)
end

function ThreadsDataHandler(multibody::AbstractMultibodySetup, solver::AbstractTimeSolver,
                            point_decomp::PointDecomposition)
    error("MultibodySetup not yet implemented!\n")
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

function calc_stable_timestep(dh::ThreadsDataHandler, safety_factor::Float64)
    Δt = zeros(length(dh.chunks))
    @threads :static for chunk_id in eachindex(dh.chunks)
        Δt[chunk_id] = calc_timestep(dh.chunks[chunk_id])
    end
    return minimum(Δt) * safety_factor
end

get_cells(n::Int) = [MeshCell(VTKCellTypes.VTK_VERTEX, (i,)) for i in 1:n]

function exchange_loc_to_halo!(dh::ThreadsDataHandler, chunk_id::Int)
    fields = loc_to_halo_fields(dh.chunks[chunk_id].storage)
    isempty(fields) && return nothing
    exchange_loc_to_halo!(dh, chunk_id, fields)
    return nothing
end

function exchange_loc_to_halo!(dh::ThreadsDataHandler, chunk_id::Int,
                               fields::NTuple{N,Symbol}) where {N}
    for field in fields
        exchange_loc_to_halo!(dh, chunk_id, field)
    end
    return nothing
end

function exchange_loc_to_halo!(dh::ThreadsDataHandler, chunk_id::Int, field::Symbol)
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

function exchange_loc_to_halo!(get_field_function::F, dh::ThreadsDataHandler,
                               chunk_id::Int) where {F<:Function}
    for ex in dh.lth_exs[chunk_id]
        dest_field = get_field_function(dh.chunks[ex.dest_chunk_id])
        src_field = get_field_function(dh.chunks[ex.src_chunk_id])
        exchange!(dest_field, src_field, ex.dest_idxs, ex.src_idxs)
    end
    return nothing
end

function exchange_halo_to_loc!(dh::ThreadsDataHandler, chunk_id::Int)
    fields = halo_to_loc_fields(dh.chunks[chunk_id].storage)
    isempty(fields) && return nothing
    exchange_halo_to_loc!(dh, chunk_id, fields)
    return nothing
end

function exchange_halo_to_loc!(dh::ThreadsDataHandler, chunk_id::Int,
                               fields::NTuple{N,Symbol}) where {N}
    for field in fields
        exchange_halo_to_loc!(dh, chunk_id, field)
    end
    return nothing
end

function exchange_halo_to_loc!(dh::ThreadsDataHandler, chunk_id::Int, field::Symbol)
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

function exchange_halo_to_loc!(get_field_function::F, dh::ThreadsDataHandler,
                               chunk_id::Int) where {F<:Function}
    for ex in dh.htl_exs[chunk_id]
        dest_field = get_field_function(dh.chunks[ex.dest_chunk_id])
        src_field = get_field_function(dh.chunks[ex.src_chunk_id])
        exchange_add!(dest_field, src_field, ex.dest_idxs, ex.src_idxs)
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

function initialize!(::AbstractThreadsDataHandler, ::AbstractTimeSolver)
    return nothing
end

function log_data_handler(options::AbstractOptions,
                          dh::AbstractThreadsDataHandler{Sys,M,P,S}) where {Sys,M,P,S}
    msg = "THREADS DATA HANDLER\n"
    msg *= log_msg("system type", Sys)
    msg *= log_msg("material type", M)
    # msg *= log_msg("parameter type", P)
    # msg *= log_msg("storage type", S)
    log_it(options, msg)
    return nothing
end
