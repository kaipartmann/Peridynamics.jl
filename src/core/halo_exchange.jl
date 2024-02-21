
struct HaloExchange
    field::Symbol
    src_chunk_id::Int
    dest_chunk_id::Int
    src_idxs::Vector{Int}
    dest_idxs::Vector{Int}
end

function find_halo_exchanges(body_chunks::Vector{B}) where {B<:BodyChunk}
    halo_exchanges_read = find_read_exchanges(body_chunks)
    halo_exchanges_write = find_write_exchanges(body_chunks)
    halo_exchanges = vcat(halo_exchanges_read, halo_exchanges_write)
    return halo_exchanges
end

function find_read_exchanges(body_chunks::Vector{B}) where {B<:BodyChunk}
    read_fields = reads_from_halo(first(body_chunks).mat)
    _halo_exchanges = Vector{Vector{HaloExchange}}(undef, length(body_chunks))
    @threads :static for chunk_id in eachindex(body_chunks)
        _halo_exchanges[chunk_id] = _find_read_exchanges(body_chunks, read_fields, chunk_id)
    end
    halo_exchanges = reduce(vcat, _halo_exchanges)
    return halo_exchanges
end

function _find_read_exchanges(body_chunks, read_fields, chunk_id)
    _halo_exchanges = Vector{HaloExchange}()
    to_chunk = body_chunks[chunk_id]
    for (from_chunk_id, to_loc_idxs) in to_chunk.ch.halo_by_src
        from_chunk = body_chunks[from_chunk_id]
        from_loc_idxs = to_chunk.ch.point_ids[to_loc_idxs] # init with global idxs
        localize!(from_loc_idxs, from_chunk.ch.localizer)
        for field in read_fields
            he = HaloExchange(field, from_chunk_id, chunk_id, from_loc_idxs, to_loc_idxs)
            push!(_halo_exchanges, he)
        end
    end
    return _halo_exchanges
end

function find_write_exchanges(body_chunks::Vector{B}) where {B<:BodyChunk}
    halo_exchanges = Vector{HaloExchange}()
    # TODO
    return halo_exchanges
end
