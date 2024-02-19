
struct HaloExchange
    field::Symbol
    from_chunk_id::Int
    to_chunk_id::Int
    from_loc_idxs::Vector{Int}
    to_loc_idxs::Vector{Int}
end

function find_halo_exchanges(body_chunks::Vector{B}) where {B<:BodyChunk}
    halo_exchanges = Vector{HaloExchange}()
    find_read_exchanges!(halo_exchanges, body_chunks)
    find_write_exchanges!(halo_exchanges, body_chunks)
    return halo_exchanges
end

function find_read_exchanges!(halo_exchanges, body_chunks)
    halo_read_fields = reads_from_halo(first(body_chunks).mat)
    for read_field in halo_read_fields
        _find_read_exchanges!(halo_exchanges, body_chunks, read_field)
    end
    return nothing
end

function _find_read_exchanges!(halo_exchanges, body_chunks, read_field)
    for to_chunk_id in eachindex(body_chunks)
        to_chunk = body_chunks[to_chunk_id]
        for (from_chunk_id, to_loc_idxs) in to_chunk.ch.halo_by_src
            from_chunk = body_chunks[from_chunk_id]
            from_loc_idxs = to_chunk.ch.point_ids[to_loc_idxs] # init with global idxs
            localize!(from_loc_idxs, from_chunk.ch.localizer)
            he = HaloExchange(read_field, from_chunk_id, to_chunk_id, from_loc_idxs,
                              to_loc_idxs)
            push!(halo_exchanges, he)
        end
    end
    return nothing
end

function find_write_exchanges!(halo_exchanges, body_chunks)
    # TODO
    # halo_write_fields = writes_to_halo(body_chunk.mat)
    # for write_field in halo_write_fields
    # end
    return nothing
end
