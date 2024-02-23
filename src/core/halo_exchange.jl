
struct HaloExchange
    field::Symbol
    src_chunk_id::Int
    dest_chunk_id::Int
    src_idxs::Vector{Int}
    dest_idxs::Vector{Int}
end

function find_halo_exchanges(chunks::Vector{B}) where {B<:BodyChunk}
    material_type = get_material_type(chunks)
    read_fields = reads_from_halo(material_type)
    write_fields = reads_from_halo(material_type)

    read_halo_exs = Vector{Vector{HaloExchange}}(undef, length(chunks))
    write_halo_exs = Vector{Vector{HaloExchange}}(undef, length(chunks))
    @threads :static for chunk_id in eachindex(chunks)
        read_halo_exs[chunk_id] = find_read_exchanges(chunks, read_fields, chunk_id)
        write_halo_exs[chunk_id] = find_write_exchanges(chunks, write_fields, chunk_id)
    end
    return read_halo_exs, write_halo_exs
end

function find_read_exchanges(chunks, read_fields, dest_chunk_id)
    halo_exs = Vector{HaloExchange}()
    dest_chunk = chunks[dest_chunk_id]
    for (src_chunk_id, dest_idxs) in dest_chunk.ch.halo_by_src
        src_chunk = chunks[src_chunk_id]
        src_idxs = dest_chunk.ch.point_ids[dest_idxs] # init with global idxs
        localize!(src_idxs, src_chunk.ch.localizer)
        for field in read_fields
            he = HaloExchange(field, src_chunk_id, dest_chunk_id, src_idxs, dest_idxs)
            push!(halo_exs, he)
        end
    end
    return halo_exs
end

function find_write_exchanges(chunks, write_fields, dest_chunk_id)
    halo_exchanges = Vector{HaloExchange}()
    # TODO
    return halo_exchanges
end

function filter_exchanges_by_dest(halo_exs::Vector{HaloExchange}, n_chunks::Int)
    return nothing
end
