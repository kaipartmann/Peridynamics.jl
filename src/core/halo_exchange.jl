
struct HaloExchange
    field::Symbol
    src_chunk_id::Int
    dest_chunk_id::Int
    src_idxs::Vector{Int}
    dest_idxs::Vector{Int}
end

function find_halo_exchanges(chunks::Vector{B}) where {B<:AbstractBodyChunk}
    read_fields, write_fields = get_halo_fields(chunks)
    read_halo_exs = Vector{Vector{HaloExchange}}(undef, length(chunks))
    write_halo_exs = Vector{Vector{HaloExchange}}(undef, length(chunks))
    @threads :static for chunk_id in eachindex(chunks)
        read_exs = Vector{HaloExchange}()
        write_exs = Vector{HaloExchange}()
        chunk = chunks[chunk_id]
        for (halo_chunk_id, idxs) in chunk.ch.halo_by_src
            halo_chunk = chunks[halo_chunk_id]
            halo_idxs = chunk.ch.point_ids[idxs]
            localize!(halo_idxs, halo_chunk.ch.localizer)
            find_read_exs!(read_exs, read_fields, chunk_id, idxs, halo_chunk_id, halo_idxs)
            find_write_exs!(write_exs, write_fields, chunk_id, idxs, halo_chunk_id,
                            halo_idxs)
        end
        read_halo_exs[chunk_id] = read_exs
        write_halo_exs[chunk_id] = write_exs
    end
    reorder_write_exs!(write_halo_exs)
    return read_halo_exs, write_halo_exs
end

# function find_read_exchanges(chunks, read_fields, dest_chunk_id)
#     halo_exs = Vector{HaloExchange}()
#     dest_chunk = chunks[dest_chunk_id]
#     for (src_chunk_id, dest_idxs) in dest_chunk.ch.halo_by_src
#         src_chunk = chunks[src_chunk_id]
#         src_idxs = dest_chunk.ch.point_ids[dest_idxs] # init with global idxs
#         localize!(src_idxs, src_chunk.ch.localizer)
#         for field in read_fields
#             he = HaloExchange(field, src_chunk_id, dest_chunk_id, src_idxs, dest_idxs)
#             push!(halo_exs, he)
#         end
#     end
#     return halo_exs
# end

function find_read_exs!(exs::Vector{HaloExchange}, read_fields::NTuple{N,Symbol},
                        chunk_id::Int, idxs::AbstractVector{Int}, halo_chunk_id::Int,
                        halo_idxs::AbstractVector{Int}) where {N}
    for field in read_fields
        push!(exs, HaloExchange(field, halo_chunk_id, chunk_id, halo_idxs, idxs))
    end
    return nothing
end

function find_write_exs!(exs::Vector{HaloExchange}, write_fields::NTuple{N,Symbol},
                         chunk_id::Int, idxs::AbstractVector{Int}, halo_chunk_id::Int,
                         halo_idxs::AbstractVector{Int}) where {N}
    for field in write_fields
        push!(exs, HaloExchange(field, chunk_id, halo_chunk_id, idxs, halo_idxs))
    end
    return nothing
end

# function find_write_exchanges(chunks, write_fields, src_chunk_id)
#     halo_exs = Vector{HaloExchange}()
#     src_chunk = chunk[src_chunk_id]
#     for (dest_chunk_id, src_idxs) in src_chunk.ch.halo_by_src
#         dest_chunk = chunks[dest_chunk_id]
#         dest_idxs = src_chunk.ch.point_ids[src_idxs]
#         localize!(dest_idxs, dest_chunk.ch.localizer)
#     end
#     return halo_exs
# end

function get_halo_fields(chunks::Vector{B}) where {B<:AbstractBodyChunk}
    material_type = get_material_type(chunks)
    read_fields = reads_from_halo(material_type)
    write_fields = writes_to_halo(material_type)
    return read_fields, write_fields
end

function reorder_write_exs!(write_halo_exs::Vector{Vector{HaloExchange}})
    all_write_exs = reduce(vcat, write_halo_exs)
    @threads :static for chunk_id in eachindex(write_halo_exs)
        write_halo_exs[chunk_id] = filter(x -> x.dest_chunk_id == chunk_id, all_write_exs)
    end
    return nothing
end

# function filter_exchanges_by_dest(halo_exs::Vector{HaloExchange}, n_chunks::Int)
#     return nothing
# end
