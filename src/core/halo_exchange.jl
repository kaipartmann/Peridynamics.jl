
struct HaloExchange
    src_chunk_id::Int
    dest_chunk_id::Int
    src_idxs::Vector{Int}
    dest_idxs::Vector{Int}
end

function find_halo_exchanges(chunks::Vector{B}) where {B<:AbstractBodyChunk}
    has_read = has_read_halos(chunks)
    has_write = has_write_halos(chunks)
    read_halo_exs = Vector{Vector{HaloExchange}}(undef, length(chunks))
    write_halo_exs = Vector{Vector{HaloExchange}}(undef, length(chunks))
    @threads :static for chunk_id in eachindex(chunks)
        read_exs, write_exs = find_exs(chunks, chunk_id, has_read, has_write)
        read_halo_exs[chunk_id] = read_exs
        write_halo_exs[chunk_id] = write_exs
    end
    reorder_write_exs!(write_halo_exs)
    return read_halo_exs, write_halo_exs
end

function has_read_halos(chunk::AbstractBodyChunk)
    hrfs = get_halo_read_fields(chunk.store)
    isempty(hrfs) && return false
    return true
end

function has_read_halos(chunks::Vector{B}) where {B<:AbstractBodyChunk}
    return has_read_halos(first(chunks))
end

function has_write_halos(chunk::AbstractBodyChunk)
    hwfs = get_halo_write_fields(chunk.store)
    isempty(hwfs) && return false
    return true
end

function has_write_halos(chunks::Vector{B}) where {B<:AbstractBodyChunk}
    return has_write_halos(first(chunks))
end

function find_exs(chunks::Vector{B}, chunk_id::Int, hasread::Bool,
                  haswrite::Bool) where {B<:AbstractBodyChunk}
    readexs = Vector{HaloExchange}()
    writeexs = Vector{HaloExchange}()
    chunk = chunks[chunk_id]
    for (halo_chunk_id, idxs) in chunk.ch.halo_by_src
        halo_chunk = chunks[halo_chunk_id]
        halo_idxs = chunk.ch.point_ids[idxs]
        localize!(halo_idxs, halo_chunk.ch.localizer)
        hasread && push!(readexs, HaloExchange(halo_chunk_id, chunk_id, halo_idxs, idxs))
        haswrite && push!(writeexs, HaloExchange(chunk_id, halo_chunk_id, idxs, halo_idxs))
    end
    return readexs, writeexs
end

function reorder_write_exs!(write_halo_exs::Vector{Vector{HaloExchange}})
    all_write_exs = reduce(vcat, write_halo_exs)
    @threads :static for chunk_id in eachindex(write_halo_exs)
        write_halo_exs[chunk_id] = filter(x -> x.dest_chunk_id == chunk_id, all_write_exs)
    end
    return nothing
end
