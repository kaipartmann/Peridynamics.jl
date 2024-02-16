
struct ThreadsHaloPull
    field::Symbol
    source::Int
    from_ids::Vector{Int}
    to_ids::Vector{Int}
end

function find_halo_pulls(body_chunk::BodyChunk, pd::PointDecomposition)
    halo_read_fields = reads_from_halo(body_chunk.mat)
    halo_write_fields = writes_to_halo(body_chunk.mat)

    #TODO
    return nothing
end

function find_halos_by_source(body_chunk::BodyChunk, pd::PointDecomposition)

    #TODO
    return nothing
end
