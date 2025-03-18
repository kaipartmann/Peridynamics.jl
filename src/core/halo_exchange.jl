"""
    HaloExchange

$(internal_api_warning())

A type used for communication between body chunks in MPI.

# Fields

- `tag::Int`: Tag used for the MPI sending and receiving commands
- `src_chunk_id::Int`: Index of the chunk that sends information
- `dest_chunk_id::Int`: Index of the chunk that receives information
- `src_idxs::Vector{Int}`: Indices of the points in the source chunk that send information
- `dest_idxs::Vector{Int}`: Indices of the points in the destination chunk that receive
    information
"""
struct HaloExchange
    tag::Int
    src_chunk_id::Int
    dest_chunk_id::Int
    src_idxs::Vector{Int}
    dest_idxs::Vector{Int}
end

function exchange!(dest::Matrix{T}, src::Matrix{T}, dest_idxs::Vector{Int},
                   src_idxs::Vector{Int}) where {T}
    for i in eachindex(dest_idxs)
        for d in axes(dest, 1)
            @inbounds dest[d, dest_idxs[i]] = src[d, src_idxs[i]]
        end
    end
    return nothing
end

function exchange_from_buf!(dest::T, buf::T, dest_idxs::Vector{Int}) where {T<:Matrix}
    for i in eachindex(dest_idxs)
        for d in axes(dest, 1)
            dest[d, dest_idxs[i]] = buf[d, i]
        end
    end
    return nothing
end

function exchange_to_buf!(buf::T, src::T, src_idxs::Vector{Int}) where {T<:Matrix}
    for i in eachindex(src_idxs)
        for d in axes(buf, 1)
            buf[d, i] = src[d, src_idxs[i]]
        end
    end
    return nothing
end

function exchange!(dest::Vector{T}, src::Vector{T}, dest_idxs::Vector{Int},
                   src_idxs::Vector{Int}) where {T}
    for i in eachindex(dest_idxs)
        @inbounds dest[dest_idxs[i]] = src[src_idxs[i]]
    end
    return nothing
end

function exchange_from_buf!(dest::T, buf::T, dest_idxs::Vector{Int}) where {T<:Vector}
    for i in eachindex(dest_idxs)
        dest[dest_idxs[i]] = buf[i]
    end
    return nothing
end

function exchange_to_buf!(buf::T, src::T, src_idxs::Vector{Int}) where {T<:Vector}
    for i in eachindex(src_idxs)
        buf[i] = src[src_idxs[i]]
    end
    return nothing
end

function exchange_add!(dest::Matrix{T}, src::Matrix{T}, dest_idxs::Vector{Int},
                       src_idxs::Vector{Int}) where {T<:Number}
    for i in eachindex(dest_idxs)
        for d in axes(dest, 1)
            @inbounds dest[d, dest_idxs[i]] += src[d, src_idxs[i]]
        end
    end
    return nothing
end

function exchange_from_buf_add!(dest::T, buf::T, dest_idxs::Vector{Int}) where {T<:Matrix}
    for i in eachindex(dest_idxs)
        for d in axes(dest, 1)
            dest[d, dest_idxs[i]] += buf[d, i]
        end
    end
    return nothing
end

function exchange_add!(dest::Vector{T}, src::Vector{T}, dest_idxs::Vector{Int},
                       src_idxs::Vector{Int}) where {T<:Number}
    for i in eachindex(dest_idxs)
        @inbounds dest[dest_idxs[i]] += src[src_idxs[i]]
    end
    return nothing
end

function exchange_from_buf_add!(dest::T, buf::T, dest_idxs::Vector{Int}) where {T<:Vector}
    for i in eachindex(dest_idxs)
        dest[dest_idxs[i]] += buf[i]
    end
    return nothing
end
