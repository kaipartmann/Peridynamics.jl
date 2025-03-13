"""
    PointDecomposition

$(internal_api_warning())

A type that describes how a body is divided into multiple body chunks

# Fields

- `n_chunks::Int`: Number of body chunks
- `decomp::Vector{UnitRange{Int}}`: Indices of the points belonging to each chunk
- `point_src::Dict{Int,Int}`: Dict that assigns all point indices to the chunk they belong
    to
"""
struct PointDecomposition
    n_chunks::Int
    decomp::Vector{UnitRange{Int}}
    point_src::Dict{Int,Int}

    function PointDecomposition(decomp::Vector{UnitRange{Int}})
        n_chunks = length(decomp)
        point_src = Dict{Int,Int}()
        for chunk_id in eachindex(decomp)
            for point_id in decomp[chunk_id]
                point_src[point_id] = chunk_id
            end
        end
        return new(n_chunks, decomp, point_src)
    end
end

function distribute_equally(n_elems::Int, n_chunks::Int)
    if n_elems >= n_chunks
        chunk_size = div(n_elems, n_chunks)
        remainder = rem(n_elems, n_chunks)
        sidx = zeros(Int64, n_chunks + 1)
        for i in 1:(n_chunks + 1)
            sidx[i] += (i - 1) * chunk_size + 1
            if i <= remainder
                sidx[i] += i - 1
            else
                sidx[i] += remainder
            end
        end
        decomp = fill(0:0, n_chunks)
        for i in 1:n_chunks
            decomp[i] = sidx[i]:(sidx[i + 1] - 1)
        end
        return decomp
    else
        sidx = [1:(n_elems + 1);]
        decomp = fill(0:-1, n_chunks)
        for i in 1:n_elems
            decomp[i] = sidx[i]:(sidx[i + 1] - 1)
        end
        return decomp
    end
end

function PointDecomposition(b::AbstractBody, n_chunks::Int)
    decomp = distribute_equally(b.n_points, n_chunks)
    return PointDecomposition(decomp)
end

function PointDecomposition(::AbstractMultibodySetup, ::Int)
    error("spatial decomposition for MultibodySetup not yet implemented!\n")
    return nothing
end
