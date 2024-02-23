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

function distribute_points(n_points::Int, n_chunks::Int)
    if n_points >= n_chunks
        chunk_size = div(n_points, n_chunks)
        remainder = rem(n_points, n_chunks)
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
        sidx = [1:(n_points + 1);]
        decomp = fill(0:-1, n_chunks)
        for i in 1:n_points
            decomp[i] = sidx[i]:(sidx[i + 1] - 1)
        end
        return decomp
    end
end

function PointDecomposition(b::Body, n_chunks::Int)
    decomp = distribute_points(b.n_points, n_chunks)
    return PointDecomposition(decomp)
end

function PointDecomposition(::MultibodySetup, ::Int)
    error("spatial decomposition for MultibodySetup not yet implemented!\n")
    return nothing
end
