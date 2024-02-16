function default_decomp(n_points::Int, n_chunks::Int)
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

function point_decomposition(b::Body, n_chunks::Int)
    return default_decomp(b.n_points, n_chunks)
end

function point_decomposition(::MultibodySetup, ::Int)
    error("spatial decomposition for MultibodySetup not yet implemented!\n")
    return nothing
end
