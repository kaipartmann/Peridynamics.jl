function defaultdist(sz::Int, nc::Int)
    if sz >= nc
        chunk_size = div(sz, nc)
        remainder = rem(sz, nc)
        sidx = zeros(Int64, nc + 1)
        for i in 1:(nc + 1)
            sidx[i] += (i - 1) * chunk_size + 1
            if i <= remainder
                sidx[i] += i - 1
            else
                sidx[i] += remainder
            end
        end
        grid = fill(0:0, nc)
        for i in 1:nc
            grid[i] = sidx[i]:(sidx[i + 1] - 1)
        end
        return grid
    else
        sidx = [1:(sz + 1);]
        grid = fill(0:-1, nc)
        for i in 1:sz
            grid[i] = sidx[i]:(sidx[i + 1] - 1)
        end
        return grid
    end
end

@timeit TO function find_tids(sum_tids::Vector{Vector{Int}})
    single_tids_ids = findall(x -> length(x) == 1 && x[1] != 1, sum_tids)
    single_tids = Vector{Tuple{Int, Int}}(undef, length(single_tids_ids))
    for i in eachindex(single_tids_ids)
        single_tids[i] = (single_tids_ids[i], sum_tids[single_tids_ids[i]][1])
    end
    multi_tids_ids = findall(x -> length(x) > 1, sum_tids)
    multi_tids = Vector{Tuple{Int, Vector{Int}}}(undef, length(multi_tids_ids))
    for i in eachindex(multi_tids_ids)
        point_id = multi_tids_ids[i]
        tids = filter(x -> x != 1, sum_tids[point_id])
        multi_tids[i] = (point_id, tids)
    end
    return single_tids, multi_tids
end

@timeit TO function update_thread_cache!(body::AbstractPDBody)
    @threads for (i, tid) in body.single_tids
        body.b_int[1, i, 1] = body.b_int[1, i, tid]
        body.b_int[2, i, 1] = body.b_int[2, i, tid]
        body.b_int[3, i, 1] = body.b_int[3, i, tid]
        body.n_active_family_members[i, 1] = body.n_active_family_members[i, tid]
    end
    @threads for (i, tids) in body.multi_tids
        for tid in tids
            body.b_int[1, i, 1] += body.b_int[1, i, tid]
            body.b_int[2, i, 1] += body.b_int[2, i, tid]
            body.b_int[3, i, 1] += body.b_int[3, i, tid]
            body.n_active_family_members[i, 1] += body.n_active_family_members[i, tid]
        end
    end
end

is_logging(io) = isa(io, Base.TTY) == false || (get(ENV, "CI", nothing) == "true")
