const Bond = Tuple{Int, Int, Float64, Bool}

function find_bonds(pc::PointCloud, mat::PDMaterial, owned_points::Vector{UnitRange{Int}})
    n_threads = nthreads()
    _bond_data = fill([(0, 0, 0.0, true)], n_threads)
    n_family_members = zeros(Int, pc.n_points)
    p = Progress(pc.n_points; dt = 1, desc = "Search bonds...     ", barlen = 30,
                 color = :normal, enabled = !is_logging(stderr))
    @threads for tid in 1:n_threads
        local_bond_data = Vector{Tuple{Int, Int, Float64, Bool}}(undef, 0)
        sizehint!(local_bond_data, pc.n_points * 500)
        idst = 0.0
        num = 0
        fail = false
        for a in owned_points[tid]
            num = 0
            δ = mat[a].δ
            for i in 1:pc.n_points
                if a !== i
                    idst = sqrt((pc.position[1, i] - pc.position[1, a])^2 +
                                (pc.position[2, i] - pc.position[2, a])^2 +
                                (pc.position[3, i] - pc.position[3, a])^2)
                    if idst <= δ
                        num += 1
                        fail = pc.failure_flag[a] & pc.failure_flag[i]
                        push!(local_bond_data, (a, i, idst, fail))
                    end
                end
            end
            n_family_members[a] = num
            next!(p)
        end
        _bond_data[tid] = local_bond_data
    end
    finish!(p)
    bond_data = reduce(append!, _bond_data)
    return bond_data, n_family_members
end

#TODO!!!
function find_bonds!(bond_buf::Vector{Bond}, i::Int, pc::PointCloud, δ::Float64)
    for j in 1:pc.n_points
        if a !== i
            idst = sqrt((pc.position[1, j] - pc.position[1, i])^2 +
                        (pc.position[2, j] - pc.position[2, i])^2 +
                        (pc.position[3, j] - pc.position[3, i])^2)
            if idst <= δ
                num += 1
                fail = pc.failure_flag[a] & pc.failure_flag[i]
                push!(local_bond_data, (a, i, idst, fail))
            end
        end
    end
end
