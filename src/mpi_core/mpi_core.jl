struct MPIHaloExchangeBuffer
    field::Symbol
    source::Int
    dest::Int
    halo_i::Vector{Int}
    buf::Matrix{Float64}
end

struct MPIDataHandler{M<:AbstractMaterial,S<:AbstractStorage,D<:AbstractDiscretization}
    mat::M
    loc_points::UnitRange{Int}
    n_loc_points::Int
    n_halo_points::Int
    g2l::Dict{Int,Int}
    hebs::Vector{MPIHaloExchangeBuffer}
    d::D
    s::S
end

# struct MPISimParameters{M<:AbstractMaterial} <: AbstractSimParameters
#     mat::M
#     pm::MPIManager
#     bcm::BCManager
#     refposition::Matrix{Float64}
#     volume::Vector{Float64}
#     bonds::Vector{Bond}
#     bond_range::Vector{UnitRange{Int}}
#     n_family_members::Vector{Int}
# end

# function MPISimParameters(job::J) where {J<:AbstractJob}
#     point_dist = defaultdist(job.pc.n_points, mpi_nranks())
#     rank_of_point = Dict{Int,Int}()
#     for jrid in eachindex(point_dist)
#         rank_id = jrid - 1
#         for point_id in point_dist[jrid]
#             rank_of_point[point_id] = rank_id
#         end
#     end
#     # rank local points
#     rl_points = point_dist[mpi_rank()+1]
#     not_rl_points = filter(x -> !in(x, rl_points), 1:pc.n_points)
#     n_local_points = length(rl_points)
#     # Dicts with the conversion of local index ⇆ global index
#     global_g2l = Dict{Int,Dict{Int,Int}}()
#     global_l2g = Dict{Int,Dict{Int,Int}}()
#     for (i,local_rl_points) in enumerate(point_dist)
#         rank_id = i - 1
#         local_g2l = Dict{Int,Int}()
#         local_l2g = Dict{Int,Int}()
#         point_cnt = 0
#         for gi in local_rl_points
#             point_cnt += 1
#             local_g2l[gi] = point_cnt
#             local_l2g[point_cnt] = gi
#         end
#         global_g2l[rank_id] = local_g2l
#         global_l2g[rank_id] = local_l2g
#     end
#     g2l = global_g2l[rank]
#     l2g = global_l2g[rank]
#     # neighbor search for each point
#     bonds, n_family_members= find_bonds_mpi(pc, mat, rl_points)
#     n_bonds = length(bonds)
#     hood_range = find_hood_range(n_family_members, n_local_points)

#     # find halo points
#     _halo_points = filter(x -> !in(x, rl_points), unique(bonds))
#     n_halo_points = length(_halo_points)
#     _rank_of_halo_points = [rank_of_point[i] for i in _halo_points]
#     hi_sorted = sortperm(_rank_of_halo_points)
#     halo_points = _halo_points[hi_sorted]
#     rank_of_halo_points = _rank_of_halo_points[hi_sorted]
#     ghost_points = filter(x -> !in(x, halo_points), not_rl_points)
#     # @info "rank: $rank" n_halo_points
#     # @info "rank: $rank" halo_points

#     local_pointmap = zeros(Int, pc.n_points)
#     local_pointmap[1:n_local_points] .= rl_points
#     local_pointmap[(n_local_points+1):(n_local_points+n_halo_points)] .= halo_points
#     # @info "rank: $rank" local_pointmap
#     _global_pointmap = MPI.Allgather(local_pointmap, mpi_comm())
#     global_pointmap = reshape(_global_pointmap, pc.n_points, mpi_nranks())
#     # @info "rank: $rank" global_pointmap

#     # add halos to conversion Dicts
#     # for (hi, gi) in enumerate(halo_points)
#     #     li = n_local_points + hi
#     #     g2l[gi] = li
#     #     l2g[li] = gi
#     # end

#     # Dicts with the conversion of local index ⇆ global index
#     all_halo_points = Dict{Int,Vector{Int}}()
#     for (i,local_rl_points) in enumerate(point_dist)
#         rank_id = i - 1
#         n_rl_points = length(local_rl_points)
#         _halo_points = Vector{Int}()
#         _pointmap = global_pointmap[:,i]
#         point_cnt = n_rl_points
#         for gi in _pointmap[n_rl_points+1:end]
#             if gi != 0
#                 point_cnt += 1
#                 global_g2l[rank_id][gi] = point_cnt
#                 global_l2g[rank_id][point_cnt] = gi
#                 push!(_halo_points, gi)
#             else
#                 break
#             end
#         end
#         all_halo_points[rank_id] = _halo_points
#     end
#     g2l = global_g2l[rank]
#     l2g = global_l2g[rank]

#     all_halos_by_src_rank = Dict{Int, Dict{Int, UnitRange{Int}}}()
#     for (rank_id, halos) in all_halo_points
#         rank_of_halos = [rank_of_point[i] for i in halos]
#         @assert rank_of_halos == sort(rank_of_halos)
#         all_halos_by_src_rank[rank_id] = get_halos_by_rank(halos, rank_of_halos)
#     end

#     hebs = Vector{HaloExchangeBuf}()
#     for (parent_rank, halos_by_src_rank) in all_halos_by_src_rank
#         for (src_rank, halo_ids) in halos_by_src_rank
#             halos_gi = all_halo_points[parent_rank][halo_ids]
#             buf = zeros(3, length(halos_gi))
#             heb = HaloExchangeBuf(src_rank, parent_rank, halos_gi, buf)
#             push!(hebs, heb)
#         end
#     end

#     cells = get_cells(n_local_points)
#     lbcs = [localize_bc(bc, rl_points) for bc in bcs]
#     # println("\nrank $(mpi_rank()) lbcs:\n$(lbcs)\n")
#     sp = SimulationParameters(rank, com_size, n_local_points, rl_points, not_rl_points,
#                             halo_points, ghost_points, global_g2l, global_l2g, g2l, l2g,
#                             all_halo_points, all_halos_by_src_rank,
#                             mat, pc, cells, n_bonds, bonds, init_dists, failure_allowed,
#                             n_family_members, hood_range, lbcs, hebs, n_timesteps,
#                             export_freq, export_path)

#                             #TODO
# end

function _init_sim_mpi(job::J) where {J<:AbstractJob}
end

function _force_density!(dh::MPIDataHandler, mat::AbstractMaterial)
    dh.s.b_int .= 0
    dh.s.n_active_family_members .= 0
    for i in dh.loc_points
        force_density!(dh.s, d, mat, i)
    end
end
