struct ThreadsDataHandler{C<:AbstractBodyChunk} <: AbstractDataHandler
    n_chunks::Int
    chunks::Vector{C}
    halo_exchanges::Vector{Vector{HaloExchange}}
end

function ThreadsDataHandler(body::Body, time_solver::AbstractTimeSolver,
                            point_decomp::PointDecomposition)
    n_param = length(body.point_params)
    v = length(n_param) == 1 ? Val{1}() : Val{2}()
    return ThreadsDataHandler(body, time_solver, point_decomp, v)
end

function ThreadsDataHandler(body::Body, time_solver::AbstractTimeSolver,
                            point_decomp::PointDecomposition, v::Val{N}) where {N}
    chunks = chop_body_threads(body, time_solver, point_decomp, v)
    n_chunks = length(chunks)
    _halo_exchanges = find_halo_exchanges(chunks)
    halo_exchanges = [Vector{HaloExchange}() for _ in eachindex(chunks)]
    @threads :static for chunk_id in eachindex(chunks)
        halo_exchanges[chunk_id] = filter(x -> x.dest_chunk_id == chunk_id, _halo_exchanges)
    end
    return ThreadsDataHandler(n_chunks, chunks, halo_exchanges)
end

function ThreadsDataHandler(multibody::MultibodySetup, time_solver::AbstractTimeSolver,
                            point_decomp::PointDecomposition)
    error("MultibodySetup not yet implemented!\n")
end

function calc_stable_timestep(dh::ThreadsDataHandler, safety_factor::Float64)
    Δt = zeros(length(dh.chunks))
    @threads :static for chunk_id in eachindex(dh.chunks)
        Δt[chunk_id] = calc_timestep(dh.chunks[chunk_id])
    end
    return minimum(Δt) * safety_factor
end

get_cells(n::Int) = [MeshCell(VTKCellTypes.VTK_VERTEX, (i,)) for i in 1:n]

function halo_exchange!(dh::ThreadsDataHandler, chunk_id::Int)
    for he in dh.halo_exchanges[chunk_id]
        src_field = get_exchange_field(dh.chunks[he.src_chunk_id], he.field)
        dest_field = get_exchange_field(dh.chunks[he.dest_chunk_id], he.field)
        exchange!(dest_field, src_field, he.dest_idxs, he.src_idxs)
    end
    return nothing
end

@inline function get_exchange_field(b::AbstractBodyChunk, fieldname::Symbol)
    return getfield(b.store, fieldname)::Matrix{Float64}
end

function exchange!(dest::Matrix{Float64}, src::Matrix{Float64}, dest_idxs::Vector{Int},
                   src_idxs::Vector{Int})
    for i in eachindex(dest_idxs, src_idxs)
        dest[1, dest_idxs[i]] = src[1, src_idxs[i]]
        dest[2, dest_idxs[i]] = src[2, src_idxs[i]]
        dest[3, dest_idxs[i]] = src[3, src_idxs[i]]
    end
    return nothing
end
