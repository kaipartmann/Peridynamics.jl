struct ThreadsDataHandler{C<:AbstractBodyChunk}
    chunks::Vector{C}
    halo_exchanges::Vector{Vector{HaloExchange}}
end

function ThreadsDataHandler(body::Body, time_solver::AbstractTimeSolver,
                            point_decomp::PointDecomposition)
    body_chunks = chop_body_threads(body, time_solver, point_decomp)
    _halo_exchanges = find_halo_exchanges(body_chunks)
    halo_exchanges = [Vector{HaloExchange}() for _ in eachindex(body_chunks)]
    @threads :static for chunk_id in eachindex(body_chunks)
        halo_exchanges[chunk_id] = filter(x -> x.to_chunk_id == chunk_id, _halo_exchanges)
    end
    return ThreadsDataHandler(body_chunks, halo_exchanges)
end

function calc_stable_timestep(dh::ThreadsDataHandler, safety_factor::Float64)
    Δt = zeros(length(dh.chunks))
    @threads :static for chunk_id in eachindex(dh.chunks)
        Δt[chunk_id] = calc_timestep(dh.chunks[chunk_id])
    end
    return minimum(Δt) * safety_factor
end
