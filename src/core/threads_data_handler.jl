struct ThreadsDataHandler{C<:AbstractBodyChunk}
    chunks::Vector{C}
    halo_exchanges::Vector{HaloExchange}
end

function ThreadsDataHandler(body::Body, time_solver::AbstractTimeSolver,
                            point_decomp::PointDecomposition)
    body_chunks = chop_body_threads(body, time_solver, point_decomp)
    halo_exchanges = find_halo_exchanges(body_chunks)
    return ThreadsDataHandler(body_chunks, halo_exchanges)
end
