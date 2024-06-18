struct ThreadsMultibodyDataHandler{Sys,M,P,S} <: AbstractThreadsMultibodyDataHandler{Sys,M,P,S}
    n_bodies::Int
    body_dhs::Vector{ThreadsBodyDataHandler{Sys,M,P,S}}
    srf_contacts::Vector{ShortRangeForceContact}
end

function threads_data_handler(ms::AbstractMultibodySetup, solver::AbstractTimeSolver,
                              n_chunks::Int)
    body_dhs = [threads_data_handler(body, solver, n_chunks) for body in each_body(ms)]
    n_bodies = length(body_dhs)
    return ThreadsMultibodyDataHandler(n_bodies, body_dhs, ms.srf_contacts)
end

@inline each_body_dh(dh::ThreadsMultibodyDataHandler) = dh.body_dhs
