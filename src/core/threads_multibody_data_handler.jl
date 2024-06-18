struct ThreadsMultibodyDataHandler{Sys,M,P,S} <: AbstractThreadsMultibodyDataHandler{Sys,M,P,S}
    n_bodies::Int
    body_dhs::Vector{ThreadsBodyDataHandler{Sys,M,P,S}}
    body_names::Vector{Symbol}
    body_idxs::Dict{Symbol,Int}
    srf_contacts::Vector{ShortRangeForceContact}
end

function threads_data_handler(ms::AbstractMultibodySetup, solver::AbstractTimeSolver,
                              n_chunks::Int)
    body_dhs = [threads_data_handler(body, solver, n_chunks) for body in each_body(ms)]
    n_bodies = length(body_dhs)
    dh = ThreadsMultibodyDataHandler(n_bodies, body_dhs, ms.body_names, ms.body_idxs, ms.srf_contacts)
    return dh
end

@inline function get_body_dh(dh::ThreadsMultibodyDataHandler, name::Symbol)
    return dh.body_dhs[dh.body_idxs[name]]
end
@inline get_body_dh(dh::ThreadsMultibodyDataHandler, idx::Int) = dh.body_dhs[idx]

@inline function get_body_name(dh::ThreadsMultibodyDataHandler, idx::Int)
    return string(dh.body_names[idx])
end

@inline each_body_dh(dh::ThreadsMultibodyDataHandler) = dh.body_dhs
@inline each_body_name(dh::ThreadsMultibodyDataHandler) = dh.body_names
@inline each_body_idx(dh::ThreadsMultibodyDataHandler) = eachindex(dh.body_dhs)

function export_reference_results(dh::ThreadsMultibodyDataHandler, options::AbstractOptions)
    options.exportflag || return nothing
    for body_idx in each_body_idx(dh)
        body_dh = get_body_dh(dh, body_idx)
        name = get_body_name(dh, body_idx)
        export_reference_results(body_dh, options; prefix=name)
    end
    return nothing
end
