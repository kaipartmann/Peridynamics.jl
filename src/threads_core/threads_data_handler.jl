
struct ThreadsHaloExchange
    field::Symbol
    source::Int
    dest::Int
    point_ids::Vector{Int}
end

struct ThreadsDataHandler{M<:AbstractMaterialConfig,D<:AbstractDiscretization,
                          S<:AbstractStorage}
    mat::M
    vclh::Vector{ChunkLocalHandler{D,S}}
    vthe::Vector{ThreadsHaloExchange}
end

# TODO: dispatch on SingleBodyJob!
function init_threads_data_handler(job::AbstractJob, n_chunks::Int)
    vclh = chunk_local_handlers(job, n_chunks)

    #TODO: ThreadsDataHandler!!!
end

function _force_density!(dh::ThreadsDataHandler, mat::Material)
    @threads :static for cid in 1:dh.n_chunks
        s = dh.vs[cid]
        d = dh.vd[cid]
        s.b_int .= 0
        s.n_active_family_members .= 0
        for i in dh.loc_points
            force_density!(s, d, mat, i)
        end
    end
end
