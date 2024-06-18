struct MPIHaloInfo
    point_ids::Dict{Int,Vector{Int}}
    halos_points::Dict{Int,Vector{Int}}
    localizers::Dict{Int,Dict{Int,Int}}
    hidxs_by_src::Dict{Int,Dict{Int,Vector{Int}}}
end

struct MPIBodyDataHandler{Sys,M,P,S,Bufs} <: AbstractMPIBodyDataHandler{Sys,M,P,S}
    chunk::BodyChunk{Sys,M,P,S}
    n_halo_fields::Int
    lth_exs_send::Vector{HaloExchange}
    lth_exs_recv::Vector{HaloExchange}
    htl_exs_send::Vector{HaloExchange}
    htl_exs_recv::Vector{HaloExchange}
    lth_send_bufs::Vector{Bufs}
    lth_recv_bufs::Vector{Bufs}
    htl_send_bufs::Vector{Bufs}
    htl_recv_bufs::Vector{Bufs}
    field_to_buf::Dict{Symbol,Int}
    lth_reqs::Vector{Vector{MPI.Request}}
    htl_reqs::Vector{Vector{MPI.Request}}
end

function mpi_data_handler(body::AbstractBody, solver::AbstractTimeSolver)
    point_decomp = PointDecomposition(body, mpi_nranks())
    param_spec = get_param_spec(body)
    @timeit_debug TO "chop chunk" begin
        chunk = chop_body_mpi(body, solver, point_decomp, param_spec)
    end
    @timeit_debug TO "find halo exchanges" begin
        halo_infos = get_halo_infos(chunk, point_decomp, body.n_points)
        lth_exs, htl_exs = find_halo_exchanges(halo_infos)
        lth_send, lth_recv = filter_send_recv(lth_exs)
        htl_send, htl_recv = filter_send_recv(htl_exs)
        bufs = find_bufs(chunk.storage, lth_send, lth_recv, htl_send, htl_recv)
        n_halo_fields = get_n_halo_fields(chunk.storage)
        lth_reqs = init_reqs(lth_send, n_halo_fields)
        htl_reqs = init_reqs(htl_send, n_halo_fields)
    end
    mdh = MPIBodyDataHandler(chunk, n_halo_fields, lth_send, lth_recv, htl_send, htl_recv,
                             bufs..., lth_reqs, htl_reqs)
    return mdh
end

function chop_body_mpi(body::AbstractBody, solver::AbstractTimeSolver,
                       point_decomp::PointDecomposition, param_spec::AbstractParamSpec)
    chunk = BodyChunk(body, solver, point_decomp, mpi_chunk_id(), param_spec)
    apply_precracks!(chunk, body)
    apply_initial_conditions!(chunk, body)
    initialize!(chunk)
    return chunk
end

function get_halo_infos(chunk::AbstractBodyChunk, pd::PointDecomposition, n_points::Int)
    point_ids = get_all_point_ids(chunk, n_points)
    halo_points = Dict{Int,Vector{Int}}()
    localizers = Dict{Int,Dict{Int,Int}}()
    hidxs_by_src = Dict{Int,Dict{Int,UnitRange{Int}}}()
    for (chunk_id, loc_point_ids) in point_ids
        n_loc_points = length(pd.decomp[chunk_id])
        halos = loc_point_ids[n_loc_points+1:end]
        halo_sources = [pd.point_src[i] for i in halos]
        halo_points[chunk_id] = halos
        localizers[chunk_id] = find_localizer(loc_point_ids)
        hidxs_by_src[chunk_id] = get_hidxs_by_src(halo_sources, n_loc_points)
    end
    return MPIHaloInfo(point_ids, halo_points, localizers, hidxs_by_src)
end

function get_all_point_ids(chunk::AbstractBodyChunk, n_points::Int)
    loc_point_ids = zeros(Int, n_points)
    for i in eachindex(chunk.ch.point_ids)
        loc_point_ids[i] = chunk.ch.point_ids[i]
    end
    _point_ids = MPI.Allgather(loc_point_ids, mpi_comm())
    point_ids = reshape(_point_ids, n_points, mpi_nranks())
    all_point_ids = Dict{Int,Vector{Int}}()
    for chunk_id in axes(point_ids, 2)
        all_point_ids[chunk_id] = filter(x -> x > 0, point_ids[:, chunk_id])
    end
    return all_point_ids
end

function find_halo_exchanges(halo_infos::MPIHaloInfo)
    lth_exs = Vector{HaloExchange}()
    htl_exs = Vector{HaloExchange}()
    ccid = mpi_chunk_id()
    for (cid, hidxs_by_src) in halo_infos.hidxs_by_src
        for (hcid, idxs) in hidxs_by_src
            _find_exs!(lth_exs, htl_exs, halo_infos, ccid, cid, hcid, idxs)
        end
    end
    return lth_exs, htl_exs
end

function _find_exs!(lth_exs::Vector, htl_exs::Vector, hi::MPIHaloInfo, ccid::Int, cid::Int,
                    hcid::Int, idxs::Vector{Int})
    (ccid == cid || ccid == hcid) || return nothing
    hidxs = hi.point_ids[cid][idxs]
    localize!(hidxs, hi.localizers[hcid])
    push!(lth_exs, HaloExchange(cid * mpi_nranks(), hcid, cid, hidxs, idxs))
    push!(htl_exs, HaloExchange(cid * mpi_nranks(), cid, hcid, idxs, hidxs))
    return nothing
end

function find_bufs(s::AbstractStorage, lth_exs_send::V, lth_exs_recv::V, htl_exs_send::V,
                   htl_exs_recv::V) where {V<:Vector{HaloExchange}}
    lth_send = [find_bufs(s, ex) for ex in lth_exs_send]
    lth_recv = [find_bufs(s, ex) for ex in lth_exs_recv]
    htl_send = [find_bufs(s, ex) for ex in htl_exs_send]
    htl_recv = [find_bufs(s, ex) for ex in htl_exs_recv]
    field_to_buf = find_field_to_buf(s)
    return lth_send, lth_recv, htl_send, htl_recv, field_to_buf
end

@inline function find_bufs(s::AbstractStorage, ex::HaloExchange)
    return _find_bufs(s, length(ex.dest_idxs))
end

@inline function _find_bufs(s::AbstractStorage, n_halo_points::Int)
    return Tuple(find_buf(a, n_halo_points) for a in get_halo_fields(s))
end

@inline find_buf(a::Matrix{T}, n::Int) where {T} = zeros(T, size(a, 1), n)
@inline find_buf(::Vector{T}, n::Int) where {T} = zeros(T, n)

@inline function find_buf(a::VecOrMat, ex::HaloExchange)
    return find_buf(a, length(ex.dest_idxs))
end

function find_field_to_buf(s::AbstractStorage)
    fields = get_halo_fieldnames(s)
    field_to_buf = Dict{Symbol,Int}()
    for (i, field) in enumerate(fields)
        field_to_buf[field] = i
    end
    return field_to_buf
end

@inline function get_buf(bufs::T, i::Int) where {T}
    return _get_buf(bufs, Val(i))
end

@inline function _get_buf(bufs::T, ::Val{N}) where {T,N}
    return bufs[N]
end

function filter_send_recv(exs::Vector{HaloExchange})
    send_exs = Vector{HaloExchange}()
    recv_exs = Vector{HaloExchange}()
    for ex in exs
        if mpi_chunk_id() == ex.src_chunk_id
            push!(send_exs, ex)
        end
        if mpi_chunk_id() == ex.dest_chunk_id
            push!(recv_exs, ex)
        end
    end
    return send_exs, recv_exs
end

function init_reqs(exs::Vector{HaloExchange}, n_halo_fields::Int)
    n_reqs = length(exs)
    return [Vector{MPI.Request}(undef, n_reqs) for _ in 1:n_halo_fields]
end

@inline function get_req_id(ex_id::Int, n_halo_fields::Int, field_id::Int)
    return ex_id * n_halo_fields + field_id - n_halo_fields
end

function get_htl_reqests(dh::MPIBodyDataHandler, field::Symbol)
    field_id = dh.field_to_buf[field]
    return dh.htl_reqs[field_id]
end

function get_lth_reqests(dh::MPIBodyDataHandler, field::Symbol)
    field_id = dh.field_to_buf[field]
    return dh.lth_reqs[field_id]
end

function exchange_loc_to_halo!(dh::MPIBodyDataHandler)
    fields = loc_to_halo_fields(dh.chunk.storage)
    isempty(fields) && return nothing
    exchange_loc_to_halo!(dh, fields)
    return nothing
end

function exchange_loc_to_halo!(dh::MPIBodyDataHandler, fields::NTuple{N,Symbol}) where {N}
    for field in fields
        exchange_loc_to_halo!(dh, field)
    end
    return nothing
end

function exchange_loc_to_halo!(dh::MPIBodyDataHandler, field::Symbol)
    for (ex_id, ex) in enumerate(dh.lth_exs_send)
        send_lth!(dh, ex, ex_id, field)
    end
    for (ex_id, ex) in enumerate(dh.lth_exs_recv)
        recv_lth!(dh, ex, ex_id, field)
    end
    MPI.Waitall(get_lth_reqests(dh, field))
    return nothing
end

function exchange_loc_to_halo!(get_field_function::F,
                               dh::MPIBodyDataHandler) where {F<:Function}
    reqs = Vector{MPI.Request}(undef, length(dh.lth_exs_send))
    for (ex_id, ex) in enumerate(dh.lth_exs_send)
        reqs[ex_id] = send_lth(get_field_function, dh, ex)
    end
    for ex in dh.lth_exs_recv
        recv_lth!(get_field_function, dh, ex)
    end
    MPI.Waitall(reqs)
    return nothing
end

function send_lth!(dh::MPIBodyDataHandler, ex::HaloExchange, ex_id::Int, field::Symbol)
    src_field = get_point_data(dh.chunk.storage, field)
    field_id = dh.field_to_buf[field]
    buf = get_buf(dh.lth_send_bufs[ex_id], field_id)
    exchange_to_buf!(buf, src_field, ex.src_idxs)
    tag = ex.tag + field_id
    dest = ex.dest_chunk_id - 1
    req = MPI.Isend(buf, mpi_comm(); dest=dest, tag=tag)
    dh.lth_reqs[field_id][ex_id] = req
    return nothing
end

function send_lth(get_field_function::F, dh::MPIBodyDataHandler,
                  ex::HaloExchange) where {F<:Function}
    src_field = get_field_function(dh.chunk)
    buf = find_buf(src_field, ex)
    exchange_to_buf!(buf, src_field, ex.src_idxs)
    dest = ex.dest_chunk_id - 1
    req = MPI.Isend(buf, mpi_comm(); dest=dest, tag=ex.tag)
    return req
end

function recv_lth!(dh::MPIBodyDataHandler, ex::HaloExchange, ex_id::Int, field::Symbol)
    field_id = dh.field_to_buf[field]
    buf = get_buf(dh.lth_recv_bufs[ex_id], field_id)
    tag = ex.tag + field_id
    source = ex.src_chunk_id - 1
    MPI.Recv!(buf, mpi_comm(); source=source, tag=tag)
    dest_field = get_point_data(dh.chunk.storage, field)
    exchange_from_buf!(dest_field, buf, ex.dest_idxs)
    return nothing
end

function recv_lth!(get_field_function::F, dh::MPIBodyDataHandler,
                   ex::HaloExchange) where {F<:Function}
    dest_field = get_field_function(dh.chunk)
    buf = find_buf(dest_field, ex)
    source = ex.src_chunk_id - 1
    MPI.Recv!(buf, mpi_comm(); source=source, tag=ex.tag)
    exchange_from_buf!(dest_field, buf, ex.dest_idxs)
    return nothing
end

function exchange_halo_to_loc!(dh::MPIBodyDataHandler)
    fields = halo_to_loc_fields(dh.chunk.storage)
    isempty(fields) && return nothing
    exchange_halo_to_loc!(dh, fields)
    return nothing
end

function exchange_halo_to_loc!(dh::MPIBodyDataHandler, fields::NTuple{N,Symbol}) where {N}
    for field in fields
        exchange_halo_to_loc!(dh, field)
    end
    return nothing
end

function exchange_halo_to_loc!(dh::MPIBodyDataHandler, field::Symbol)
    for (ex_id, ex) in enumerate(dh.htl_exs_send)
        send_htl!(dh, ex, ex_id, field)
    end
    for (ex_id, ex) in enumerate(dh.htl_exs_recv)
        recv_htl!(dh, ex, ex_id, field)
    end
    MPI.Waitall(get_htl_reqests(dh, field))
    return nothing
end

function exchange_halo_to_loc!(get_field_function::F,
                               dh::MPIBodyDataHandler) where {F<:Function}
    reqs = Vector{MPI.Request}(undef, length(dh.htl_exs_send))
    for (ex_id, ex) in enumerate(dh.htl_exs_send)
        reqs[ex_id] = send_htl(get_field_function, dh, ex)
    end
    for ex in dh.htl_exs_recv
        recv_htl!(get_field_function, dh, ex)
    end
    MPI.Waitall(reqs)
    return nothing
end

function send_htl!(dh::MPIBodyDataHandler, ex::HaloExchange, ex_id::Int, field::Symbol)
    src_field = get_point_data(dh.chunk.storage, field)
    field_id = dh.field_to_buf[field]
    buf = get_buf(dh.htl_send_bufs[ex_id], field_id)
    exchange_to_buf!(buf, src_field, ex.src_idxs)
    tag = ex.tag + field_id
    dest = ex.dest_chunk_id - 1
    req = MPI.Isend(buf, mpi_comm(); dest=dest, tag=tag)
    dh.htl_reqs[field_id][ex_id] = req
    return nothing
end

function send_htl(get_field_function::F, dh::MPIBodyDataHandler,
                  ex::HaloExchange) where {F<:Function}
    src_field = get_field_function(dh.chunk)
    buf = find_buf(src_field, ex)
    exchange_to_buf!(buf, src_field, ex.src_idxs)
    dest = ex.dest_chunk_id - 1
    req = MPI.Isend(buf, mpi_comm(); dest=dest, tag=ex.tag)
    return req
end

function recv_htl!(dh::MPIBodyDataHandler, ex::HaloExchange, ex_id::Int, field::Symbol)
    field_id = dh.field_to_buf[field]
    buf = get_buf(dh.htl_recv_bufs[ex_id], field_id)
    tag = ex.tag + field_id
    source = ex.src_chunk_id - 1
    MPI.Recv!(buf, mpi_comm(); source=source, tag=tag)
    dest_field = get_point_data(dh.chunk.storage, field)
    exchange_from_buf_add!(dest_field, buf, ex.dest_idxs)
    return nothing
end

function recv_htl!(get_field_function::F, dh::MPIBodyDataHandler,
                   ex::HaloExchange) where {F<:Function}
    dest_field = get_field_function(dh.chunk)
    buf = find_buf(dest_field, ex)
    source = ex.src_chunk_id - 1
    MPI.Recv!(buf, mpi_comm(); source=source, tag=ex.tag)
    exchange_from_buf_add!(dest_field, buf, ex.dest_idxs)
    return nothing
end

function export_results(dh::MPIBodyDataHandler, options::AbstractOptions, n::Int,
                        t::Float64; prefix="")
    options.exportflag || return nothing
    if mod(n, options.freq) == 0
        _export_results(options, dh.chunk, mpi_chunk_id(), mpi_nranks(), prefix, n, t)
    end
    return nothing
end

function export_reference_results(dh::MPIBodyDataHandler, options::AbstractOptions;
                                  prefix="")
    options.exportflag || return nothing
    _export_results(options, dh.chunk, mpi_chunk_id(), mpi_nranks(), prefix, 0, 0.0)
    return nothing
end

function initialize!(::AbstractMPIBodyDataHandler, ::AbstractTimeSolver)
    return nothing
end

function log_data_handler(options::AbstractOptions,
                          dh::AbstractMPIBodyDataHandler{Sys}) where {Sys<:BondSystem}
    msg = "BOND SYSTEM\n"

    n_bonds = MPI.Reduce(length(dh.chunk.system.bonds), MPI.SUM, mpi_comm())
    if mpi_isroot()
        msg *= log_qty("number of bonds", n_bonds)
    end
    log_it(options, msg)
    return nothing
end
