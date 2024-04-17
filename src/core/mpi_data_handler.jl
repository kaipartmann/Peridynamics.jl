struct HaloExchangeBuf{T,K}
    n_fields::Int
    src_chunk_id::Int
    dest_chunk_id::Int
    src_idxs::Vector{Int}
    dest_idxs::Vector{Int}
    bufs::T
    tags::K
    function HaloExchangeBuf(n::Int, src_chunk_id::Int, dest_chunk_id::Int,
                             src_idxs::AbstractVector{Int}, dest_idxs::AbstractVector{Int},
                             bufs::T, tags::K) where {T,K}
        return new{T,K}(n, src_chunk_id, dest_chunk_id, src_idxs, dest_idxs, bufs, tags)
    end
end

struct MPIHaloInfo
    point_ids::Dict{Int,Vector{Int}}
    halos_points::Dict{Int,Vector{Int}}
    localizers::Dict{Int,Dict{Int,Int}}
    hidxs_by_src::Dict{Int,Dict{Int,Vector{Int}}}
end

struct MPIDataHandler{C<:AbstractBodyChunk,RB,WB} <: AbstractMPIDataHandler
    chunk::C
    lth_exs_send::Vector{HaloExchangeBuf{RB}}
    lth_exs_recv::Vector{HaloExchangeBuf{RB}}
    htl_exs_send::Vector{HaloExchangeBuf{WB}}
    htl_exs_recv::Vector{HaloExchangeBuf{WB}}
    lth_reqs::Vector{MPI.Request}
    htl_reqs::Vector{MPI.Request}
end

function MPIDataHandler(body::AbstractBody, time_solver::AbstractTimeSolver,
                        point_decomp::PointDecomposition)
    n_param = length(body.point_params)
    v = n_param == 1 ? Val{1}() : Val{2}()
    return MPIDataHandler(body, time_solver, point_decomp, v)
end

function MPIDataHandler(body::AbstractBody, time_solver::AbstractTimeSolver,
                        point_decomp::PointDecomposition, v::Val{N}) where {N}
    chunk = chop_body_mpi(body, time_solver, point_decomp, v)
    @timeit_debug TO "find halo exchanges" begin
        halo_infos = get_halo_infos(chunk, point_decomp, body.n_points)
        lth_exs, htl_exs = find_halo_exchange_bufs(chunk, halo_infos)
        lth_exs_send, lth_exs_recv = filter_send_recv(lth_exs)
        htl_exs_send, htl_exs_recv = filter_send_recv(htl_exs)
        lth_reqs, htl_reqs = init_reqs(lth_exs_send), init_reqs(htl_exs_send)
    end
    return MPIDataHandler(chunk, lth_exs_send, lth_exs_recv, htl_exs_send,
                          htl_exs_recv, lth_reqs, htl_reqs)
end

function MPIDataHandler(multibody::AbstractMultibodySetup,
                        time_solver::AbstractTimeSolver,
                        point_decomp::PointDecomposition)
    error("MultibodySetup not yet implemented!\n")
end

function chop_body_mpi(body::Body{M,P}, ts::T, pd::PointDecomposition,
                       v::Val{N}) where {M,P,T,N}
    body_chunks = _chop_body_mpi(body, ts, pd, v)
    return body_chunks
end

function _chop_body_mpi(body::Body{M,P}, ts::T, pd::PointDecomposition,
                        ::Val{1}) where {M,P,T}
    @timeit_debug TO "chop chunk" chunk = BodyChunk(body, ts, pd, mpi_chunk_id())
    init_chunk!(chunk)
    apply_precracks!(chunk, body)
    apply_initial_conditions!(chunk, body)
    return chunk
end

function _chop_body_mpi(body::Body{M,P}, ts::T, pd::PointDecomposition,
                        ::Val{N}) where {M,P,T,N}
    @timeit_debug TO "chop chunk" chunk = MultiParamBodyChunk(body, ts, pd, mpi_chunk_id())
    init_chunk!(chunk)
    apply_precracks!(chunk, body)
    apply_initial_conditions!(chunk, body)
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

function find_halo_exchange_bufs(chunk::AbstractBodyChunk, halo_infos::MPIHaloInfo)
    has_lth = has_lth_exs(chunk)
    has_htl = has_htl_exs(chunk)
    lth_exs, htl_exs = _find_halo_exchange_bufs(chunk, halo_infos, has_lth, has_htl)
    return lth_exs, htl_exs
end

function has_lth_exs(chunk::AbstractBodyChunk)
    hrfs = loc_to_halo_fields(chunk.storage)
    isempty(hrfs) && return false
    return true
end

function has_htl_exs(chunk::AbstractBodyChunk)
    hwfs = halo_to_loc_fields(chunk.storage)
    isempty(hwfs) && return false
    return true
end

function _find_halo_exchange_bufs(chunk::AbstractBodyChunk, halo_info::MPIHaloInfo,
                                  has_lth::Bool, has_htl::Bool)
    RB = get_loc_to_halo_buftype(chunk.storage)
    WB = get_halo_to_loc_buftype(chunk.storage)
    lth_exs = Vector{HaloExchangeBuf{RB}}()
    htl_exs = Vector{HaloExchangeBuf{WB}}()
    ccid = mpi_chunk_id()
    for (cid, hidxs_by_src) in halo_info.hidxs_by_src
        for (hcid, idxs) in hidxs_by_src
            _find_hexs!(lth_exs, htl_exs, chunk, halo_info, ccid, cid, hcid, idxs,
                        has_lth, has_htl)
        end
    end
    return lth_exs, htl_exs
end

@inline function get_loc_to_halo_buftype(s::AbstractStorage)
    return typeof(get_loc_to_halo_fields(s))
end

@inline function get_halo_to_loc_buftype(s::AbstractStorage)
    return typeof(get_halo_to_loc_fields(s))
end

function _find_hexs!(rhexs::Vector, whexs::Vector, chunk::AbstractBodyChunk,
                     hi::MPIHaloInfo, ccid::Int, cid::Int, hcid::Int, idxs::Vector{Int},
                     has_lth::Bool, has_htl::Bool)
    (ccid == cid || ccid == hcid) || return nothing
    hidxs = hi.point_ids[cid][idxs]
    localize!(hidxs, hi.localizers[hcid])
    if has_lth
        bufs = get_loc_to_halo_bufs(chunk.storage, length(idxs))
        tags = get_tags(bufs, cid * mpi_nranks())
        n_fields = length(bufs)
        hex = HaloExchangeBuf(n_fields, hcid, cid, hidxs, idxs, bufs, tags)
        push!(rhexs, hex)
    end
    if has_htl
        bufs = get_loc_to_halo_bufs(chunk.storage, length(idxs))
        tags = get_tags(bufs, cid * mpi_nranks())
        n_fields = length(bufs)
        hex = HaloExchangeBuf(n_fields, cid, hcid, idxs, hidxs, bufs, tags)
        push!(whexs, hex)
    end
    return nothing
end

@inline function get_loc_to_halo_bufs(s::AbstractStorage, n_halo_points::Int)
    return Tuple(get_buf(a, n_halo_points) for a in get_loc_to_halo_fields(s))
end

@inline function get_halo_to_loc_bufs(s::AbstractStorage, n_halo_points::Int)
    return Tuple(get_buf(a, n_halo_points) for a in get_halo_to_loc_fields(s))
end

@inline get_buf(a::Matrix{T}, n::Int) where {T} = zeros(T, size(a, 1), n)
@inline get_buf(::Vector{T}, n::Int) where {T} = zeros(T, n)

@inline function get_tags(bufs, i)
    return Tuple(i + j for j in eachindex(bufs))
end

function filter_send_recv(hexs::Vector{HaloExchangeBuf{B}}) where {B}
    send_hexs = Vector{HaloExchangeBuf{B}}()
    recv_hexs = Vector{HaloExchangeBuf{B}}()
    for hex in hexs
        if mpi_chunk_id() == hex.src_chunk_id
            push!(send_hexs, hex)
        end
        if mpi_chunk_id() == hex.dest_chunk_id
            push!(recv_hexs, hex)
        end
    end
    return send_hexs, recv_hexs
end

function init_reqs(hexs::Vector{HaloExchangeBuf{B}}) where {B}
    n_reqs = 0
    for hex in hexs
        if mpi_chunk_id() == hex.src_chunk_id
            n_reqs += hex.n_fields
        end
    end
    return Vector{MPI.Request}(undef, n_reqs)
end

function calc_stable_timestep(dh::MPIDataHandler, safety_factor::Float64)
    _Δt = calc_timestep(dh.chunk)
    Δt = MPI.Allreduce(_Δt, MPI.MIN, mpi_comm())
    return Δt * safety_factor
end

function exchange_loc_to_halo!(dh::MPIDataHandler)
    for (hex_id, hex) in enumerate(dh.lth_exs_send)
        send_lth!(dh, hex, hex_id)
    end
    for hex in dh.lth_exs_recv
        recv_lth!(dh, hex)
    end
    MPI.Waitall(dh.lth_reqs)
    return nothing
end

function send_lth!(dh::MPIDataHandler, he::HaloExchangeBuf, hex_id::Int)
    fields = loc_to_halo_fields(dh.chunk.storage)
    for i in eachindex(fields, he.bufs)
        src_field = get_point_data(dh.chunk.storage, fields[i])
        exchange_to_buf!(he.bufs[i], src_field, he.src_idxs)
        req = MPI.Isend(he.bufs[i], mpi_comm(); dest=he.dest_chunk_id - 1, tag=he.tags[i])
        req_id = hex_id * he.n_fields + i - he.n_fields
        dh.lth_reqs[req_id] = req
    end
    return nothing
end

function recv_lth!(dh::MPIDataHandler, he::HaloExchangeBuf)
    fields = loc_to_halo_fields(dh.chunk.storage)
    for i in eachindex(fields, he.bufs)
        MPI.Recv!(he.bufs[i], mpi_comm(); source=he.src_chunk_id - 1, tag=he.tags[i])
        dest_field = get_point_data(dh.chunk.storage, fields[i])
        exchange_from_buf!(dest_field, he.bufs[i], he.dest_idxs)
    end
    return nothing
end

function exchange_halo_to_loc!(dh::MPIDataHandler)
    for (hex_id, hex) in enumerate(dh.htl_exs_send)
        send_htl!(dh, hex, hex_id)
    end
    for hex in dh.htl_exs_recv
        recv_htl!(dh, hex)
    end
    MPI.Waitall(dh.htl_reqs)
    return nothing
end

function send_htl!(dh::MPIDataHandler, he::HaloExchangeBuf, hex_id::Int)
    fields = halo_to_loc_fields(dh.chunk.storage)
    for i in eachindex(fields, he.bufs)
        src_field = get_point_data(dh.chunk.storage, fields[i])
        exchange_to_buf!(he.bufs[i], src_field, he.src_idxs)
        req = MPI.Isend(he.bufs[i], mpi_comm(); dest=he.dest_chunk_id - 1, tag=he.tags[i])
        req_id = hex_id * he.n_fields + i - he.n_fields
        dh.htl_reqs[req_id] = req
    end
    return nothing
end

function recv_htl!(dh::MPIDataHandler, he::HaloExchangeBuf)
    fields = halo_to_loc_fields(dh.chunk.storage)
    for i in eachindex(fields, he.bufs)
        MPI.Recv!(he.bufs[i], mpi_comm(); source=he.src_chunk_id - 1, tag=he.tags[i])
        dest_field = get_point_data(dh.chunk.storage, fields[i])
        exchange_from_buf_add!(dest_field, he.bufs[i], he.dest_idxs)
    end
    return nothing
end

function export_results(dh::MPIDataHandler, options::AbstractOptions, n::Int, t::Float64)
    options.exportflag || return nothing
    if mod(n, options.freq) == 0
        _export_results(dh.chunk, mpi_chunk_id(), mpi_nranks(), options, n, t)
    end
    return nothing
end

function export_reference_results(dh::MPIDataHandler, options::AbstractOptions)
    options.exportflag || return nothing
    _export_results(dh.chunk, mpi_chunk_id(), mpi_nranks(), options, 0, 0.0)
    return nothing
end
