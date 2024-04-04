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
    read_hexs_send::Vector{HaloExchangeBuf{RB}}
    read_hexs_recv::Vector{HaloExchangeBuf{RB}}
    write_hexs_send::Vector{HaloExchangeBuf{WB}}
    write_hexs_recv::Vector{HaloExchangeBuf{WB}}
    read_reqs::Vector{MPI.Request}
    write_reqs::Vector{MPI.Request}
end

function MPIDataHandler(body::Body, time_solver::AbstractTimeSolver,
                        point_decomp::PointDecomposition)
    n_param = length(body.point_params)
    v = n_param == 1 ? Val{1}() : Val{2}()
    return MPIDataHandler(body, time_solver, point_decomp, v)
end

function MPIDataHandler(body::Body, time_solver::AbstractTimeSolver,
                        point_decomp::PointDecomposition, v::Val{N}) where {N}
    chunk = chop_body_mpi(body, time_solver, point_decomp, v)
    @timeit_debug TO "find halo exchanges" begin
        halo_infos = get_halo_infos(chunk, point_decomp, body.n_points)
        read_hexs, write_hexs = find_halo_exchange_bufs(chunk, halo_infos)
        read_hexs_send, read_hexs_recv = filter_send_recv(read_hexs)
        write_hexs_send, write_hexs_recv = filter_send_recv(write_hexs)
        read_reqs, write_reqs = init_reqs(read_hexs_send), init_reqs(write_hexs_send)
    end
    return MPIDataHandler(chunk, read_hexs_send, read_hexs_recv, write_hexs_send,
                          write_hexs_recv, read_reqs, write_reqs)
end

function MPIDataHandler(multibody::MultibodySetup, time_solver::AbstractTimeSolver,
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
    apply_precracks!(chunk, body)
    apply_initial_conditions!(chunk, body)
    return chunk
end

function _chop_body_mpi(body::Body{M,P}, ts::T, pd::PointDecomposition,
                        ::Val{N}) where {M,P,T,N}
    @timeit_debug TO "chop chunk" chunk = MultiParamBodyChunk(body, ts, pd, mpi_chunk_id())
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
    has_read = has_read_halos(chunk)
    has_write = has_write_halos(chunk)
    read_hexs, write_hexs = _find_halo_exchange_bufs(chunk, halo_infos, has_read, has_write)
    return read_hexs, write_hexs
end

function has_read_halos(chunk::AbstractBodyChunk)
    hrfs = get_halo_read_fields(chunk.store)
    isempty(hrfs) && return false
    return true
end

function has_write_halos(chunk::AbstractBodyChunk)
    hwfs = get_halo_write_fields(chunk.store)
    isempty(hwfs) && return false
    return true
end

function _find_halo_exchange_bufs(chunk::AbstractBodyChunk, halo_info::MPIHaloInfo,
                                  hasread::Bool, haswrite::Bool)
    RB = get_halo_read_buftype(chunk.store)
    WB = get_halo_write_buftype(chunk.store)
    read_hexs = Vector{HaloExchangeBuf{RB}}()
    write_hexs = Vector{HaloExchangeBuf{WB}}()
    ccid = mpi_chunk_id()
    for (cid, hidxs_by_src) in halo_info.hidxs_by_src
        for (hcid, idxs) in hidxs_by_src
            _find_hexs!(read_hexs, write_hexs, chunk, halo_info, ccid, cid, hcid, idxs,
                        hasread, haswrite)
        end
    end
    return read_hexs, write_hexs
end

@inline function get_halo_read_buftype(s::AbstractStorage)
    return typeof(get_halo_read_fields(s))
end

@inline function get_halo_write_buftype(s::AbstractStorage)
    return typeof(get_halo_write_fields(s))
end

function _find_hexs!(rhexs::Vector, whexs::Vector, chunk::AbstractBodyChunk,
                     hi::MPIHaloInfo, ccid::Int, cid::Int, hcid::Int, idxs::Vector{Int},
                     hasread::Bool, haswrite::Bool)
    (ccid == cid || ccid == hcid) || return nothing
    hidxs = hi.point_ids[cid][idxs]
    localize!(hidxs, hi.localizers[hcid])
    if hasread
        bufs = get_halo_read_bufs(chunk.store, length(idxs))
        tags = get_tags(bufs, cid * mpi_nranks())
        n_fields = length(bufs)
        hex = HaloExchangeBuf(n_fields, hcid, cid, hidxs, idxs, bufs, tags)
        push!(rhexs, hex)
    end
    if haswrite
        bufs = get_halo_read_bufs(chunk.store, length(idxs))
        tags = get_tags(bufs, cid * mpi_nranks())
        n_fields = length(bufs)
        hex = HaloExchangeBuf(n_fields, cid, hcid, idxs, hidxs, bufs, tags)
        push!(whexs, hex)
    end
    return nothing
end

@inline function get_halo_read_bufs(s::AbstractStorage, n_halo_points::Int)
    return Tuple(get_buf(a, n_halo_points) for a in get_halo_read_fields(s))
end

@inline function get_halo_write_bufs(s::AbstractStorage, n_halo_points::Int)
    return Tuple(get_buf(a, n_halo_points) for a in get_halo_write_fields(s))
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

function exchange_read_fields!(dh::MPIDataHandler)
    for (hex_id, hex) in enumerate(dh.read_hexs_send)
        send_read_he!(dh, hex, hex_id)
    end
    for hex in dh.read_hexs_recv
        recv_read_he!(dh, hex)
    end
    MPI.Waitall(dh.read_reqs)
    return nothing
end

function send_read_he!(dh::MPIDataHandler, he::HaloExchangeBuf, hex_id::Int)
    src_fields = get_halo_read_fields(dh.chunk.store)
    for i in eachindex(src_fields, he.bufs)
        exchange_to_buf!(he.bufs[i], src_fields[i], he.src_idxs)
        req = MPI.Isend(he.bufs[i], mpi_comm(); dest=he.dest_chunk_id - 1, tag=he.tags[i])
        req_id = hex_id * he.n_fields + i - he.n_fields
        dh.read_reqs[req_id] = req
    end
    return nothing
end

function recv_read_he!(dh::MPIDataHandler, he::HaloExchangeBuf)
    dest_fields = get_halo_read_fields(dh.chunk.store)
    for i in eachindex(dest_fields, he.bufs)
        MPI.Recv!(he.bufs[i], mpi_comm(); source=he.src_chunk_id - 1, tag=he.tags[i])
        exchange_from_buf!(dest_fields[i], he.bufs[i], he.dest_idxs)
    end
    return nothing
end

function exchange_write_fields!(dh::MPIDataHandler)
    for (hex_id, hex) in enumerate(dh.write_hexs_send)
        send_write_he!(dh, hex, hex_id)
    end
    for hex in dh.write_hexs_recv
        recv_write_he!(dh, hex)
    end
    MPI.Waitall(dh.write_reqs)
    return nothing
end

function send_write_he!(dh::MPIDataHandler, he::HaloExchangeBuf, hex_id::Int)
    src_fields = get_halo_write_fields(dh.chunk.store)
    for i in eachindex(src_fields, he.bufs)
        exchange_to_buf!(he.bufs[i], src_fields[i], he.src_idxs)
        req = MPI.Isend(he.bufs[i], mpi_comm(); dest=he.dest_chunk_id - 1, tag=he.tags[i])
        req_id = hex_id * he.n_fields + i - he.n_fields
        dh.write_reqs[req_id] = req
    end
    return nothing
end

function recv_write_he!(dh::MPIDataHandler, he::HaloExchangeBuf)
    dest_fields = get_halo_write_fields(dh.chunk.store)
    for i in eachindex(dest_fields, he.bufs)
        MPI.Recv!(he.bufs[i], mpi_comm(); source=he.src_chunk_id - 1, tag=he.tags[i])
        exchange_from_buf_add!(dest_fields[i], he.bufs[i], he.dest_idxs)
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
