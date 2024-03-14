struct HaloExchangeBuf{T,K}
    src_chunk_id::Int
    dest_chunk_id::Int
    src_idxs::Vector{Int}
    dest_idxs::Vector{Int}
    bufs::T
    tags::K
    function HaloExchangeBuf(src_chunk_id::Int, dest_chunk_id::Int,
                             src_idxs::AbstractVector{Int}, dest_idxs::AbstractVector{Int},
                             bufs::T, tags::K) where {T,K}
        return new{T,K}(src_chunk_id, dest_chunk_id, src_idxs, dest_idxs, bufs, tags)
    end
end

struct MPIHaloInfo
    point_ids::Dict{Int,Vector{Int}}
    halos_points::Dict{Int,Vector{Int}}
    localizers::Dict{Int,Dict{Int,Int}}
    hidxs_by_src::Dict{Int,Dict{Int,Vector{Int}}}
end

struct MPIDataHandler{C<:AbstractBodyChunk,RB,WB} <: AbstractDataHandler
    chunk::C
    read_halo_exs::Vector{HaloExchangeBuf{RB}}
    write_halo_exs::Vector{HaloExchangeBuf{WB}}
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
    halo_infos = get_halo_infos(chunk, point_decomp, body.n_points)
    read_halo_exs, write_halo_exs = find_halo_exchange_bufs(chunk, halo_infos)
    return MPIDataHandler(chunk, read_halo_exs, write_halo_exs)
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
    chunk = BodyChunk(body, ts, pd, mpi_chunk_id())
    apply_precracks!(chunk, body)
    apply_initial_conditions!(chunk, body)
    return chunk
end

function _chop_body_mpi(body::Body{M,P}, ts::T, pd::PointDecomposition,
                        ::Val{N}) where {M,P,T,N}
    chunk = MultiParamBodyChunk(body, ts, pd, mpi_chunk_id())
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

function get_halos_by_chunk_id(halo_chunk_ids::Vector{Int})
    halos_by_chunk_id = Dict{Int,Vector{Int}}()
    unique_chunk_ids = unique(halo_chunk_ids)
    for chunk_id in unique_chunk_ids
        ids = findall(x -> x == chunk_id, halo_chunk_ids)
        halos_by_chunk_id[chunk_id] = ids
    end
    return halos_by_chunk_id
end

function find_halo_exchange_bufs(chunk::AbstractBodyChunk, halo_infos::MPIHaloInfo)
    has_read = has_read_halos(chunk)
    has_write = has_write_halos(chunk)
    read_exs, write_exs = find_exs(chunk, halo_infos, has_read, has_write)
    return read_exs, write_exs
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

function find_exs(chunk::AbstractBodyChunk, halo_info::MPIHaloInfo, hasread::Bool,
                  haswrite::Bool)
    RB = get_halo_read_buftype(chunk.store)
    WB = get_halo_write_buftype(chunk.store)
    readexs = Vector{HaloExchangeBuf{RB}}()
    writeexs = Vector{HaloExchangeBuf{WB}}()
    # current_chunk_id = mpi_chunk_id()
    for (chunk_id, hidxs_by_src) in halo_info.hidxs_by_src
        for (halo_chunk_id, idxs) in hidxs_by_src
            halo_idxs = halo_info.point_ids[chunk_id][idxs]
            localize!(halo_idxs, halo_info.localizers[halo_chunk_id])
            if hasread
                src_id = halo_chunk_id
                dest_id = chunk_id
                src_idxs = halo_idxs
                dest_idxs = idxs
                bufs = get_halo_read_bufs(chunk.store, length(idxs))
                tags = get_tags(bufs, chunk_id * mpi_nranks())
                heb = HaloExchangeBuf(src_id, dest_id, src_idxs, dest_idxs, bufs, tags)
                push!(readexs, heb)
            end
            if haswrite
                src_id = chunk_id
                dest_id = halo_chunk_id
                src_idxs = idxs
                dest_idxs = halo_idxs
                bufs = get_halo_read_bufs(chunk.store, length(idxs))
                tags = get_tags(bufs, chunk_id * mpi_nranks())
                heb = HaloExchangeBuf(src_id, dest_id, src_idxs, dest_idxs, bufs, tags)
                push!(writeexs, heb)
            end
        end
    end
    return readexs, writeexs
end

@inline function get_halo_read_buftype(s::AbstractStorage)
    return typeof(get_halo_read_fields(s))
end

@inline function get_halo_write_buftype(s::AbstractStorage)
    return typeof(get_halo_write_fields(s))
end

@inline function get_halo_read_bufs(s::AbstractStorage, n_halo_points::Int)
    return Tuple(get_buf(a, n_halo_points) for a in get_halo_read_fields(s))
end

@inline function get_halo_write_bufs(s::AbstractStorage, n_halo_points::Int)
    return Tuple(get_buf(a, n_halo_points) for a in get_halo_write_fields(s))
end

@inline get_buf(a::Matrix{T}, n::Int) where {T} = zeros(T, size(a,1), n)
@inline get_buf(::Vector{T}, n::Int) where {T} = zeros(T, n)

@inline function get_tags(bufs, i)
    return Tuple(i + j for j in eachindex(bufs))
end

function calc_stable_timestep(dh::MPIDataHandler, safety_factor::Float64)
    _Δt = calc_timestep(dh.chunk)
    Δt = MPI.Allreduce(_Δt, MPI.MIN, mpi_comm())
    return Δt * safety_factor
end

function exchange_read_fields!(dh::MPIDataHandler)
    for he in dh.read_halo_exs
        if mpi_chunk_id() == he.src_chunk_id
            src_fields = get_halo_read_fields(dh.chunk.store)
            for i in eachindex(src_fields, he.bufs)
                exchange_to_buf!(he.bufs[i], src_fields[i], he.src_idxs)
                MPI.Isend(he.bufs[i], mpi_comm(); dest=he.dest_chunk_id-1, tag=he.tags[i])
            end
        end
    end
    for he in dh.read_halo_exs
        if mpi_chunk_id() == he.dest_chunk_id
            dest_fields = get_halo_read_fields(dh.chunk.store)
            for i in eachindex(dest_fields, he.bufs)
                MPI.Recv!(he.bufs[i], mpi_comm(); source=he.src_chunk_id-1, tag=he.tags[i])
                exchange_from_buf!(dest_fields[i], he.bufs[i], he.dest_idxs)
            end
        end
    end
    return nothing
end

function exchange_write_fields!(dh::MPIDataHandler)
    for he in dh.write_halo_exs
        if mpi_chunk_id() == he.src_chunk_id
            src_fields = get_halo_write_fields(dh.chunk.store)
            for i in eachindex(src_fields, he.bufs)
                exchange_to_buf!(he.bufs[i], src_fields[i], he.src_idxs)
                MPI.Isend(he.bufs[i], mpi_comm(); dest=he.dest_chunk_id-1, tag=he.tags[i])
            end
        end
    end
    for he in dh.write_halo_exs
        if mpi_chunk_id() == he.dest_chunk_id
            dest_fields = get_halo_write_fields(dh.chunk.store)
            for i in eachindex(dest_fields, he.bufs)
                MPI.Recv!(he.bufs[i], mpi_comm(); source=he.src_chunk_id-1, tag=he.tags[i])
                exchange_from_buf_add!(dest_fields[i], he.bufs[i], he.dest_idxs)
            end
        end
    end
    return nothing
end

function export_results(dh::MPIDataHandler, options::ExportOptions, n::Int, t::Float64)
    options.exportflag || return nothing
    if mod(n, options.freq) == 0
        _export_results(dh.chunk, mpi_chunk_id(), mpi_nranks(), options, n, t)
    end
    return nothing
end

function export_reference_results(dh::MPIDataHandler, options::ExportOptions)
    options.exportflag || return nothing
    _export_results(dh.chunk, mpi_chunk_id(), mpi_nranks(), options, 0, 0.0)
    return nothing
end
