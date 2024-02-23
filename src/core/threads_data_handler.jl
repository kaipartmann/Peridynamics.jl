struct ThreadsDataHandler{C<:AbstractBodyChunk} <: AbstractDataHandler
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
    body_chunks = chop_body_threads(body, time_solver, point_decomp, v)
    _halo_exchanges = find_halo_exchanges(body_chunks)
    halo_exchanges = [Vector{HaloExchange}() for _ in eachindex(body_chunks)]
    @threads :static for chunk_id in eachindex(body_chunks)
        halo_exchanges[chunk_id] = filter(x -> x.dest_chunk_id == chunk_id, _halo_exchanges)
    end
    return ThreadsDataHandler(body_chunks, halo_exchanges)
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

function _export_results(dh::ThreadsDataHandler, options::ExportOptions, n::Int, t::Float64)
    filename = joinpath(options.vtk, @sprintf("timestep_%05d", n))
    n_chunks = length(dh.chunks)
    @threads :static for chunk_id in eachindex(dh.chunks)
        chunk = dh.chunks[chunk_id]
        n_local_points = length(chunk.ch.loc_points)
        position = @views chunk.store.position[:, 1:n_local_points]
        cells = get_cells(n_local_points)
        pvtk_grid(filename, position, cells; part=chunk_id, nparts=n_chunks) do vtk
            for fld in options.fields
                vtk[string(fld), VTKPointData()] = getfield(chunk.store, fld)
            end
            vtk["time", VTKFieldData()] = t
        end
    end
    return nothing
end

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
