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
        halo_exchanges[chunk_id] = filter(x -> x.to_chunk_id == chunk_id, _halo_exchanges)
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

function halo_exchange!(dh::ThreadsDataHandler, chunk_id)
    for he in dh.halo_exchanges[chunk_id]
        _halo_exchange!(dh, he)
    end
    return nothing
end

function _halo_exchange!(dh::ThreadsDataHandler, he::HaloExchange)
    from_chunk = dh.chunks[he.from_chunk_id]
    to_chunk = dh.chunks[he.to_chunk_id]
    field = he.field
    from_loc_idxs = he.from_loc_idxs
    to_loc_idxs = he.to_loc_idxs
    for i in eachindex(from_loc_idxs, to_loc_idxs)
        from_idx = from_loc_idxs[i]
        to_idx = to_loc_idxs[i]
        from_store_entry = getfield(from_chunk.store, field)
        to_store_entry = getfield(to_chunk.store, field)
        to_store_entry[1, to_idx] = from_store_entry[1, from_idx]
        to_store_entry[2, to_idx] = from_store_entry[2, from_idx]
        to_store_entry[3, to_idx] = from_store_entry[3, from_idx]
    end
    return nothing
end
