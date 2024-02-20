struct ThreadsDataHandler{C<:AbstractBodyChunk}
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

function _export_results(dh::ThreadsDataHandler, options::ExportOptions, n::Int, t::Float64)
    filename = joinpath(options.vtk, @sprintf("timestep_%05d", n))
    n_chunks = length(dh.chunks)
    @threads :static for chunk_id in eachindex(dh.chunks)
        chunk = dh.chunks[chunk_id]
        position = @views chunk.dscr.position[:, chunk.ch.loc_points]
        pvtk_grid(filename, position, sp.cells; part=chunk_id, nparts=n_chunks) do vtk
            for fld in options.fields
                vtk[string(fld)] = getfield(chunk.store, fld)
            end
            vtk["time", VTKFieldData()] = t
        end
    end
    return nothing
end
