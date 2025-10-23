"""
    ThreadsMultibodyDataHandler

$(internal_api_warning())

A type for handling all data of multiple bodies in multithreading simulations.

# Type Parameters

- `BDH`: Body data handler type.
- `PC`: Position cache type.
- `VC`: Volume cache type.

# Fields

- `n_bodies::Int`: Number of bodies in the simulation.
- `body_dhs::BDH`: Tuple containing all body data handlers.
- `body_names::Vector{Symbol}`: Names of the bodies.
- `body_idxs::Dict{Symbol,Int}`: Names of bodies assigned to their indices.
- `srf_contacts::Vector{ShortRangeForceContact}`: All short range force contacts of this
    simulation.
- `position_caches::PC`: Positions of all points of all bodies (should be of type
    `Vector{Matrix{Float64}}`).
- `volume_caches::VC`: Volumes of all points of all bodies (should be of type
    `Vector{Vector{Float64}}`).
"""
struct ThreadsMultibodyDataHandler{BDH,PC,VC} <: AbstractThreadsMultibodyDataHandler
    n_bodies::Int
    body_dhs::BDH
    body_names::Vector{Symbol}
    body_idxs::Dict{Symbol,Int}
    srf_contacts::Vector{ShortRangeForceContact}
    position_caches::PC #Vector{Matrix{Float64}}
    volume_caches::VC #Vector{Vector{Float64}}
end

function threads_data_handler(ms::AbstractMultibodySetup, solver::AbstractTimeSolver,
                              n_chunks::Int)
    body_dhs = Tuple(threads_data_handler(body, solver, n_chunks) for body in each_body(ms))
    position_caches = [body.position for body in each_body(ms)]
    volume_caches = [body.volume for body in each_body(ms)]
    n_bodies = length(body_dhs)
    dh = ThreadsMultibodyDataHandler(n_bodies, body_dhs, ms.body_names,
                                     ms.body_idxs, ms.srf_contacts, position_caches,
                                     volume_caches)
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

function calc_force_density!(dh::ThreadsMultibodyDataHandler, t, Δt)
    for body_idx in each_body_idx(dh)
        body_dh = get_body_dh(dh, body_idx)
        @threads :static for chunk_id in eachindex(body_dh.chunks)
            exchange_loc_to_halo!(body_dh, chunk_id)
            calc_force_density!(body_dh.chunks[chunk_id], t, Δt)
        end
    end
    return nothing
end

function update_caches!(dh::ThreadsMultibodyDataHandler)
    for body_idx in each_body_idx(dh)
        body_dh = get_body_dh(dh, body_idx)
        position_cache = dh.position_caches[body_idx]
        @threads :static for chunk in body_dh.chunks
            @inbounds for (li, i) in enumerate(get_loc_points(chunk.system))
                position_cache[1, i] = chunk.storage.position[1, li]
                position_cache[2, i] = chunk.storage.position[2, li]
                position_cache[3, i] = chunk.storage.position[3, li]
            end
        end
    end
    return nothing
end

function initialize!(dh::AbstractThreadsMultibodyDataHandler, solver::AbstractTimeSolver)
    init_contact_nhs!(dh)
    calc_force_density!(dh, 0.0, solver.Δt)
    return nothing
end

@inline function init_contact_nhs!(dh::AbstractThreadsMultibodyDataHandler)
    init_srf_contacts_nhs!(dh)
    return nothing
end

function calc_contact_force_densities!(dh)
    calc_short_range_force_contacts!(dh)
    return nothing
end

function export_reference_results(dh::ThreadsMultibodyDataHandler,
                                  options::AbstractJobOptions)
    options.export_allowed || return nothing
    for body_idx in each_body_idx(dh)
        body_dh = get_body_dh(dh, body_idx)
        export_reference_results(body_dh, options)
    end
    return nothing
end

function log_data_handler(options::AbstractJobOptions,
                          dh::AbstractThreadsMultibodyDataHandler)
    for body_dh in each_body_dh(dh)
        log_system(options, body_dh)
    end
    return nothing
end
