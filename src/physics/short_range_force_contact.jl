const CONTACT_KWARGS = (:radius, :penalty_factor)

"""
    ShortRangeForceContact

$(internal_api_warning())

A type for contact simulations.

# Type Parameters

- `N`: Neighborhood search object used by `PointNeighbors.jl`.

# Fields

- `body_id_a::Symbol`: Index of a body of the contact.
- `body_id_b::Symbol`: Index of a body of the contact.
- `radius::Float64`: Search radius for contact.
- `penalty_factor::Float64`: Penalty factor for the contact simulation.
- `nhs::N`: Neighborhood search object used by `PointNeighbors.jl`.
"""
struct ShortRangeForceContact{N}
    body_id_a::Symbol
    body_id_b::Symbol
    radius::Float64
    penalty_factor::Float64
    nhs::N
end

function get_contact_params(p::Dict{Symbol,Any})
    if haskey(p, :radius)
        radius::Float64 = float(p[:radius])
    else
        throw(UndefKeywordError(:radius))
    end
    if radius ≤ 0
        throw(ArgumentError("`radius` should be larger than zero!\n"))
    end

    if haskey(p, :penalty_factor)
        penalty_factor::Float64 = float(p[:penalty_factor])
    else
        penalty_factor = 1e12
    end
    if penalty_factor ≤ 0
        throw(ArgumentError("`sc` should be larger than zero!\n"))
    end

    return radius, penalty_factor
end

"""
    contact!(multibody_setup, name_body_a, name_body_b; kwargs...)

Define a short range force contact between body `name_body_a` and `name_body_b` in the
[`MultibodySetup`](@ref) `multibody_setup`.

# Arguments

- `multibody_setup::MultibodySetup`: [`MultibodySetup`](@ref).
- `name_body_a::Symbol`: The name of a body in this multibody setup.
- `name_body_b::Symbol`: The name of a body in this multibody setup.

# Keywords

- `radius::Float64`: Contact search radius. If a the distance of a point in body
    `name_body_a` and a point in body `name_body_b` is lower than this radius, a contact
    force is calculated. This radius should be in the order of the point spacing of a
    point cloud.
- `penalty_factor::Float64`: Penalty factor for the short range force contact algorithm.
    (default: `1e12`)

# Throws
- Error if `multibody_setup` does not contain bodies with name `name_body_a` and
    `name_body_b`.
- Error if the keyword `radius` is not specified or `radius ≤ 0`.
- Error if `penalty_factor ≤ 0`.

# Examples
```julia-repl
julia> ms = MultibodySetup(:a => body_a, :b => body_b)
2000-point MultibodySetup:
  1000-point Body{BBMaterial{NoCorrection}} with name `b`
  1000-point Body{BBMaterial{NoCorrection}} with name `b`

julia> contact!(ms, :a, :b; radius=0.001)

julia> ms
2000-point MultibodySetup:
  1000-point Body{BBMaterial{NoCorrection}} with name `b`
  1000-point Body{BBMaterial{NoCorrection}} with name `b`
  2 short range force contact(s)
```
"""
function contact!(multibody_setup::AbstractMultibodySetup, name_body_a::Symbol,
                  name_body_b::Symbol; kwargs...)
    check_if_bodyname_is_defined(multibody_setup, name_body_a)
    check_if_bodyname_is_defined(multibody_setup, name_body_b)

    p = Dict{Symbol,Any}(kwargs)
    check_kwargs(p, CONTACT_KWARGS)
    radius, penalty_factor = get_contact_params(p)

    body_a = get_body(multibody_setup, name_body_a)
    body_b = get_body(multibody_setup, name_body_b)
    nhs_a = GridNeighborhoodSearch{3}(search_radius=radius, n_points=body_a.n_points,
                                      update_strategy=SemiParallelUpdate())
    nhs_b = GridNeighborhoodSearch{3}(search_radius=radius, n_points=body_b.n_points,
                                      update_strategy=SemiParallelUpdate())

    srfc_a = ShortRangeForceContact(name_body_a, name_body_b, radius, penalty_factor, nhs_b)
    srfc_b = ShortRangeForceContact(name_body_b, name_body_a, radius, penalty_factor, nhs_a)
    push!(multibody_setup.srf_contacts, srfc_a, srfc_b)

    return nothing
end

function init_srf_contacts_nhs!(dh::AbstractThreadsMultibodyDataHandler)
    for contact in dh.srf_contacts
        body_b_idx = dh.body_idxs[contact.body_id_b]
        posc_b = dh.position_caches[body_b_idx]
        initialize_grid!(contact.nhs, posc_b)
    end
    return nothing
end

function calc_short_range_force_contacts!(dh::AbstractThreadsMultibodyDataHandler)
    for contact in dh.srf_contacts
        body_a_idx = dh.body_idxs[contact.body_id_a]
        body_b_idx = dh.body_idxs[contact.body_id_b]

        # update neighborhoodsearches
        posc_a = dh.position_caches[body_a_idx]
        posc_b = dh.position_caches[body_b_idx]
        update_grid!(contact.nhs, posc_b; parallelization_backend=ThreadsStaticBackend())

        # calc contact
        body_dh_a = get_body_dh(dh, body_a_idx)
        volc_b = dh.volume_caches[body_b_idx]
        @threads :static for chunk_a in body_dh_a.chunks
            calc_contact_force_density!(chunk_a, contact, posc_a, posc_b, volc_b)
        end
    end
    return nothing
end

function calc_contact_force_density!(chunk_a, contact, pos_cache_a, pos_cache_b,
                                     vol_cache_b)
    storage_a = chunk_a.storage
    nhs = contact.nhs
    r = contact.radius
    C0 = 9 * contact.penalty_factor / π
    for (li, ii) in each_point_idx_pair(chunk_a)
        params_a = get_params(chunk_a, li)
        C = C0 / params_a.δ^4
        foreach_neighbor(pos_cache_a, pos_cache_b, nhs, ii) do i, j, Δxij, L
            temp = C * (r - L) / L * vol_cache_b[j]
            storage_a.b_int[1, li] += temp * Δxij[1]
            storage_a.b_int[2, li] += temp * Δxij[2]
            storage_a.b_int[3, li] += temp * Δxij[3]
        end
    end
    return nothing
end
