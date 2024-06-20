const CONTACT_KWARGS = (:radius, :penalty_factor)

struct ShortRangeForceContact{N}
    body_id_a::Symbol
    body_id_b::Symbol
    radius::Float64
    penalty_factor::Float64
    nhs::N
end

function get_contact_params(p::Dict{Symbol,Any})
    local radius::Float64
    local sc::Float64

    if haskey(p, :radius)
        radius = float(p[:radius])
    else
        throw(UndefKeywordError(:radius))
    end
    if radius ≤ 0
        throw(ArgumentError("`radius` should be larger than zero!\n"))
    end

    if haskey(p, :penalty_factor)
        penalty_factor = float(p[:penalty_factor])
    else
        penalty_factor = 1e12
    end
    if penalty_factor ≤ 0
        throw(ArgumentError("`sc` should be larger than zero!\n"))
    end

    return radius, penalty_factor
end

"""
    contact!(ms, body_a, body_b; kwargs...)

Defines contact between multiple bodies

# Arguments

- `ms::AbstractMultibodySetup`: The multibody setup defined to simulate the contact
- `body_a::Symbol`: First body in contact
- `body_b::Symbol`: Second body in contact

# Keywords

- `radius::Float64`:
- `sc::Float64`:

# Throws

- Error if a called body is not defined in the multibody setup
- Error if keyword is not allowed

TODO kwargs
"""
function contact!(ms::AbstractMultibodySetup, name_body_a::Symbol, name_body_b::Symbol;
                  kwargs...)
    check_if_bodyname_is_defined(ms, name_body_a)
    check_if_bodyname_is_defined(ms, name_body_b)

    p = Dict{Symbol,Any}(kwargs)
    check_kwargs(p, CONTACT_KWARGS)
    radius, penalty_factor = get_contact_params(p)

    body_a = get_body(ms, name_body_a)
    body_b = get_body(ms, name_body_b)
    nhs_a = GridNeighborhoodSearch{3}(radius, body_a.n_points)
    nhs_b = GridNeighborhoodSearch{3}(radius, body_b.n_points)

    srfc_a = ShortRangeForceContact(name_body_a, name_body_b, radius, penalty_factor, nhs_b)
    srfc_b = ShortRangeForceContact(name_body_b, name_body_a, radius, penalty_factor, nhs_a)
    push!(ms.srf_contacts, srfc_a, srfc_b)

    return nothing
end

function calc_short_range_force_contacts!(dh)
    for contact in dh.srf_contacts
        body_a_idx = dh.body_idxs[contact.body_id_a]
        body_b_idx = dh.body_idxs[contact.body_id_b]

        # update neighborhoodsearches
        posc_a = dh.position_caches[body_a_idx]
        posc_b = dh.position_caches[body_b_idx]
        PointNeighbors.initialize!(contact.nhs, posc_a, posc_b)

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
            @assert i == ii
            temp = C * (r - L) / L
            storage_a.b_int[1, li] += temp * Δxij[1] * vol_cache_b[j]
            storage_a.b_int[2, li] += temp * Δxij[2] * vol_cache_b[j]
            storage_a.b_int[3, li] += temp * Δxij[3] * vol_cache_b[j]
        end
    end
    return nothing
end
