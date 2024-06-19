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

    body_b = get_body(ms, name_body_b)
    nhs = GridNeighborhoodSearch{3}(radius, body_b.n_points)

    srfc = ShortRangeForceContact(name_body_a, name_body_b, radius, penalty_factor, nhs)
    push!(ms.srf_contacts, srfc)

    return nothing
end
