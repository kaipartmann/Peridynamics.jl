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
