const CONTACT_KWARGS = (:radius, :sc)

struct Contact
    body_id_a::Symbol
    body_id_b::Symbol
    radius::Float64
    sc::Float64
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

    if haskey(p, :sc)
        sc = float(p[:sc])
    else
        sc = 1e12
    end
    if sc ≤ 0
        throw(ArgumentError("`sc` should be larger than zero!\n"))
    end

    return radius, sc
end
