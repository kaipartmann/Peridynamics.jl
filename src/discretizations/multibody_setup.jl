mutable struct MultibodySetup
    bodies::Dict{Symbol,Body}
    contacts::Vector{Contact}
    function MultibodySetup(body_pairs::B...) where {B<:Pair{Symbol,Body}}
        bodies = Dict{Symbol,Body}()
        for body_pair in body_pairs
            bodies[body_pair.first] = body_pair.second
        end
        contacts = Vector{Contact}()
        return new(bodies, contacts)
    end
end

function check_if_bodyname_is_defined(ms::MultibodySetup, name::Symbol)
    if !haskey(ms.bodies, name)
        throw(ArgumentError("there is no body with name $(name)!"))
    end
    return nothing
end

function contact!(ms::MultibodySetup, body_a::Symbol, body_b::Symbol, contact_radius::Real,
                  contact_constant::Real=1e12)
    check_if_bodyname_is_defined(ms, body_a)
    check_if_bodyname_is_defined(ms, body_b)
    if contact_radius ≤ 0
        throw(ArgumentError("contact radius should be larger than zero!\n"))
    end
    if contact_constant ≤ 0
        throw(ArgumentError("contact spring constant should be larger than zero!\n"))
    end
    push!(ms.contacts, Contact(body_a, body_b, contact_radius, contact_constant))
end
