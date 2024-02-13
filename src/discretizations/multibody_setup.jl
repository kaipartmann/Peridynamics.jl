mutable struct MultibodySetup{M}
    bodies::Dict{Symbol,Body{M}}
    contacts::Vector{Contact}

    function MultibodySetup(body_pairs::Vararg{Pair{Symbol,Body{M,P}},N}) where {M,P,N}
        bodies = Dict{Symbol,Body}()
        for body_pair in body_pairs
            bodies[body_pair.first] = body_pair.second
        end
        check_if_same_material(bodies)
        contacts = Vector{Contact}()
        return new{M}(bodies, contacts)
    end
end

@inline material_type(::MultibodySetup{M}) where {M} = M

function check_if_same_material(bodies)
    mat_types = material_type.(values(bodies))
    if !allequal(mat_types)
        msg = "bodies have different material types!\n"
        msg *= "Only bodies with the same material types can be used for MultibodySetup!\n"
        throw(ArgumentError(msg))
    end
    return nothing
end

function check_if_bodyname_is_defined(ms::MultibodySetup, name::Symbol)
    if !haskey(ms.bodies, name)
        throw(ArgumentError("there is no body with name $(name)!"))
    end
    return nothing
end

function contact!(ms::MultibodySetup, body_a::Symbol, body_b::Symbol; kwargs...)
    check_if_bodyname_is_defined(ms, body_a)
    check_if_bodyname_is_defined(ms, body_b)

    p = Dict{Symbol,Any}(kwargs)
    check_kwargs(p, CONTACT_KWARGS)
    radius, sc = get_contact_params(p)

    push!(ms.contacts, Contact(body_a, body_b, radius, sc))
    return nothing
end
