mutable struct MultibodySetup{M<:AbstractMaterial,P<:AbstractPointParameters}
    bodies::Dict{Symbol,Body{M,P}}
    contacts::Vector{Contact}

    function MultibodySetup(bodies::Dict{Symbol,Body{M,P}}) where {M,P}
        contacts = Vector{Contact}()
        return new{M,P}(bodies, contacts)
    end
end

function MultibodySetup(::Dict{Symbol,Body})
    msg = "bodies have different material types!\n"
    msg *= "Only bodies with the same material types can be used for MultibodySetup!\n"
    throw(ArgumentError(msg))
    return nothing
end

MultibodySetup(body_pairs...) = MultibodySetup(Dict(body_pairs...))

@inline material_type(::MultibodySetup{M}) where {M} = M

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
