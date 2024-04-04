"""
    MultibodySetup{M<:AbstractMaterial,P<:AbstractPointParameters}

setup for a peridynamic simulation with multiple bodies

# Fields

- `bodies::Dict{Symbol,Body{M,P}}`:
- `contacts::Vector{Contact}`:

---

Constructors:

    MultibodySetup(bodies::Dict{Symbol,Body{M,P}}) where {M,P}

# Arguments

- `bodies::Dict{Symbol,Body{M,P}}`:

# Throws

- error if less than 2 bodies are defined
- error if defined bodies have different material types

# Example

TODO
"""
mutable struct MultibodySetup{M<:AbstractMaterial,
                              P<:AbstractPointParameters} <: AbstractMultibodySetup
    bodies::Dict{Symbol,Body{M,P}}
    contacts::Vector{Contact}

    function MultibodySetup(bodies::Dict{Symbol,Body{M,P}}) where {M,P}
        if length(bodies) < 2
            msg = "not enough bodies given!\n"
            msg *= "Please specify 2 or more!\n"
            throw(ArgumentError(msg))
        end
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

"""
    contact!(ms::MultibodySetup, body_a::Symbol, body_b::Symbol; kwargs...)

defines contact between multiple bodies

# Arguments

- `ms::MultibodySetup`: the multibody setup defined to simulate the contact
- `body_a::Symbol`: first body in contact
- `body_b::Symbol`: second body in contact

# Keywords

- `:radius::Float64`:
- `:sc::Float64`:

# Throws

- error if a called body is not defined in the multibody setup
- error if keyword is not allowed

TODO kwargs
"""
function contact!(ms::MultibodySetup, body_a::Symbol, body_b::Symbol; kwargs...)
    check_if_bodyname_is_defined(ms, body_a)
    check_if_bodyname_is_defined(ms, body_b)

    p = Dict{Symbol,Any}(kwargs)
    check_kwargs(p, CONTACT_KWARGS)
    radius, sc = get_contact_params(p)

    push!(ms.contacts, Contact(body_a, body_b, radius, sc))
    return nothing
end

function pre_submission_check(ms::MultibodySetup)
    #TODO: check if everything is defined for job submission!
    return nothing
end

@inline function storage_type(ms::MultibodySetup, ts::AbstractTimeSolver)
    return storage_type(first(ms.bodies).mat, ts)
end
