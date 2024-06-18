"""
    MultibodySetup(bodies)

Setup for a peridynamic simulation with multiple bodies

# Arguments

- `bodies::Dict{Symbol,Body{M,P}}`:

# Throws

- Error if less than 2 bodies are defined
- Error if defined bodies have different material types

# Example

---

!!! warning "Internal use only"
    Please note that the fields are intended for internal use only. They are *not* part of
    the public API of Peridynamics.jl, and thus can be altered (or removed) at any time
    without it being considered a breaking change.

```julia
MultibodySetup{Material,PointParameters}
```

# Type Parameters

- `Material <: AbstractMaterial`: Type of the specified material model
- `PointParameters <: AbstractPointParameters`: Type of the point parameters

# Fields

- `bodies::Dict{Symbol,Body{M,P}}`:
- `srf_contacts::Vector{ShortRangeForceContact}`:

TODO
"""
mutable struct MultibodySetup{M<:AbstractMaterial,
                              P<:AbstractPointParameters} <: AbstractMultibodySetup{M}
    bodies::Dict{Symbol,Body{M,P}}
    srf_contacts::Vector{ShortRangeForceContact}

    function MultibodySetup(bodies::Dict{Symbol,Body{M,P}}) where {M,P}
        if length(bodies) < 2
            msg = "not enough bodies given!\n"
            msg *= "Please specify 2 or more!\n"
            throw(ArgumentError(msg))
        end
        srf_contacts = Vector{ShortRangeForceContact}()
        return new{M,P}(bodies, srf_contacts)
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

function check_if_bodyname_is_defined(ms::AbstractMultibodySetup, name::Symbol)
    if !haskey(ms.bodies, name)
        throw(ArgumentError("there is no body with name $(name)!"))
    end
    return nothing
end

@inline get_body(ms::AbstractMultibodySetup, name::Symbol) = ms.bodies[name]

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

function pre_submission_check(ms::AbstractMultibodySetup)
    #TODO: check if everything is defined for job submission!
    return nothing
end

@inline function storage_type(ms::AbstractMultibodySetup, ts::AbstractTimeSolver)
    return storage_type(first(values(ms.bodies)).mat, ts)
end
