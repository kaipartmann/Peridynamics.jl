"""
    MultibodySetup(body_pairs...)

Setup for a peridynamic simulation with multiple bodies.

# Arguments

- `body_pairs::Pair{Symbol,<:AbstractBody}`: Pairs of `:body_name => body_object`.
    The name of the body has to be specified as a Symbol.

# Throws

- Errors if less than 2 bodies are defined

# Examples

```julia-repl
julia> sphere = Body(BBMaterial(), pos_sphere, vol_sphere)
280-point Body{BBMaterial{NoCorrection}}:
  1 point set(s):
    280-point set `all_points`

julia> plate = Body(BBMaterial(), pos_plate, vol_plate)
25600-point Body{BBMaterial{NoCorrection}}:
  1 point set(s):
    25600-point set `all_points`

julia> ms = MultibodySetup(:sphere => sphere, :plate => plate)
25880-point MultibodySetup:
  280-point Body{BBMaterial{NoCorrection}} with name `sphere`
  25600-point Body{BBMaterial{NoCorrection}} with name `plate`
```

---

!!! warning "Internal use only"
    Please note that the fields are intended for internal use only. They are *not* part of
    the public API of Peridynamics.jl, and thus can be altered (or removed) at any time
    without it being considered a breaking change.

```julia
MultibodySetup{Bodies}
```

# Type Parameters

- `Bodies <: Tuple`: All types of the different bodies in the multibody setup.

# Fields

- `bodies::Bodies`: A Tuple containing all the bodies.
- `body_names::Vector{Symbol}`: All body names.
- `body_idxs::Dict{Symbol,Int}`: A Dict to get the body index with the body name.
- `srf_contacts::Vector{ShortRangeForceContact}`: All short range force contacts.
"""
struct MultibodySetup{B} <: AbstractMultibodySetup
    bodies::B
    body_names::Vector{Symbol}
    body_idxs::Dict{Symbol,Int}
    srf_contacts::Vector{ShortRangeForceContact}
end

function MultibodySetup(body_pairs::Vararg{Pair{Symbol,<:AbstractBody},N}) where {N}
    if N < 2
        msg = "not enough bodies specified!\n"
        msg *= "Specify at least 2 bodies for a `MultibodySetup`!\n"
        throw(ArgumentError(msg))
    end
    bodies = Tuple(body.second for body in body_pairs)
    body_names = [body_pair.first for body_pair in body_pairs]
    for (body_id, body) in enumerate(bodies)
        change_name!(body, body_names[body_id])
    end
    body_idxs = Dict{Symbol,Int}()
    for (i, name) in enumerate(body_names)
        body_idxs[name] = i
    end
    srf_contacts = Vector{ShortRangeForceContact}()
    return MultibodySetup(bodies, body_names, body_idxs, srf_contacts)
end

function Base.show(io::IO, ms::MultibodySetup)
    print(io, n_points(ms), "-point MultibodySetup")
    return nothing
end

function Base.show(io::IO, ::MIME"text/plain", ms::MultibodySetup)
    if get(io, :compact, false)
        show(io, ms)
        return nothing
    end
    print(io, n_points(ms), "-point MultibodySetup:")
    for body in each_body(ms)
        print(io, "\n  ")
        show(io, body)
    end
    n_contact = length(ms.srf_contacts)
    if n_contact > 0
        print(io, "\n  ", n_contact, " short range force contact(s)")
    end
    return nothing
end

function check_if_bodyname_is_defined(ms::AbstractMultibodySetup, name::Symbol)
    if !haskey(ms.body_idxs, name)
        throw(ArgumentError("there is no body with name $(name)!"))
    end
    return nothing
end

@inline get_body(ms::AbstractMultibodySetup, name::Symbol) = ms.bodies[ms.body_idxs[name]]
@inline get_body(ms::AbstractMultibodySetup, idx::Int) = ms.bodies[idx]
@inline get_body_name(ms::AbstractMultibodySetup, idx::Int) = string(ms.body_names[idx])

@inline each_body(ms::AbstractMultibodySetup) = ms.bodies
@inline each_body_idx(ms::AbstractMultibodySetup) = eachindex(ms.bodies)
@inline each_body_name(ms::AbstractMultibodySetup) = ms.body_names

function pre_submission_check(ms::AbstractMultibodySetup)
    if mpi_run()
        @mpiroot begin
            @error "Multibody simulations with MPI are not yet implemented!\n"
            MPI.Abort(mpi_comm(), 1)
        end
        MPI.Barrier(mpi_comm())
    end

    # check all bodies separately
    n_all_bcs, n_all_ics = 0, 0
    for body in each_body(ms)
        pre_submission_check(body; body_in_multibody_setup=true)
        n_all_bcs += n_bcs(body)
        n_all_ics += n_ics(body)
    end

    # there has to be any condition in one of the bodies
    if n_all_bcs + n_all_ics == 0
        msg = "no initial or boundary condition specified!\n"
        msg *= "At least one initial or boundary condition is needed in one of "
        msg *= "the bodies!\nOtherwise nothing will happen during the simulation!\n"
        error(msg)
    end

    return nothing
end

function log_spatial_setup(options::AbstractJobOptions, ms::MultibodySetup)
    for body in each_body(ms)
        log_spatial_setup(options, body)
    end
    return nothing
end

@inline n_points(ms::AbstractMultibodySetup) = sum(x -> n_points(x), each_body(ms))
