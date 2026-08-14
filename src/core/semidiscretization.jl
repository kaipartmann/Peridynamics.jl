"""
    semidiscretize(job; n_chunks=Threads.nthreads())

Create an `ODEProblem` from a [`Job`](@ref) for integration with OrdinaryDiffEq.jl.

The returned problem uses the initial conditions, material data, boundary conditions, and
time span stored in `job`. Its parameters contain both the job and the initialized data
handler, available as `ode.p.job` and `ode.p.data_handler`, respectively.

Currently, only jobs using [`VelocityVerlet`](@ref) are supported. The existing
[`submit`](@ref) interface is unaffected and continues to use Peridynamics.jl's built-in
time integration.

# Keywords
- `n_chunks::Integer=Threads.nthreads()`: Number of chunks used by the threaded data
    handler. This must be positive.

# Returns
A `DynamicalODEProblem`, represented by SciMLBase as an `ODEProblem`, with velocity and
displacement as its two state partitions.
"""
function semidiscretize(job::Job; n_chunks::Integer=nthreads())
    job.time_solver isa VelocityVerlet ||
        throw(ArgumentError("`semidiscretize` currently only supports jobs using " *
                            "`VelocityVerlet`"))
    mpi_run() && throw(ArgumentError("`semidiscretize` does not support MPI runs"))
    n_chunks > 0 || throw(ArgumentError("`n_chunks` must be positive"))

    dh = threads_data_handler(job.spatial_setup, job.time_solver, Int(n_chunks))
    init_time_solver!(job.time_solver, dh)
    initialize!(dh, job.time_solver)

    n_dofs = ode_n_dofs(dh)
    v0 = Vector{Float64}(undef, n_dofs)
    u0 = Vector{Float64}(undef, n_dofs)
    copy_initial_state!(v0, u0, dh)

    p = (; job, data_handler=dh)
    tspan = (0.0, job.time_solver.end_time)
    return SciMLBase.DynamicalODEProblem(ode_acceleration!, ode_velocity!, v0, u0, tspan,
                                        p)
end

function ode_n_dofs(dh::AbstractThreadsBodyDataHandler)
    return sum(get_n_loc_dof(chunk) for chunk in dh.chunks)
end

function ode_n_dofs(dh::AbstractThreadsMultibodyDataHandler)
    return sum(ode_n_dofs(body_dh) for body_dh in each_body_dh(dh))
end

function copy_initial_state!(v0, u0, dh::AbstractThreadsBodyDataHandler, offset=0)
    for chunk in dh.chunks
        n_dofs = get_n_loc_dof(chunk)
        copyto!(v0, offset + 1, chunk.storage.velocity, 1, n_dofs)
        copyto!(u0, offset + 1, chunk.storage.displacement, 1, n_dofs)
        offset += n_dofs
    end
    return offset
end

function copy_initial_state!(v0, u0, dh::AbstractThreadsMultibodyDataHandler, offset=0)
    for body_dh in each_body_dh(dh)
        offset = copy_initial_state!(v0, u0, body_dh, offset)
    end
    return offset
end

function update_ode_state!(v, u, dh::AbstractThreadsBodyDataHandler, t, offset=0)
    for chunk in dh.chunks
        (; system, storage) = chunk
        for dof in each_loc_dof(chunk)
            i = offset + dof
            @inbounds storage.velocity[dof] = v[i]
            @inbounds storage.velocity_half[dof] = v[i]
            @inbounds storage.displacement[dof] = u[i]
            @inbounds storage.position[dof] = system.position[dof] + u[i]
        end
        apply_boundary_conditions!(chunk, t)
        offset += get_n_loc_dof(chunk)
    end
    return offset
end

function update_ode_state!(v, u, dh::AbstractThreadsMultibodyDataHandler, t, offset=0)
    for body_dh in each_body_dh(dh)
        offset = update_ode_state!(v, u, body_dh, t, offset)
    end
    return offset
end

function ode_velocity!(du, v, u, p, t)
    dh = p.data_handler
    update_ode_state!(v, u, dh, t)
    copy_ode_velocity!(du, dh)
    return du
end

function copy_ode_velocity!(du, dh::AbstractThreadsBodyDataHandler, offset=0)
    for chunk in dh.chunks
        n_dofs = get_n_loc_dof(chunk)
        copyto!(du, offset + 1, chunk.storage.velocity_half, 1, n_dofs)
        offset += n_dofs
    end
    return offset
end

function copy_ode_velocity!(du, dh::AbstractThreadsMultibodyDataHandler, offset=0)
    for body_dh in each_body_dh(dh)
        offset = copy_ode_velocity!(du, body_dh, offset)
    end
    return offset
end

function ode_acceleration!(dv, v, u, p, t)
    dh = p.data_handler
    update_ode_state!(v, u, dh, t)
    calc_force_density!(dh, t, p.job.time_solver.Δt)
    prepare_ode_acceleration!(dh)
    copy_ode_acceleration!(dv, dh)
    return dv
end

prepare_ode_acceleration!(::AbstractThreadsBodyDataHandler) = nothing

function prepare_ode_acceleration!(dh::AbstractThreadsMultibodyDataHandler)
    update_caches!(dh)
    calc_contact_force_densities!(dh)
    for body_dh in each_body_dh(dh)
        for chunk_id in eachindex(body_dh.chunks)
            exchange_halo_to_loc!(body_dh, chunk_id)
        end
    end
    return nothing
end

function copy_ode_acceleration!(dv, dh::AbstractThreadsBodyDataHandler, offset=0)
    for chunk in dh.chunks
        update_acceleration!(chunk)
        n_dofs = get_n_loc_dof(chunk)
        copyto!(dv, offset + 1, chunk.storage.acceleration, 1, n_dofs)
        offset += n_dofs
    end
    return offset
end

function copy_ode_acceleration!(dv, dh::AbstractThreadsMultibodyDataHandler, offset=0)
    for body_dh in each_body_dh(dh)
        offset = copy_ode_acceleration!(dv, body_dh, offset)
    end
    return offset
end

function update_acceleration!(chunk::AbstractBodyChunk)
    _update_acceleration!(chunk, chunk.paramsetup)
    return nothing
end

function _update_acceleration!(chunk::AbstractBodyChunk,
                               paramhandler::AbstractParameterHandler)
    (; acceleration, b_int, b_ext) = chunk.storage
    for i in each_point_idx(chunk)
        params = get_params(paramhandler, i)
        for dim in each_dim(chunk)
            dof = get_dof(chunk.system, dim, i)
            _update_acc!(acceleration, b_int, b_ext, params.rho, dof)
        end
    end
    return nothing
end

function _update_acceleration!(chunk::AbstractBodyChunk, params::AbstractPointParameters)
    (; acceleration, b_int, b_ext) = chunk.storage
    for dof in each_loc_dof(chunk)
        _update_acc!(acceleration, b_int, b_ext, params.rho, dof)
    end
    return nothing
end
