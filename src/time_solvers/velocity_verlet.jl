"""
    VelocityVerlet(; kwargs...)

Time integration solver for the Velocity Verlet algorithm. Specify either the number of
steps or the time the simulation should cover.

# Keywords
- `time::Real`: The total time the simulation will cover. If this keyword is specified, the
    keyword `steps` is no longer allowed. (optional)
- `steps::Int`: Number of calculated time steps. If this keyword is specified, the keyword
    `time` is no longer allowed. (optional)
- `stepsize::Real`: Manually specify the size of the time step. (optional)
- `safety_factor::Real`: Safety factor for step size to ensure stability. (default: `0.7`)

!!! warning "Specification of the time step"
    Keep in mind that manually specifying the critical time step is dangerous! If the
    specified time step is too high and the CFL condition no longer holds, the simulation
    will give wrong results and maybe crash!

# Throws
- Errors if both `time` and `steps` are specified as keywords.
- Errors if neither `time` nor `steps` are specified as keywords.
- Errors if `safety_factor < 0` or `safety_factor > 1`.

# Example

```julia-repl
julia> VelocityVerlet(steps=2000)
VelocityVerlet:
  n_steps        2000
  safety_factor  0.7

julia> VelocityVerlet(time=0.001)
VelocityVerlet:
  end_time       0.001
  safety_factor  0.7

julia> VelocityVerlet(steps=2000, stepsize=0.0001)
┌ Warning: stepsize specified! Please be sure that the CFD-condition holds!
└ @ Peridynamics ~/Code/Peridynamics.jl/src/time_solvers/velocity_verlet.jl:66
VelocityVerlet:
  n_steps        2000
  Δt             0.0001
  safety_factor  0.7
```
"""
mutable struct VelocityVerlet <: AbstractTimeSolver
    end_time::Float64
    n_steps::Int
    Δt::Float64
    safety_factor::Float64

    function VelocityVerlet(; time::Real=-1, steps::Int=-1, stepsize::Real=-1,
                            safety_factor::Real=0.7)
        if time > 0 && steps > 0
            msg = "specify either time or number of steps, not both!"
            throw(ArgumentError(msg))
        elseif time < 0 && steps < 0
            msg = "specify either time or number of steps!"
            throw(ArgumentError(msg))
        end
        if !(0 < safety_factor < 1)
            msg = "wrong safety factor specified! condition: 0 < safety_factor < 1"
            throw(ArgumentError(msg))
        end
        if stepsize > 0
            @warn "stepsize specified! Please be sure that the CFD-condition holds!"
        end
        new(time, steps, stepsize, safety_factor)
    end
end

function Base.show(io::IO, @nospecialize(vv::VelocityVerlet))
    print(io, typeof(vv))
    fields = Vector{Symbol}()
    for field in fieldnames(typeof(vv))
        value = getfield(vv, field)
        if value > 0
            push!(fields, field)
        end
    end
    print(io, msg_fields_in_brackets(vv, Tuple(fields)))
    return nothing
end

function Base.show(io::IO, ::MIME"text/plain", @nospecialize(vv::VelocityVerlet))
    if get(io, :compact, false)
        show(io, vv)
    else
        println(io, typeof(vv), ":")
        fields = Vector{Symbol}()
        for field in fieldnames(typeof(vv))
            value = getfield(vv, field)
            if value > 0
                push!(fields, field)
            end
        end
        print(io, msg_fields(vv, Tuple(fields)))
    end
    return nothing
end

function init_time_solver!(vv::VelocityVerlet, dh::AbstractDataHandler)
    if vv.Δt < 0
        vv.Δt = calc_stable_timestep(dh, vv.safety_factor)
    end
    if vv.end_time < 0
        vv.end_time = vv.n_steps * vv.Δt
    elseif vv.n_steps < 0
        vv.n_steps = vv.end_time ÷ vv.Δt + 1
    end
    velocity_verlet_check(vv)
    return nothing
end

function velocity_verlet_check(vv::VelocityVerlet)
    if vv.end_time < 0
        error("`end_time` of VelocityVerlet smaller than zero!\n")
    end
    if vv.n_steps < 0
        error("`n_steps` of VelocityVerlet smaller than zero!\n")
    end
    if vv.Δt < 0
        error("`Δt` of VelocityVerlet smaller than zero!\n")
    end
    if !(0 < vv.safety_factor < 1)
        error("`safety_factor` of VelocityVerlet has invalid value!\n")
    end
    return nothing
end

function calc_stable_timestep(dh::AbstractDataHandler, safety_factor::Float64)
    throw(MethodError(calc_stable_timestep, dh, safety_factor))
end

function calc_stable_timestep(dh::ThreadsBodyDataHandler, safety_factor::Float64)
    Δt = zeros(length(dh.chunks))
    @threads :static for chunk_id in eachindex(dh.chunks)
        Δt[chunk_id] = calc_timestep(dh.chunks[chunk_id])
    end
    return minimum(Δt) * safety_factor
end

function calc_stable_timestep(dh::ThreadsMultibodyDataHandler, safety_factor::Float64)
    Δt = minimum(calc_stable_timestep(bdh, safety_factor) for bdh in each_body_dh(dh))
    return Δt
end

function calc_stable_timestep(dh::MPIBodyDataHandler, safety_factor::Float64)
    _Δt = calc_timestep(dh.chunk)
    Δt = MPI.Allreduce(_Δt, MPI.MIN, mpi_comm())
    return Δt * safety_factor
end

function calc_timestep(chunk::AbstractBodyChunk)
    isempty(each_point_idx(chunk)) && return Inf
    Δt = fill(Inf, length(each_point_idx(chunk)))
    for point_id in each_point_idx(chunk)
        pp = get_params(chunk, point_id)
        Δt[point_id] = calc_timestep_point(chunk.system, pp, point_id)
    end
    return minimum(Δt)
end

function solve!(dh::AbstractDataHandler, vv::VelocityVerlet, options::AbstractJobOptions)
    export_reference_results(dh, options)
    Δt = vv.Δt
    Δt½ = 0.5 * vv.Δt
    if mpi_isroot()
        p = Progress(vv.n_steps; dt=1, desc="TIME INTEGRATION LOOP", color=:normal,
                     barlen=40, enabled=progress_bars())
    end
    for n in 1:vv.n_steps
        verlet_timestep!(dh, options, Δt, Δt½, n)
        mpi_isroot() && next!(p)
    end
    mpi_isroot() && finish!(p)
    return dh
end

function verlet_timestep!(dh::AbstractThreadsBodyDataHandler, options::AbstractJobOptions,
                          Δt::Float64, Δt½::Float64, n::Int)
    t = n * Δt
    @threads :static for chunk_id in eachindex(dh.chunks)
        chunk = dh.chunks[chunk_id]
        update_vel_half!(chunk, Δt½)
        apply_boundary_conditions!(chunk, t)
        update_disp_and_pos!(chunk, Δt)
    end
    @threads :static for chunk_id in eachindex(dh.chunks)
        exchange_loc_to_halo!(dh, chunk_id)
        calc_force_density!(dh.chunks[chunk_id])
    end
    @threads :static for chunk_id in eachindex(dh.chunks)
        exchange_halo_to_loc!(dh, chunk_id)
        chunk = dh.chunks[chunk_id]
        calc_damage!(chunk)
        update_acc_and_vel!(chunk, Δt½)
        export_results(dh, options, chunk_id, n, t)
    end
    return nothing
end

function verlet_timestep!(dh::AbstractThreadsMultibodyDataHandler,
                          options::AbstractJobOptions, Δt::Float64, Δt½::Float64, n::Int)
    t = n * Δt
    for body_idx in each_body_idx(dh)
        body_dh = get_body_dh(dh, body_idx)
        @threads :static for chunk_id in eachindex(body_dh.chunks)
            chunk = body_dh.chunks[chunk_id]
            update_vel_half!(chunk, Δt½)
            apply_boundary_conditions!(chunk, t)
            update_disp_and_pos!(chunk, Δt)
        end
        @threads :static for chunk_id in eachindex(body_dh.chunks)
            exchange_loc_to_halo!(body_dh, chunk_id)
            calc_force_density!(body_dh.chunks[chunk_id])
        end
    end
    update_caches!(dh)
    calc_contact_force_densities!(dh)
    for body_idx in each_body_idx(dh)
        body_dh = get_body_dh(dh, body_idx)
        body_name = get_body_name(dh, body_idx)
        @threads :static for chunk_id in eachindex(body_dh.chunks)
            exchange_halo_to_loc!(body_dh, chunk_id)
            chunk = body_dh.chunks[chunk_id]
            calc_damage!(chunk)
            update_acc_and_vel!(chunk, Δt½)
            export_results(body_dh, options, chunk_id, n, t)
        end
    end
    return nothing
end

function verlet_timestep!(dh::AbstractMPIBodyDataHandler, options::AbstractJobOptions,
                          Δt::Float64, Δt½::Float64, n::Int)
    t = n * Δt
    chunk = dh.chunk
    @timeit_debug TO "update_vel_half!" update_vel_half!(chunk, Δt½)
    @timeit_debug TO "apply_boundary_conditions!" apply_boundary_conditions!(chunk, t)
    @timeit_debug TO "update_disp_and_pos!" update_disp_and_pos!(chunk, Δt)
    @timeit_debug TO "exchange_loc_to_halo!" exchange_loc_to_halo!(dh)
    @timeit_debug TO "calc_force_density!" calc_force_density!(chunk)
    @timeit_debug TO "exchange_halo_to_loc!" exchange_halo_to_loc!(dh)
    @timeit_debug TO "calc_damage!" calc_damage!(chunk)
    @timeit_debug TO "update_acc_and_vel!" update_acc_and_vel!(chunk, Δt½)
    @timeit_debug TO "export_results" export_results(dh, options, n, t)
    return nothing
end

function update_vel_half!(b::AbstractBodyChunk, Δt½::Float64)
    _update_vel_half!(b.storage.velocity_half, b.storage.velocity, b.storage.acceleration,
                      Δt½, each_point_idx(b))
    return nothing
end

function _update_vel_half!(velocity_half, velocity, acceleration, Δt½, each_point)
    for i in each_point
        velocity_half[1, i] = velocity[1, i] + acceleration[1, i] * Δt½
        velocity_half[2, i] = velocity[2, i] + acceleration[2, i] * Δt½
        velocity_half[3, i] = velocity[3, i] + acceleration[3, i] * Δt½
    end
    return nothing
end

function update_disp_and_pos!(b::AbstractBodyChunk, Δt::Float64)
    _update_disp_and_pos!(b.storage.displacement, b.storage.position,
                          b.storage.velocity_half, Δt, each_point_idx(b))
    return nothing
end

function _update_disp_and_pos!(displacement, position, velocity_half, Δt, each_point)
    for i in each_point
        u_x = velocity_half[1, i] * Δt
        u_y = velocity_half[2, i] * Δt
        u_z = velocity_half[3, i] * Δt
        displacement[1, i] += u_x
        displacement[2, i] += u_y
        displacement[3, i] += u_z
        position[1, i] += u_x
        position[2, i] += u_y
        position[3, i] += u_z
    end
    return nothing
end

function update_acc_and_vel!(b::AbstractBodyChunk, Δt½::Float64)
    _update_acc_and_vel!(b.storage, b.paramsetup, Δt½, each_point_idx(b))
    return nothing
end

function _update_acc_and_vel!(s::AbstractStorage, paramhandler::AbstractParameterHandler,
                              Δt½::Float64, each_point)
    for point_id in each_point
        param = get_params(paramhandler, point_id)
        _update_acc!(s.acceleration, s.b_int, s.b_ext, param.rho, point_id)
        _update_vel!(s.velocity, s.velocity_half, s.acceleration, Δt½, point_id)
    end
    return nothing
end

function _update_acc_and_vel!(s::AbstractStorage, params::AbstractPointParameters,
                              Δt½::Float64, each_point)
    for point_id in each_point
        _update_acc!(s.acceleration, s.b_int, s.b_ext, params.rho, point_id)
        _update_vel!(s.velocity, s.velocity_half, s.acceleration, Δt½, point_id)
    end
    return nothing
end

function _update_acc!(acceleration, b_int, b_ext, rho, i)
    acceleration[1, i] = (b_int[1, i] + b_ext[1, i]) / rho
    acceleration[2, i] = (b_int[2, i] + b_ext[2, i]) / rho
    acceleration[3, i] = (b_int[3, i] + b_ext[3, i]) / rho
    return nothing
end

function _update_vel!(velocity, velocity_half, acceleration, Δt½, i)
    velocity[1, i] = velocity_half[1, i] + acceleration[1, i] * Δt½
    velocity[2, i] = velocity_half[2, i] + acceleration[2, i] * Δt½
    velocity[3, i] = velocity_half[3, i] + acceleration[3, i] * Δt½
    return nothing
end

function init_field_solver(::VelocityVerlet, system::AbstractSystem, ::Val{:position})
    return copy(system.position)
end

function init_field_solver(::VelocityVerlet, system::AbstractSystem, ::Val{:displacement})
    return zeros(3, get_n_loc_points(system))
end

function init_field_solver(::VelocityVerlet, system::AbstractSystem, ::Val{:velocity})
    return zeros(3, get_n_loc_points(system))
end

function init_field_solver(::VelocityVerlet, system::AbstractSystem, ::Val{:velocity_half})
    return zeros(3, get_n_loc_points(system))
end

function init_field_solver(::VelocityVerlet, system::AbstractSystem, ::Val{:acceleration})
    return zeros(3, get_n_loc_points(system))
end

function init_field_solver(::VelocityVerlet, system::AbstractSystem, ::Val{:b_int})
    return zeros(3, get_n_loc_points(system))
end

function init_field_solver(::VelocityVerlet, system::AbstractSystem, ::Val{:b_ext})
    return zeros(3, get_n_loc_points(system))
end

function req_point_data_fields_timesolver(::Type{<:VelocityVerlet})
    fields = (:position, :displacement, :velocity, :velocity_half, :acceleration, :b_int,
              :b_ext)
    return fields
end


function req_bond_data_fields_timesolver(::Type{<:VelocityVerlet})
    return ()
end

function req_data_fields_timesolver(::Type{<:VelocityVerlet})
    return ()
end

function log_timesolver(options::AbstractJobOptions, vv::VelocityVerlet)
    msg = "VELOCITY VERLET TIME SOLVER\n"
    msg *= msg_qty("number of time steps", vv.n_steps)
    msg *= msg_qty("time step size", vv.Δt)
    msg *= msg_qty("time step safety factor", vv.safety_factor)
    msg *= msg_qty("simulation time", vv.end_time)
    log_it(options, msg)
    return nothing
end

register_solver!(VelocityVerlet)
