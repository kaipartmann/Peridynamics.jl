"""
    VelocityVerlet <: AbstractTimeSolver

Procedure for calculating discrete time steps

# Fields

- `end_time::Float64`: time covered by the simulation
- `n_steps::Int`: number of time steps calculated
- `Δt::Float64`: size of discrete time steps
- `safety_factor::Float64`: safety factor for step size to ensure stability

---

Constructors:

    VelocityVerlet <: AbstractTimeSolver

# Keywords

- `time::Real=-1`: time covered by the simulation
- `steps::Int=-1`: number of time steps calculated
- `stepsize::Real=-1`: size of discrete time steps
- `safety_factor::Real=0.7`: safety factor for step size to ensure stability

# Throws

- error if time and number of steps are specified
- error if no time or number of steps are specified
- error if safety factor is not between 0 and 1

# Example

```julia-repl
julia> vv = VelocityVerlet(steps=2000)
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

function calc_timestep(b::AbstractBodyChunk)
    Δt = fill(Inf, length(each_point_idx(b.ch)))
    for point_id in each_point_idx(b.ch)
        pp = get_param(b, point_id)
        Δt[point_id] = _calc_timestep(b.system, pp, point_id)
    end
    return minimum(Δt)
end

function _calc_timestep(bd::BondSystem, pp::AbstractPointParameters, point_id::Int)
    dtsum = 0.0
    for bond_id in each_bond_idx(bd, point_id)
        bond = bd.bonds[bond_id]
        dtsum += bd.volume[bond.neighbor] * pp.bc / bond.length
    end
    return sqrt(2 * pp.rho / dtsum)
end

function update_vel_half!(b::AbstractBodyChunk, Δt½::Float64)
    _update_vel_half!(b.storage.velocity_half, b.storage.velocity, b.storage.acceleration, Δt½,
                      each_point_idx(b))
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
    _update_disp_and_pos!(b.storage.displacement, b.storage.position, b.storage.velocity_half, Δt,
                          each_point_idx(b))
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
    _update_acc_and_vel!(b.storage.acceleration, b.storage.b_int, b.storage.b_ext,
                         b.storage.velocity_half, b.storage.velocity, b.param.rho, Δt½,
                         each_point_idx(b))
    return nothing
end

function update_acc_and_vel!(b::MultiParamBodyChunk, Δt½::Float64)
    s = b.storage
    for i in each_point
        param = get_param(b, i)
        _update_acc!(s.acceleration, s.b_int, s.b_ext, param.rho, i)
        _update_vel!(s.velocity, s.velocity_half, s.acceleration, Δt½, i)
    end
    return nothing
end

function _update_acc_and_vel!(acc, b_int, b_ext, vel_half, vel, rho, Δt½, each_point)
    for i in each_point
        _update_acc!(acc, b_int, b_ext, rho, i)
        _update_vel!(vel, vel_half, acc, Δt½, i)
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

function solve!(dh::AbstractThreadsDataHandler, vv::VelocityVerlet,
                options::AbstractOptions)
    export_reference_results(dh, options)
    Δt = vv.Δt
    Δt½ = 0.5 * vv.Δt
    p = Progress(vv.n_steps; dt=1, desc="solve...", color=:normal, barlen=20,
                 enabled=progress_enabled())
    for n in 1:vv.n_steps
        solve_timestep!(dh, options, Δt, Δt½, n)
        next!(p)
    end
    finish!(p)
    return dh
end

function solve_timestep!(dh::AbstractThreadsDataHandler, options::AbstractOptions,
                         Δt::Float64, Δt½::Float64, n::Int)
    t = n * Δt
    @threads :static for chunk_id in eachindex(dh.chunks)
        chunk = dh.chunks[chunk_id]
        update_vel_half!(chunk, Δt½)
        apply_bcs!(chunk, t)
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

function solve!(dh::AbstractMPIDataHandler, vv::VelocityVerlet, options::AbstractOptions)
    export_reference_results(dh, options)
    Δt = vv.Δt
    Δt½ = 0.5 * vv.Δt
    if mpi_isroot()
        p = Progress(vv.n_steps; dt=1, desc="solve...", color=:normal, barlen=20)
    end
    for n in 1:vv.n_steps
        solve_timestep!(dh, options, Δt, Δt½, n)
        mpi_isroot() && next!(p)
    end
    mpi_isroot() && finish!(p)
    return dh
end

function solve_timestep!(dh::AbstractMPIDataHandler, options::AbstractOptions, Δt::Float64,
                         Δt½::Float64, n::Int)
    t = n * Δt
    chunk = dh.chunk
    @timeit_debug TO "update_vel_half!" update_vel_half!(chunk, Δt½)
    @timeit_debug TO "apply_bcs!" apply_bcs!(chunk, t)
    @timeit_debug TO "update_disp_and_pos!" update_disp_and_pos!(chunk, Δt)
    @timeit_debug TO "exchange_loc_to_halo!" exchange_loc_to_halo!(dh)
    @timeit_debug TO "calc_force_density!" calc_force_density!(chunk)
    @timeit_debug TO "exchange_halo_to_loc!" exchange_halo_to_loc!(dh)
    @timeit_debug TO "calc_damage!" calc_damage!(chunk)
    @timeit_debug TO "update_acc_and_vel!" update_acc_and_vel!(chunk, Δt½)
    @timeit_debug TO "export_results" export_results(dh, options, n, t)
    return nothing
end

function req_point_data_fields_timesolver(::Type{VelocityVerlet})
    fields = (:position, :displacement, :velocity, :velocity_half, :acceleration, :b_int,
             :b_ext)
    return fields
end

function req_data_fields_timesolver(::Type{VelocityVerlet})
    return ()
end
