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

function update_acc_and_vel!(b::BodyChunk, Δt½::Float64)
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

function solve!(dh::ThreadsDataHandler, vv::VelocityVerlet, options::ExportOptions)
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

function solve_timestep!(dh::ThreadsDataHandler, options::ExportOptions, Δt::Float64,
                         Δt½::Float64, n::Int)
    t = n * Δt
    @threads :static for chunk_id in eachindex(dh.chunks)
        chunk = dh.chunks[chunk_id]
        update_vel_half!(chunk, Δt½)
        apply_bcs!(chunk, t)
        update_disp_and_pos!(chunk, Δt)
    end
    @threads :static for chunk_id in eachindex(dh.chunks)
        exchange_read_fields!(dh, chunk_id)
        calc_force_density!(dh.chunks[chunk_id])
    end
    @threads :static for chunk_id in eachindex(dh.chunks)
        exchange_write_fields!(dh, chunk_id)
        chunk = dh.chunks[chunk_id]
        calc_damage!(chunk)
        update_acc_and_vel!(chunk, Δt½)
        export_results(dh, options, chunk_id, n, t)
    end
    return nothing
end

function solve!(dh::MPIDataHandler, vv::VelocityVerlet, options::ExportOptions)
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

function solve_timestep!(dh::MPIDataHandler, options::ExportOptions, Δt::Float64,
                         Δt½::Float64, n::Int)
    t = n * Δt
    chunk = dh.chunk
    @timeit_debug TO "update_vel_half!" update_vel_half!(chunk, Δt½)
    @timeit_debug TO "apply_bcs!" apply_bcs!(chunk, t)
    @timeit_debug TO "update_disp_and_pos!" update_disp_and_pos!(chunk, Δt)
    @timeit_debug TO "exchange_read_fields!" exchange_read_fields!(dh)
    @timeit_debug TO "calc_force_density!" calc_force_density!(chunk)
    @timeit_debug TO "exchange_write_fields!" exchange_write_fields!(dh)
    @timeit_debug TO "calc_damage!" calc_damage!(chunk)
    @timeit_debug TO "update_acc_and_vel!" update_acc_and_vel!(chunk, Δt½)
    @timeit_debug TO "export_results" export_results(dh, options, n, t)
    return nothing
end
