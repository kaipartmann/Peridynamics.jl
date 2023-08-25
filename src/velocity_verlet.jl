"""
    VelocityVerlet <: AbstractTimeDiscretization

Velocity Verlet algorithm for dynamic simulations.

# Fields
- `n_steps::Int`: number of time steps
- `Δt::Float64`: critical stable time step
- `safety_factor::Float64`: safety factor for time step

---
```julia
VelocityVerlet(n_steps::Int, Δt::Real = -1; safety_factor::Real = 0.7)
```
Only the number of time steps `n_steps` is needed. Specification of `Δt` and the
safety factor are optional.
"""
mutable struct VelocityVerlet <: AbstractTimeDiscretization
    n_steps::Int
    Δt::Float64
    safety_factor::Float64

    function VelocityVerlet(n_steps::Int, Δt::Real = -1; safety_factor::Real = 0.7)
        new(n_steps, Δt, safety_factor)
    end
end

function init_time_discretization!(vv::VelocityVerlet, body::AbstractPDBody,
                                   mat::PDMaterial)
    if vv.Δt < 0
        Δt = calc_stable_timestep(body, mat, vv.safety_factor)
        vv.Δt = Δt
    end
    return nothing
end

function calc_stable_timestep(body::AbstractPDBody, mat::PDMaterial, safety_factor::Float64)
    _Δt = zeros(Float64, body.n_threads)
    @inbounds @threads for tid in 1:body.n_threads
        timesteps = zeros(Float64, body.n_points)
        dtsum = zeros(Float64, (body.n_points, body.n_threads))
        for current_one_ni in body.owned_bonds[tid]
            (a, i, L, _) = body.bond_data[current_one_ni]
            dtsum[a, tid] += body.volume[i] * 1 / L * 18 * mat[a].K / (π * mat[a].δ^4)
        end
        for a in body.owned_points[tid]
            dtsum[a, 1] = sum(@view dtsum[a, :])
            timesteps[a] = √(2 * mat[a].rho / dtsum[a, 1])
        end
        _Δt[tid] = safety_factor * minimum(timesteps[timesteps .> 0])
    end
    Δt = minimum(_Δt)
    return Δt
end

function time_loop!(body::AbstractPDBody, vv::VelocityVerlet, mat::PDMaterial,
                    bcs::Vector{<:AbstractBC}, ics::Vector{<:AbstractIC},
                    es::ExportSettings)
    apply_ics!(body, ics)
    if es.exportflag
        export_vtk(body, es.resultfile_prefix, 0, 0.0)
    end
    p = Progress(vv.n_steps; dt = 1, desc = "Time integration... ", barlen = 30,
                 color = :normal, enabled = !is_logging(stderr))
    Δt½ = 0.5 * vv.Δt
    for t in 1:(vv.n_steps)
        time = t * vv.Δt
        update_velhalf!(body, Δt½)
        apply_bcs!(body, bcs, time)
        update_disp_and_position!(body, vv.Δt)
        compute_forcedensity!(body, mat)
        update_thread_cache!(body)
        calc_damage!(body)
        compute_equation_of_motion!(body, Δt½, mat)
        if mod(t, es.exportfreq) == 0
            export_vtk(body, es.resultfile_prefix, t, time)
        end
        next!(p)
    end
    finish!(p)
    return nothing
end

@timeit TO function update_velhalf!(body::AbstractPDBody, Δt½::Float64)
    @inbounds @threads for i in 1:body.n_points
        body.velocity_half[1, i] = body.velocity[1, i] + body.acceleration[1, i] * Δt½
        body.velocity_half[2, i] = body.velocity[2, i] + body.acceleration[2, i] * Δt½
        body.velocity_half[3, i] = body.velocity[3, i] + body.acceleration[3, i] * Δt½
    end
    return nothing
end

@timeit TO function update_disp_and_position!(body::AbstractPDBody, Δt::Float64)
    @inbounds @threads for i in 1:body.n_points
        body.displacement[1, i] += body.velocity_half[1, i] * Δt
        body.displacement[2, i] += body.velocity_half[2, i] * Δt
        body.displacement[3, i] += body.velocity_half[3, i] * Δt
        body.position[1, i] += body.velocity_half[1, i] * Δt
        body.position[2, i] += body.velocity_half[2, i] * Δt
        body.position[3, i] += body.velocity_half[3, i] * Δt
    end
    return nothing
end

@timeit TO function compute_equation_of_motion!(body::AbstractPDBody, Δt½::Float64, mat::PDMaterial)
    @inbounds @threads for i in 1:body.n_points
        body.acceleration[1, i] = (body.b_int[1, i, 1] + body.b_ext[1, i]) / mat[i].rho
        body.acceleration[2, i] = (body.b_int[2, i, 1] + body.b_ext[2, i]) / mat[i].rho
        body.acceleration[3, i] = (body.b_int[3, i, 1] + body.b_ext[3, i]) / mat[i].rho
        body.velocity[1, i] = body.velocity_half[1, i] + body.acceleration[1, i] * Δt½
        body.velocity[2, i] = body.velocity_half[2, i] + body.acceleration[2, i] * Δt½
        body.velocity[3, i] = body.velocity_half[3, i] + body.acceleration[3, i] * Δt½
    end
    return nothing
end
