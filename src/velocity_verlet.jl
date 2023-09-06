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

function init_time_discretization!(pdp::PDProblem)
    if pdp.sp.td.Δt < 0
        Δt = calc_stable_timestep(pdp)
        pdp.sp.td.Δt = Δt
    end
    return nothing
end

function calc_stable_timestep(pdp::PDProblem)
    _Δt = zeros(Float64, pdp.sp.n_threads)
    @inbounds @threads :static for tid in 1:pdp.sp.n_threads
        timesteps = zeros(Float64, pdp.sp.pc.n_points)
        dtsum = zeros(Float64, (pdp.sp.pc.n_points, pdp.sp.n_threads))
        for current_one_ni in pdp.sp.owned_bonds[tid]
            (a, i, L, _) = pdp.sp.bond_data[current_one_ni]
            dtsum[a, tid] += pdp.sp.pc.volume[i] * 1 / L * 18 * pdp.sp.mat[a].K / (π * pdp.sp.mat[a].δ^4)
        end
        for a in pdp.sp.owned_points[tid]
            dtsum[a, 1] = sum(@view dtsum[a, :])
            timesteps[a] = √(2 * pdp.sp.mat[a].rho / dtsum[a, 1])
        end
        _Δt[tid] = pdp.sp.td.safety_factor * minimum(timesteps[timesteps .> 0])
    end
    Δt = minimum(_Δt)
    return Δt
end

function time_loop!(pdp::PDProblem)
    apply_ics!(pdp)
    if pdp.sp.es.exportflag
        export_vtk(pdp, 0, 0.0)
    end
    p = Progress(pdp.sp.td.n_steps; dt = 1, desc = "Time integration... ", barlen = 30,
                 color = :normal, enabled = !is_logging(stderr))
    Δt½ = 0.5 * pdp.sp.td.Δt
    for t in 1:pdp.sp.td.n_steps
        time = t * pdp.sp.td.Δt
        update_velhalf!(pdp.gs, Δt½)
        apply_bcs!(pdp, time)
        update_disp_and_position!(pdp.gs, pdp.sp.td.Δt)
        compute_forcedensity!(pdp)
        reduce_tls_to_gs!(pdp)
        calc_damage!(pdp)
        compute_equation_of_motion!(pdp, Δt½)
        if mod(t, pdp.sp.es.exportfreq) == 0
            export_vtk(pdp, t, time)
        end
        next!(p)
    end
    finish!(p)
    return nothing
end

@timeit TO function update_velhalf!(gs::GlobalStorage, Δt½::Float64)
    @inbounds @threads :static for i in axes(gs.velocity_half, 2)
        gs.velocity_half[1, i] = gs.velocity[1, i] + gs.acceleration[1, i] * Δt½
        gs.velocity_half[2, i] = gs.velocity[2, i] + gs.acceleration[2, i] * Δt½
        gs.velocity_half[3, i] = gs.velocity[3, i] + gs.acceleration[3, i] * Δt½
    end
    return nothing
end

@timeit TO function update_disp_and_position!(gs::GlobalStorage, Δt::Float64)
    @inbounds @threads :static for i in axes(gs.displacement, 2)
        gs.displacement[1, i] += gs.velocity_half[1, i] * Δt
        gs.displacement[2, i] += gs.velocity_half[2, i] * Δt
        gs.displacement[3, i] += gs.velocity_half[3, i] * Δt
        gs.position[1, i] += gs.velocity_half[1, i] * Δt
        gs.position[2, i] += gs.velocity_half[2, i] * Δt
        gs.position[3, i] += gs.velocity_half[3, i] * Δt
    end
    return nothing
end

@timeit TO function compute_equation_of_motion!(pdp::PDProblem, Δt½::Float64)
    @inbounds @threads :static for i in 1:pdp.sp.pc.n_points
        pdp.gs.acceleration[1, i] = (pdp.gs.b_int[1, i, 1] + pdp.gs.b_ext[1, i]) / pdp.sp.mat[i].rho
        pdp.gs.acceleration[2, i] = (pdp.gs.b_int[2, i, 1] + pdp.gs.b_ext[2, i]) / pdp.sp.mat[i].rho
        pdp.gs.acceleration[3, i] = (pdp.gs.b_int[3, i, 1] + pdp.gs.b_ext[3, i]) / pdp.sp.mat[i].rho
        pdp.gs.velocity[1, i] = pdp.gs.velocity_half[1, i] + pdp.gs.acceleration[1, i] * Δt½
        pdp.gs.velocity[2, i] = pdp.gs.velocity_half[2, i] + pdp.gs.acceleration[2, i] * Δt½
        pdp.gs.velocity[3, i] = pdp.gs.velocity_half[3, i] + pdp.gs.acceleration[3, i] * Δt½
    end
    return nothing
end
