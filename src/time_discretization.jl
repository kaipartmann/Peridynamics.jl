const SUPPORTED_TD_ALGS = [:verlet, :dynrelax]

@doc raw"""
    TimeDiscretization

Time discretization type for setting the number of timesteps and the timestep `Δt`.

# Fields
- `n_timesteps::Int`: number of time steps
- `Δt::Float64`: constant time step
- `alg::Symbol`: algorithm used for time integration. Possible values:
    - `:verlet`: Velocity verlet algorithm for explicit time integration
    - `:dynrelax`: Adaptive dynamic relaxation for quasistatic time integration

---
```julia
TimeDiscretization(n_timesteps::Int[, Δt::Real]; alg::Symbol=:verlet)
```

# Arguments
- `n_timesteps::Int`: number of time steps
- `Δt::Real`: optional specified time step

# Keywords
- `alg::Symbol`: optional specification of algorithm used for time integration. Possible
  values:
    - `:verlet` (default): Velocity verlet algorithm for explicit time integration
    - `:dynrelax`: Adaptive dynamic relaxation for quasistatic time integration
"""
mutable struct TimeDiscretization
    n_timesteps::Int
    Δt::Float64
    alg::Symbol
end

function TimeDiscretization(n_timesteps::Int; alg::Symbol=:verlet)
    if alg in (SUPPORTED_TD_ALGS)
        if alg == :dynrelax
            return TimeDiscretization(n_timesteps, 1.0, alg)
        else
            return TimeDiscretization(n_timesteps, -1.0, alg)
        end
    else
        msg = "input $alg for keyword argument alg not supported!\n"
        msg *= "Supported input: $SUPPORTED_TD_ALGS"
        throw(DomainError(alg, msg))
    end
end

function TimeDiscretization(n_timesteps::Int, Δt::Real; alg::Symbol=:verlet)
    if alg in (SUPPORTED_TD_ALGS)
        if alg == :dynrelax && !(Δt ≈ 1)
            @warn "for dynamic relaxation a time step of Δt = 1 is recommended!"
        end
        return TimeDiscretization(n_timesteps, Δt, alg)
    else
        msg = "input $alg for keyword argument alg not supported!\n"
        msg *= "Supported input: $SUPPORTED_TD_ALGS"
        throw(DomainError(alg, msg))
    end
end

function Base.show(io::IO, ::MIME"text/plain", td::TimeDiscretization)
    println(io, typeof(td), ":")
    println(" number of time steps: ", td.n_timesteps)
    println(" Δt: ", td.Δt)
    if td.alg == :verlet
    print(" algorithm: dynamic (Velocity-Verlet algorithm)")
    elseif td.alg == :dynrelax
        print(" algorithm: quasi-static (adaptive dynamic relaxation)")
    end
    return nothing
end


function calc_stable_timestep(body::AbstractPDBody, mat::PDMaterial)
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
        _Δt[tid] = 0.7 * minimum(timesteps[timesteps .> 0])
    end
    Δt = minimum(_Δt)
    return Δt
end

"""
    calc_stable_user_timestep(pc::PointCloud, mat::AbstractPDMaterial, Sf::Float64=0.7)

Function to determine the stable timestep for the specified point cloud.

# Arguments
- `pc::PointCloud`: point cloud
- `mat::AbstractPDMaterial`: material model
- `Sf::Float64`: safety factor for time step, default value `Sf = 0.7`

# Returns
- `Float64`: stable user timestep `Δt`
"""
function calc_stable_user_timestep(pc::PointCloud, mat::PDMaterial, Sf::Float64=0.7)
    owned_points = defaultdist(pc.n_points, nthreads())
    bond_data, _ = find_bonds(pc, mat, owned_points)
    owned_bonds = defaultdist(length(bond_data), nthreads())
    _Δt = zeros(Float64, nthreads())
    @inbounds @threads for tid in 1:nthreads()
        timesteps = zeros(Float64, pc.n_points)
        dtsum = zeros(Float64, (pc.n_points, nthreads()))
        for current_one_ni in owned_bonds[tid]
            (a, i, L, _) = bond_data[current_one_ni]
            dtsum[a, tid] += pc.volume[i] * 1 / L * 18 * mat.K / (π * mat.δ^4)
        end
        for a in owned_points[tid]
            dtsum[a, 1] = sum(@view dtsum[a, :])
            timesteps[a] = √(2 * mat.rho / dtsum[a, 1])
        end
        _Δt[tid] = Sf * minimum(timesteps[timesteps .> 0])
    end
    Δt = minimum(_Δt)
    return Δt
end


function velocity_verlet!(body::AbstractPDBody, sim::AbstractPDAnalysis)
    p = Progress(sim.td.n_timesteps;
        dt=1,
        desc="Time integration... ",
        barlen=30,
        color=:normal,
        enabled=!is_logging(stderr),
    )
    Δt½ = 0.5 * sim.td.Δt
    for t in 1:sim.td.n_timesteps
        time = t * sim.td.Δt
        update_velhalf!(body, Δt½)
        apply_bcs!(body, sim.bcs, time)
        update_disp_and_position!(body, sim.td.Δt)
        compute_forcedensity!(body, sim.mat)
        update_thread_cache!(body)
        calc_damage!(body)
        compute_equation_of_motion!(body, Δt½, sim.mat)
        if mod(t, sim.es.exportfreq) == 0
            export_vtk(body, sim.es.resultfile_prefix, t, time)
        end
        next!(p)
    end
    finish!(p)
    return nothing
end

function dynamic_relaxation_finite_difference!(
    body::AbstractPDBody, sim::AbstractPDAnalysis
)
    damping_matrix = fill(
        6 * sim.mat.K * sim.td.Δt^2 / (1/3 * sim.mat.δ^2),
        (3, body.n_points),
    )
    velocity_half_old = zeros(Float64, (3, body.n_points))
    b_int_old = zeros(Float64, (3, body.n_points))
    p = Progress(sim.td.n_timesteps;
        dt=1,
        desc="Time integration... ",
        barlen=30,
        color=:normal,
        enabled=!is_logging(stderr),
    )
    for t in 1:sim.td.n_timesteps
        time = t * sim.td.Δt
        apply_bcs!(body, sim.bcs, time)
        compute_forcedensity!(body, sim.mat)
        update_thread_cache!(body)
        calc_damage!(body)
        cn = calc_damping(body, damping_matrix, velocity_half_old, b_int_old, sim.td.Δt)
        if t == 1
            finite_difference_first_step!(
                body, damping_matrix, velocity_half_old, b_int_old, sim.td.Δt
            )
        else
            finite_difference!(
                body, damping_matrix, velocity_half_old, b_int_old, sim.td.Δt, cn
            )
        end
        if mod(t, sim.es.exportfreq) == 0
            export_vtk(body, sim.es.resultfile_prefix, t, time)
        end
        next!(p)
    end
    finish!(p)
    return nothing
end

function calc_damping(
    body::AbstractPDBody,
    damping_matrix::Matrix{Float64},
    velocity_half_old::Matrix{Float64},
    b_int_old::Matrix{Float64},
    Δt::Float64,
)
    cn1 = 0.0
    cn2 = 0.0
    for i in 1:body.n_points
        for d in 1:3
            body.b_int[d, i, 1] = sum(@view body.b_int[d, i, :])
            if velocity_half_old[d, i] !== 0.0
                cn1 -=
                    body.displacement[d, i]^2 * (body.b_int[d, i, 1] - b_int_old[d, i]) /
                    (damping_matrix[d, i] * Δt * velocity_half_old[d, i])
            end
            cn2 += body.displacement[d, i]^2
        end
    end
    if cn2 !== 0.0
        if cn1 / cn2 > 0.0
            cn = 2.0 * sqrt(cn1 / cn2)
        else
            cn = 0.0
        end
    else
        cn = 0.0
    end
    if cn > 2.0
        cn = 1.9
    end
    return cn
end

function finite_difference_first_step!(
    body::AbstractPDBody,
    damping_matrix::Matrix{Float64},
    velocity_half_old::Matrix{Float64},
    b_int_old::Matrix{Float64},
    Δt::Float64,
)
    @threads for i in 1:body.n_points
        for d in 1:3
            body.velocity_half[d, i] =
                0.5 * Δt / damping_matrix[d, i] * (body.b_int[d, i, 1] + body.b_ext[d, i])
            body.velocity[d, i] = 0.5 * (velocity_half_old[d, i] + body.velocity_half[d, i])
            body.displacement[d, i] += body.velocity_half[d, i] * Δt
            body.position[d, i] += body.velocity_half[d, i] * Δt
            velocity_half_old[d, i] = body.velocity_half[d, i]
            b_int_old[d, i] = body.b_int[d, i, 1]
        end
    end
    return nothing
end

function finite_difference!(
    body::AbstractPDBody,
    damping_matrix::Matrix{Float64},
    velocity_half_old::Matrix{Float64},
    b_int_old::Matrix{Float64},
    Δt::Float64,
    cn::Float64,
)
    @threads for i in 1:body.n_points
        for d in 1:3
            body.velocity_half[d, i] =
                (
                    (2 - cn * Δt) * velocity_half_old[d, i] +
                    2 * Δt / damping_matrix[d, i] *
                    (body.b_int[d, i, 1] + body.b_ext[d, i])
                ) / (2 + cn * Δt)
            body.velocity[d, i] = 0.5 * (velocity_half_old[d, i] + body.velocity_half[d, i])
            body.displacement[d, i] += body.velocity_half[d, i] * Δt
            body.position[d, i] += body.velocity_half[d, i] * Δt
            velocity_half_old[d, i] = body.velocity_half[d, i]
            b_int_old[d, i] = body.b_int[d, i, 1]
        end
    end
    return nothing
end

function update_velhalf!(body::AbstractPDBody, Δt½::Float64)
    @inbounds @threads for i in 1:body.n_points
        body.velocity_half[1, i] = body.velocity[1, i] + body.acceleration[1, i] * Δt½
        body.velocity_half[2, i] = body.velocity[2, i] + body.acceleration[2, i] * Δt½
        body.velocity_half[3, i] = body.velocity[3, i] + body.acceleration[3, i] * Δt½
    end
    return nothing
end

function update_disp_and_position!(body::AbstractPDBody, Δt::Float64)
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

function compute_equation_of_motion!(body::AbstractPDBody, Δt½::Float64, mat::PDMaterial)
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

function update_thread_cache!(body::AbstractPDBody)
    @threads for (i, tid) in body.single_tids
        body.b_int[1, i, 1] = body.b_int[1, i, tid]
        body.b_int[2, i, 1] = body.b_int[2, i, tid]
        body.b_int[3, i, 1] = body.b_int[3, i, tid]
        body.n_active_family_members[i, 1] = body.n_active_family_members[i, tid]
    end
    @threads for (i, tids) in body.multi_tids
        for tid in tids
            body.b_int[1, i, 1] += body.b_int[1, i, tid]
            body.b_int[2, i, 1] += body.b_int[2, i, tid]
            body.b_int[3, i, 1] += body.b_int[3, i, tid]
            body.n_active_family_members[i, 1] += body.n_active_family_members[i, tid]
        end
    end
end
