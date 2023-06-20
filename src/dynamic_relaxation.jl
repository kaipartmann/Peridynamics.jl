mutable struct DynamicRelaxation <: AbstractTimeDiscretization
    n_steps::Int
    Δt::Float64
    Λ::Float64
    function DynamicRelaxation(n_steps::Int, Δt::Real=1; damping_factor::Float64=1.0)
        new(n_steps, Δt, damping_factor)
    end
end

function init_time_discretization!(::DynamicRelaxation, ::AbstractPDBody, ::PDMaterial)
    return nothing
end

function time_loop!(body::AbstractPDBody, dr::DynamicRelaxation, mat::PDMaterial,
                    bcs::Vector{<:AbstractBC}, ics::Vector{<:AbstractIC}, es::ExportSettings)
    apply_ics!(body, ics)
    if es.exportflag
        export_vtk(body, es.resultfile_prefix, 0, 0.0)
    end
    damping_matrix = zeros(3, body.n_points)
    for i in axes(damping_matrix, 2)
        damping_matrix[:,i] .= dr.Λ * 6 * mat[i].K * dr.Δt^2 / (1/3 * mat[i].δ^2)
    end
    velocity_half_old = zeros(Float64, (3, body.n_points))
    b_int_old = zeros(Float64, (3, body.n_points))
    p = Progress(dr.n_steps; dt=1, desc="Time integration... ", barlen=30, color=:normal,
                 enabled=!is_logging(stderr))
    for t in 1:dr.n_steps
        time = t * dr.Δt
        apply_bcs!(body, bcs, time)
        update_disp_and_position!(body, dr.Δt)
        compute_forcedensity!(body, mat)
        update_thread_cache!(body)
        calc_damage!(body)
        cn = calc_damping(body, damping_matrix, velocity_half_old, b_int_old, dr.Δt)
        if t == 1
            finite_difference_first_step!(body, damping_matrix, velocity_half_old,
                                          b_int_old, dr.Δt)
        else
            finite_difference!(body, damping_matrix, velocity_half_old, b_int_old, dr.Δt, cn)
        end
        if mod(t, es.exportfreq) == 0
            export_vtk(body, es.resultfile_prefix, t, time)
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
            if velocity_half_old[d, i] !== 0.0
                cn1 -= body.displacement[d, i]^2 * (body.b_int[d, i, 1] - b_int_old[d, i]) /
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
            body.velocity_half[d, i] = 0.5 * Δt / damping_matrix[d, i] * (body.b_int[d, i, 1] + body.b_ext[d, i])
            body.velocity[d, i] = 0.5 * (velocity_half_old[d, i] + body.velocity_half[d, i])
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
            body.velocity_half[d, i] = ((2 - cn * Δt) * velocity_half_old[d, i] + 2 * Δt / damping_matrix[d, i] *
                                       (body.b_int[d, i, 1] + body.b_ext[d, i])) / (2 + cn * Δt)
            body.velocity[d, i] = 0.5 * (velocity_half_old[d, i] + body.velocity_half[d, i])
            velocity_half_old[d, i] = body.velocity_half[d, i]
            b_int_old[d, i] = body.b_int[d, i, 1]
        end
    end
    return nothing
end