"""
    Contact

Contact definition.

# Fields
- `body_id_set::Tuple{Int,Int}`: bodies for which contact is defined
- `search_radius::Float64`: search radius for contact definition to activate
- `spring_constant::Float64`: spring constant used for contact force calculation, default
  value: `spring_constant = 1.0e12`
"""
struct Contact
    body_id_set::Tuple{Int,Int}
    search_radius::Float64
    spring_constant::Float64
    function Contact(
        body_id_set::Tuple{Int,Int}, search_radius::Float64, spring_constant::Float64=1.0e12
    )
        return new(body_id_set, search_radius, spring_constant)
    end
end

"""
    BodySetup

Setup of multiple bodies for `PDContactAnalysis`.

# Fields:
- `pc::PointCloud`: point cloud
- `mat::AbstractPDMaterial`: material model
- `precracks::Vector{PreCrack}`: predefined cracks
- `bcs::Vector{<:AbstractBC}`: boundary conditions
- `ics::Vector{<:AbstractIC}`: initial conditions
- `calc_timestep::Bool`: use body for time step calculation
"""
struct BodySetup
    pc::PointCloud
    mat::AbstractPDMaterial
    precracks::Vector{PreCrack}
    bcs::Vector{<:AbstractBC}
    ics::Vector{<:AbstractIC}
    calc_timestep::Bool
    function BodySetup(
        pc::PointCloud,
        mat::AbstractPDMaterial;
        precracks::Vector{PreCrack}=Vector{PreCrack}(),
        bcs::Vector{<:AbstractBC}=Vector{AbstractBC}(),
        ics::Vector{<:AbstractIC}=Vector{AbstractIC}(),
        calc_timestep::Bool=true,
    )
        return new(pc, mat, precracks, bcs, ics, calc_timestep)
    end
end

"""
    PDContactAnalysis <: AbstractPDAnalysis

Peridynamic contact analysis.

# Fields
- `name::String`: simulation name
- `body_setup::Vector{BodySetup}`: bodies used in this simulation
- `n_bodies::Int`: number of bodies
- `contact::Vector{Contact}`: contact definitions
- `td::TimeDiscretization`: time discretization
- `es::ExportSettings`: export settings
"""
struct PDContactAnalysis <: AbstractPDAnalysis
    name::String
    body_setup::Vector{BodySetup}
    n_bodies::Int
    contact::Vector{Contact}
    td::TimeDiscretization
    es::ExportSettings
    function PDContactAnalysis(;
        name::String,
        body_setup::Vector{BodySetup},
        contact::Vector{Contact},
        td::TimeDiscretization,
        es::ExportSettings,
    )
        es.resultfile_prefix = joinpath(es.path, name)
        es.logfile = es.resultfile_prefix * ".log"
        if !es.exportflag
            es.exportfreq = td.n_timesteps + 1
        end
        if td.alg == :dynrelax
            error("Dynamic relaxation not allowed for contact analysis! Check settings!")
        end
        return new(name, body_setup, length(body_setup), contact, td, es)
    end
end

function submit(sim::PDContactAnalysis)
    timingsummary = @timed begin
        log_header(sim.es)
        log_describe_sim(sim)
        bodies = [create_simmodel(ps.mat, ps.pc) for ps in sim.body_setup]
        log_describe_interactions(bodies, sim.es)
        for i in 1:(sim.n_bodies)
            for precrack in sim.body_setup[i].precracks
                define_precrack!(bodies[i], precrack)
            end
            calc_damage!(bodies[i])
        end
        if sim.td.Δt < 0.0
            _Δt = Vector{Float64}(undef, 0)
            for i in 1:(sim.n_bodies)
                if sim.body_setup[i].calc_timestep
                    push!(
                        _Δt,
                        calc_stable_timestep(
                            bodies[i],
                            sim.body_setup[i].mat.rho,
                            sim.body_setup[i].mat.K,
                            sim.body_setup[i].mat.δ,
                        ),
                    )
                end
            end
            sim.td.Δt = minimum(_Δt)
        end
        log_simsetup(sim)
        for i in 1:(sim.n_bodies)
            apply_ics!(bodies[i], sim.body_setup[i].ics)
        end
        if sim.es.exportflag
            for i in 1:(sim.n_bodies)
                export_results(bodies[i], string(sim.es.resultfile_prefix, "_b", i), 0, 0.0)
            end
        end
        velocity_verlet!(bodies, sim)
    end
    log_closesim(bodies, timingsummary, sim.es)
    return bodies
end

function log_describe_sim(sim::PDContactAnalysis)
    msg = "Peridynamic multibody contact analysis: " * sim.name * "\nMaterial parameters:\n"
    for i in 1:(sim.n_bodies)
        msg *= @sprintf("Body %d:\n", i)
        msg *= @sprintf(
            "  - Number of material points [-]:      %30g\n",
            sim.body_setup[i].pc.n_points
        )
        msg *= describe_mat(sim.body_setup[i].mat)
    end
    print(msg)
    if sim.es.exportflag
        open(sim.es.logfile, "a") do io
            write(io, msg)
        end
    end
    return nothing
end

function log_describe_interactions(parts::Vector{<:AbstractPDBody}, es::ExportSettings)
    msg = "Interactions:\n"
    for i in 1:length(parts)
        msg *= @sprintf "Body %d:\n" i
        msg *= describe_interactions(parts[i])
    end
    msg *= @sprintf "Total memory used for all parts [MB]:   %30g\n" Base.summarysize(
        parts
    ) * 1e-6
    print(msg)
    if es.exportflag
        open(es.logfile, "a") do io
            write(io, msg)
        end
    end
    return nothing
end

function log_closesim(
    parts::Vector{<:AbstractPDBody}, timingsummary::NamedTuple, es::ExportSettings
)
    msg = @sprintf "✓ Simulation completed after %g seconds\nResults:\n" timingsummary.time
    for i in 1:length(parts)
        msg *= @sprintf "Body %d:\n" i
        msg *= @sprintf "  - Max. abs. x-displacement [m]:       %30g\n" abs(
            maximum(parts[i].displacement[1, :])
        )
        msg *= @sprintf "  - Max. abs. y-displacement [m]:       %30g\n" abs(
            maximum(parts[i].displacement[2, :])
        )
        msg *= @sprintf "  - Max. abs. z-displacement [m]:       %30g\n" abs(
            maximum(parts[i].displacement[3, :])
        )
        msg *= @sprintf "  - Max. damage [-]:                    %30g\n" maximum(
            parts[i].damage
        )
    end
    print(msg)
    if es.exportflag
        open(es.logfile, "a") do io
            write(io, msg)
        end
    end
    return nothing
end

function velocity_verlet!(body::Vector{<:AbstractPDBody}, sim::PDContactAnalysis)
    p = Progress(
        sim.td.n_timesteps; dt=1, desc="Time integration... ", barlen=30, color=:normal
    )
    Δt½ = 0.5 * sim.td.Δt
    for t in 1:(sim.td.n_timesteps)
        time = t * sim.td.Δt
        for i in 1:(sim.n_bodies)
            update_velhalf!(body[i], Δt½)
            apply_bcs!(body[i], sim.body_setup[i].bcs, time)
            update_disp_and_position!(body[i], sim.td.Δt)
            compute_forcedensity!(body[i], sim.body_setup[i].mat)
            calc_damage!(body[i])
        end
        compute_contactforcedensity!(body, sim)
        for i in 1:(sim.n_bodies)
            compute_equation_of_motion!(body[i], Δt½, sim.body_setup[i].mat.rho)
        end
        if mod(t, sim.es.exportfreq) == 0
            for i in 1:(sim.n_bodies)
                export_results(body[i], string(sim.es.resultfile_prefix, "_b", i), t, time)
            end
        end
        next!(p)
    end
    finish!(p)
    return nothing
end

function compute_contactforcedensity!(
    bodies::Vector{<:AbstractPDBody}, sim::PDContactAnalysis
)
    for contact in sim.contact
        ii, jj = contact.body_id_set
        body_a = bodies[ii]
        body_b = bodies[jj]
        r = contact.search_radius
        C = 9 * contact.spring_constant / (π * sim.body_setup[ii].mat.δ^4)
        @threads :static for _ in 1:nthreads()
            tid = threadid()
            for i in body_a.owned_points[tid]
                for j in 1:(body_b.n_points)
                    ϑx = body_b.position[1, j] - body_a.position[1, i]
                    ϑy = body_b.position[2, j] - body_a.position[2, i]
                    ϑz = body_b.position[3, j] - body_a.position[3, i]
                    η = sqrt(ϑx * ϑx + ϑy * ϑy + ϑz * ϑz)
                    if η < r
                        temp::Float64 = C * (r - η) / η
                        body_a.b_int[1, i, tid] -= temp * ϑx * body_b.volume[j]
                        body_a.b_int[2, i, tid] -= temp * ϑy * body_b.volume[j]
                        body_a.b_int[3, i, tid] -= temp * ϑz * body_b.volume[j]
                        body_b.b_int[1, j, tid] += temp * ϑx * body_a.volume[i]
                        body_b.b_int[2, j, tid] += temp * ϑy * body_a.volume[i]
                        body_b.b_int[3, j, tid] += temp * ϑz * body_a.volume[i]
                    end
                end
            end
        end
    end
    return nothing
end
