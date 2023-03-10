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
    mat::PDMaterial
    precracks::Vector{PreCrack}
    bcs::Vector{<:AbstractBC}
    ics::Vector{<:AbstractIC}
    calc_timestep::Bool
    function BodySetup(
        pc::PointCloud,
        mat::PDMaterial;
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
        for i in 1:sim.n_bodies
            for precrack in sim.body_setup[i].precracks
                define_precrack!(bodies[i], precrack)
            end
            update_thread_cache_contact!(bodies[i])
            calc_damage!(bodies[i])
        end
        if sim.td.Δt < 0.0
            _Δt = Vector{Float64}(undef, 0)
            for i in 1:sim.n_bodies
                if sim.body_setup[i].calc_timestep
                    push!(
                        _Δt,
                        calc_stable_timestep(
                            bodies[i],
                            sim.body_setup[i].mat,
                        ),
                    )
                end
            end
            sim.td.Δt = minimum(_Δt)
        end
        log_simsetup(sim)
        for i in 1:sim.n_bodies
            apply_ics!(bodies[i], sim.body_setup[i].ics)
        end
        if sim.es.exportflag
            for i in 1:sim.n_bodies
                export_vtk(bodies[i], string(sim.es.resultfile_prefix, "_b", i), 0, 0.0)
            end
        end
        velocity_verlet!(bodies, sim)
    end
    log_closesim(bodies, timingsummary, sim.es)
    return bodies
end

function log_describe_sim(sim::PDContactAnalysis)
    msg = "Peridynamic multibody contact analysis: " * sim.name * "\n"
    for i in 1:sim.n_bodies
        msg *= @sprintf("Body %d:\n", i)
        msg *= "Geometry:\n"
        msg *= describe_geometry(sim.body_setup[i].pc)
        msg *= @sprintf("Material Parameters: %49s\n", eltype(sim.body_setup[i].mat))
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

function log_describe_interactions(bodies::Vector{<:AbstractPDBody}, es::ExportSettings)
    msg = "Interactions:\n"
    for i in 1:length(bodies)
        msg *= @sprintf "Body %d:\n" i
        msg *= describe_interactions(bodies[i])
    end
    msg *= @sprintf "Total memory used for all parts [MB]:   %30g\n" Base.summarysize(
        bodies
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
    bodies::Vector{<:AbstractPDBody}, timingsummary::NamedTuple, es::ExportSettings
)
    msg = @sprintf "✓ Simulation completed after %g seconds\nResults:\n" timingsummary.time
    for i in 1:length(bodies)
        msg *= @sprintf "Body %d:\n" i
        msg *= log_displacement_and_damage(bodies[i])
    end
    print(msg)
    if es.exportflag
        open(es.logfile, "a") do io
            write(io, msg)
        end
    end
    return nothing
end

function velocity_verlet!(bodies::Vector{<:AbstractPDBody}, sim::PDContactAnalysis)
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
        for i in 1:sim.n_bodies
            update_velhalf!(bodies[i], Δt½)
            apply_bcs!(bodies[i], sim.body_setup[i].bcs, time)
            update_disp_and_position!(bodies[i], sim.td.Δt)
            compute_forcedensity!(bodies[i], sim.body_setup[i].mat)
            update_thread_cache_contact!(bodies[i])
            calc_damage!(bodies[i])
        end
        compute_contactforcedensity!(bodies, sim)
        for i in 1:sim.n_bodies
            compute_equation_of_motion_contact!(bodies[i], Δt½, sim.body_setup[i].mat.rho)
        end
        if mod(t, sim.es.exportfreq) == 0
            for i in 1:sim.n_bodies
                export_vtk(bodies[i], string(sim.es.resultfile_prefix, "_b", i), t, time)
            end
        end
        next!(p)
    end
    finish!(p)
    return nothing
end

function compute_equation_of_motion_contact!(
    body::AbstractPDBody,
    Δt½::Float64,
    rho::Float64,
)
    @threads for i in 1:body.n_points
        body.b_int[1, i, 1] = sum(@view body.b_int[1, i, :])
        body.b_int[2, i, 1] = sum(@view body.b_int[2, i, :])
        body.b_int[3, i, 1] = sum(@view body.b_int[3, i, :])
        body.acceleration[1, i] = (body.b_int[1, i, 1] + body.b_ext[1, i]) / rho
        body.acceleration[2, i] = (body.b_int[2, i, 1] + body.b_ext[2, i]) / rho
        body.acceleration[3, i] = (body.b_int[3, i, 1] + body.b_ext[3, i]) / rho
        body.velocity[1, i] = body.velocity_half[1, i] + body.acceleration[1, i] * Δt½
        body.velocity[2, i] = body.velocity_half[2, i] + body.acceleration[2, i] * Δt½
        body.velocity[3, i] = body.velocity_half[3, i] + body.acceleration[3, i] * Δt½
    end
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
        @threads for tid in 1:nthreads()
            for i in body_a.owned_points[tid]
                for j in 1:body_b.n_points
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

function update_thread_cache_contact!(body::AbstractPDBody)
    @threads for (i, tid) in body.single_tids
        body.n_active_family_members[i, 1] = body.n_active_family_members[i, tid]
    end
    @threads for (i, tids) in body.multi_tids
        for tid in tids
            body.n_active_family_members[i, 1] += body.n_active_family_members[i, tid]
        end
    end
end
