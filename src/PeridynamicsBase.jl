@doc raw"""
    PointCloud

Peridynamic spatial discretization with material points defining a point cloud.

# Fields
- `n_points::Int`: number of material points
- `position::Matrix{Float64}`: coordinates of points in reference configuration
- `volume::Vector{Float64}`: material point volumes
- `failure_flag::BitVector`: if failure of point is possible: element=`true`
- `radius::Vector{Float64}`: radius of the material point sphere

---
```julia
PointCloud(position, volume[, point_sets])
```

Create a `PointCloud` by specifying position and volume of all material points. The radius
of the point sphere is derived from the volume with the function [`sphere_radius`](@ref) and
the `failure_flag` set to `true` for all points.

# Arguments
- `position::Matrix{Float64}`: the position of the material points ($3 \times N$-matrix for
  $N$ material points)
- `volume::Vector{Float64}`: the volume of the material points
- `point_sets::Dict{String,Vector{Int}}=Dict{String,Vector{Int}}()`: optional point sets

# Returns
- `PointCloud`: point cloud with specified material point position and volume

# Examples

Creating a `PointCloud` with 4 manually defined points:
```julia-repl
julia> position = [
           0.0 1.0 0.0 0.0
           0.0 0.0 1.0 0.0
           0.0 0.0 0.0 1.0
       ]
3×4 Matrix{Float64}:
 0.0  1.0  0.0  0.0
 0.0  0.0  1.0  0.0
 0.0  0.0  0.0  1.0

julia> volume = [1.0, 1.0, 1.0, 1.0]
4-element Vector{Float64}:
 1.0
 1.0
 1.0
 1.0

julia> pc = PointCloud(position, volume)
4-points PointCloud

julia> pc.failure_flag
4-element BitVector:
 1
 1
 1
 1

julia> pc.radius
4-element Vector{Float64}:
 0.6203504908994001
 0.6203504908994001
 0.6203504908994001
 0.6203504908994001
```

---
```julia
PointCloud(W, L, H, Δx; center_x=0, center_y=0, center_z=0)
```

Generate a uniformly distributed `PointCloud` with width `W`, length `L`, height `H` and
point spacing `Δx`. Optional keyword arguments provide the possibility to set the center.

# Arguments
- `W::Real`: width of the `PointCloud`-block in x-direction
- `L::Real`: length of the `PointCloud`-block in y-direction
- `H::Real`: height of the `PointCloud`-block in z-direction
- `Δx::Real`: point spacing in x-, y- and z- direction
- `center_x::Real=0`: x-coordinate of the `PointCloud`-block center
- `center_y::Real=0`: y-coordinate of the `PointCloud`-block center
- `center_z::Real=0`: z-coordinate of the `PointCloud`-block center

# Returns
- `PointCloud`: point cloud with with `W`, length `L`, height `H` and
  point spacing `Δx`

# Examples

Cube with side length 1 and point spacing $Δx = 0.1$:
```julia-repl
julia> PointCloud(1, 1, 1, 0.1)
1000-points PointCloud
```
"""
struct PointCloud
    n_points::Int
    position::Matrix{Float64}
    volume::Vector{Float64}
    failure_flag::BitVector
    radius::Vector{Float64}
    point_sets::Dict{String,Vector{Int}}
end

function PointCloud(
    position::Matrix{Float64},
    volume::Vector{Float64},
    point_sets::Dict{String,Vector{Int}}=Dict{String,Vector{Int}}(),
)
    if size(position, 1) !== 3 || size(position, 2) !== length(volume)
        throw(DimensionMismatch("size of position: $(size(position)) !== (3, n_points)"))
    end
    n_points = length(volume)
    radius = sphere_radius.(volume)
    failure_flag = BitVector(undef, n_points)
    failure_flag .= true
    return PointCloud(n_points, position, volume, failure_flag, radius, point_sets)
end

function PointCloud(
    W::Real,
    L::Real,
    H::Real,
    Δx::Real;
    center_x::Real=0,
    center_y::Real=0,
    center_z::Real=0,
)
    _gridX = range(; start=(-W + Δx) / 2, stop=(W - Δx) / 2, step=Δx)
    gridX = _gridX .- sum(_gridX) / length(_gridX)
    _gridY = range(; start=(-L + Δx) / 2, stop=(L - Δx) / 2, step=Δx)
    gridY = _gridY .- sum(_gridY) / length(_gridY)
    _gridZ = range(; start=(-H + Δx) / 2, stop=(H - Δx) / 2, step=Δx)
    gridZ = _gridZ .- sum(_gridZ) / length(_gridZ)
    positions = hcat(([x; y; z] for x in gridX for y in gridY for z in gridZ)...)
    if center_x !== 0
        positions[1, :] .+= center_x
    end
    if center_y !== 0
        positions[2, :] .+= center_y
    end
    if center_z !== 0
        positions[3, :] .+= center_z
    end
    n_points = size(positions, 2)
    volumes = fill(Δx^3, n_points)
    radius = sphere_radius.(volumes)
    failure_flag = BitVector(undef, n_points)
    failure_flag .= true
    point_sets = Dict{String,Vector{Int}}()
    return PointCloud(n_points, positions, volumes, failure_flag, radius, point_sets)
end

@doc raw"""
    sphere_radius(vol::T) where {T<:Real}

Calculate the radius $r$ of the sphere by equation
```math
r = \sqrt[3]{\frac{3 \; V}{4 \; \pi}}
```
with specified sphere volume $V$.

# Arguments
- `vol::T where {T<:Real}`: volume $V$ of the sphere

# Returns
- `Float64`: radius $r$ of sphere with volume `vol`
"""
function sphere_radius(vol::T) where {T<:Real}
    r = (3 * vol / (4 * π))^(1 / 3)
    return r
end

function Base.show(io::IO, ::MIME"text/plain", pc::PointCloud)
    println(io, pc.n_points, "-points ", typeof(pc))
    return nothing
end

"""
    PreCrack(point_id_set_a::Vector{Int}, point_id_set_b::Vector{Int})

Definition of an preexisting crack in the model. Points in `point_id_set_a` cannot have
interactions with points in `point_id_set_b`.

# Fields
- `point_id_set_a::Vector{Int}`: first point-id set
- `point_id_set_b::Vector{Int}`: second point-id set
"""
struct PreCrack
    point_id_set_a::Vector{Int}
    point_id_set_b::Vector{Int}
end

"""
    AbstractBC

Abstract type for boundary conditions.
"""
abstract type AbstractBC end

"""
    AbstractIC

Abstract type for initial conditions.
"""
abstract type AbstractIC end

@doc raw"""
    VelocityBC <: AbstractBC

Velocity boundary condition. The value of the velocity is calculated every time step with
the function `fun` and applied to the dimension `dim` of the points specified by
`point_id_set`.

# Fields
- `fun::Function`: function `f(t)` with current time `t` as argument, that calculates the
  value of the velocity for each timestep
- `point_id_set::Vector{Int}`: point-id set with all points for which the boundary condition
  gets applied every timestep
- `dim::Int`: dimension on which the boundary condition gets applied to. Possible values:
    - x-direction: `dim=1`
    - y-direction: `dim=2`
    - z-direction: `dim=3`

# Examples

The constant velocity $v = 0.1$ in $y$-direction gets applied to the first $10$ points:
```julia-repl
julia> VelocityBC(t -> 0.1, 1:10, 2)
VelocityBC(var"#7#8"(), [1, 2, 3, 4, 5, 6, 7, 8, 9, 10], 2)
```
"""
struct VelocityBC <: AbstractBC
    fun::Function
    point_id_set::Vector{Int}
    dim::Int
end

@doc raw"""
    ForceDensityBC <: AbstractBC

Force density boundary condition. The value of the force density is calculated every time
step with the function `fun` and applied to the dimension `dim` of the points specified by
`point_id_set`.

# Fields
- `fun::Function`: function `f(t)` with current time `t` as argument, that calculates the
  value of the force density for each timestep
- `point_id_set::Vector{Int}`: point-id set with all points for which the boundary condition
  gets applied every timestep
- `dim::Int`: dimension on which the boundary condition gets applied to. Possible values:
    - x-direction: `dim=1`
    - y-direction: `dim=2`
    - z-direction: `dim=3`
"""
struct ForceDensityBC <: AbstractBC
    fun::Function
    point_id_set::Vector{Int}
    dim::Int
end

@doc raw"""
    VelocityIC <: AbstractIC

Velocity initial condition. The value `val` of the velocity is applied as initial condition
to the dimension `dim` of the points specified by `point_id_set`.

# Fields
- `val::Float64`: value of the velocity
- `point_id_set::Vector{Int}`: point-id set with all points for which the initial condition
  gets applied to
- `dim::Int`: dimension on which the initial condition gets applied to. Possible values:
    - x-direction: `dim=1`
    - y-direction: `dim=2`
    - z-direction: `dim=3`
"""
struct VelocityIC <: AbstractIC
    val::Float64
    point_id_set::Vector{Int}
    dim::Int
end

@doc raw"""
    PosDepVelBC <: AbstractBC

Position dependent velocity boundary condition. The value of the force density is calculated
every time step with the function `fun` and applied to the dimension `dim` of the points
specified by `point_id_set`.

# Fields
- `fun::Function`: function `f(x, y, z, t)` with x-, y-, z-position of the point and current
  time `t` as arguments, calculates the value of the force density for each timestep
  dependent of the point position
- `point_id_set::Vector{Int}`: point-id set with all points for which the boundary condition
  gets applied to
- `dim::Int`: dimension on which the initial condition gets applied to. Possible values:
    - x-direction: `dim=1`
    - y-direction: `dim=2`
    - z-direction: `dim=3`
"""
struct PosDepVelBC <: AbstractBC
    fun::Function
    point_id_set::Vector{Int}
    dim::Int
end

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

"""
    ExportSettings

Export settings.

# Fields
- `path::String`: path where results will be saved
- `exportfreq::Int`: export frequency, will export every `exportfreq`-th timestep
- `resultfile_prefix::String`: prefix of the result-filename
- `logfile::String`: name of logfile
- `exportflag::Bool`: disable export for a simulation where saved results are not needed

---
```julia
ExportSettings([path::String, freq::Int])
```

Create `ExportSettings` only by `path` and `freq`. If no arguments are specified, the
`exportflag` will be set to `false` and export disabled.

# Arguments
- `path::String`: path where results will be saved
- `freq::Int`: export frequency
"""
mutable struct ExportSettings
    path::String
    exportfreq::Int
    resultfile_prefix::String
    logfile::String
    exportflag::Bool
end

ExportSettings() = ExportSettings("", 0, "", "", false)
ExportSettings(path::String, freq::Int) = ExportSettings(path, freq, "", "", true)

abstract type AbstractPDBody end

"""
    AbstractPDMaterial

Abstract type for a peridynamic material.
"""
abstract type AbstractPDMaterial end

"""
    AbstractPDAnalysis

Abstract type for a peridynamic analysis.
"""
abstract type AbstractPDAnalysis end

function find_bonds(pc::PointCloud, δ::Float64, owned_points::Vector{UnitRange{Int}})
    n_threads = nthreads()
    _bond_data = fill([(0, 0, 0.0, true)], n_threads)
    n_family_members = zeros(Int, pc.n_points)
    p = Progress(pc.n_points; dt=1, desc="Search bonds...     ", barlen=30, color=:normal)
    @threads for _ in 1:n_threads
        tid = threadid()
        local_bond_data = Vector{Tuple{Int,Int,Float64,Bool}}(undef, 0)
        sizehint!(local_bond_data, pc.n_points * 500)
        idst = 0.0
        num = 0
        fail = false
        for a in owned_points[tid]
            num = 0
            for i in 1:(pc.n_points)
                if a !== i
                    idst = sqrt(
                        (pc.position[1, i] - pc.position[1, a])^2 +
                        (pc.position[2, i] - pc.position[2, a])^2 +
                        (pc.position[3, i] - pc.position[3, a])^2,
                    )
                    if idst <= δ
                        num += 1
                        fail = pc.failure_flag[a] & pc.failure_flag[i]
                        push!(local_bond_data, (a, i, idst, fail))
                    end
                end
            end
            n_family_members[a] = num
            next!(p)
        end
        _bond_data[tid] = local_bond_data
    end
    finish!(p)
    bond_data = reduce(append!, _bond_data)
    return bond_data, n_family_members
end

function find_unique_bonds(pc::PointCloud, δ::Float64, owned_points::Vector{UnitRange{Int}})
    n_threads = nthreads()
    _bond_data = fill([(0, 0, 0.0, true)], n_threads)
    n_family_members = zeros(Int, pc.n_points)
    p = Progress(pc.n_points; dt=1, desc="Search bonds...     ", barlen=30, color=:normal)
    @threads for _ in 1:n_threads
        tid = threadid()
        local_bonds_data = Vector{Tuple{Int,Int,Float64,Bool}}(undef, 0)
        sizehint!(local_bonds_data, pc.n_points * 500)
        idst = 0.0
        num = 0
        fail = false
        for a in owned_points[tid]
            num = 0
            for i in 1:(pc.n_points)
                if a !== i
                    idst = sqrt(
                        (pc.position[1, i] - pc.position[1, a])^2 +
                        (pc.position[2, i] - pc.position[2, a])^2 +
                        (pc.position[3, i] - pc.position[3, a])^2,
                    )
                    if idst <= δ
                        num += 1
                        fail = pc.failure_flag[a] & pc.failure_flag[i]
                        push!(
                            local_bonds_data,
                            a < i ? (a, i, idst, fail) : (i, a, idst, fail),
                        )
                    end
                end
            end
            n_family_members[a] = num
            next!(p)
        end
        _bond_data[tid] = local_bonds_data
    end
    finish!(p)
    one_ni_data = unique(reduce(append!, _bond_data))
    return one_ni_data, n_family_members
end

function defaultdist(sz::Int, nc::Int)
    if sz >= nc
        chunk_size = div(sz, nc)
        remainder = rem(sz, nc)
        sidx = zeros(Int64, nc + 1)
        for i in 1:(nc + 1)
            sidx[i] += (i - 1) * chunk_size + 1
            if i <= remainder
                sidx[i] += i - 1
            else
                sidx[i] += remainder
            end
        end
        grid = fill(0:0, nc)
        for i in 1:nc
            grid[i] = sidx[i]:(sidx[i + 1] - 1)
        end
        return grid
    else
        sidx = [1:(sz + 1);]
        grid = fill(0:0, nc)
        for i in 1:sz
            grid[i] = sidx[i]:(sidx[i + 1] - 1)
        end
        return grid
    end
end

@doc raw"""
    PDSingleBodyAnalysis{T<:AbstractPDMaterial} <: AbstractPDAnalysis

Peridynamic single body analysis.

# Fields
- `name::String`: simulation name
- `pc::PointCloud`: point cloud
- `mat::T`: material model for the body
- `precracks::Vector{PreCrack}`: predefined cracks
- `bcs::Vector{<:AbstractBC}`: boundary conditions
- `ics::Vector{<:AbstractIC}`: initial conditions
- `td::TimeDiscretization`: time discretization
- `es::ExportSettings`: export settings
"""
struct PDSingleBodyAnalysis{T<:AbstractPDMaterial} <: AbstractPDAnalysis
    name::String
    pc::PointCloud
    mat::T
    precracks::Vector{PreCrack}
    bcs::Vector{<:AbstractBC}
    ics::Vector{<:AbstractIC}
    td::TimeDiscretization
    es::ExportSettings
    function PDSingleBodyAnalysis(;
        name::String,
        pc::PointCloud,
        mat::AbstractPDMaterial,
        precracks::Vector{PreCrack}=Vector{PreCrack}(),
        bcs::Vector{<:AbstractBC}=Vector{AbstractBC}(),
        ics::Vector{<:AbstractIC}=Vector{AbstractIC}(),
        td::TimeDiscretization,
        es::ExportSettings,
    )
        es.resultfile_prefix = joinpath(es.path, name)
        es.logfile = es.resultfile_prefix * ".log"
        if !es.exportflag
            es.exportfreq = td.n_timesteps + 1
        end
        return new{typeof(mat)}(name, pc, mat, precracks, bcs, ics, td, es)
    end
end

"""
    submit(sim::T) where {T<:AbstractPDAnalysis}

Submit a simulation job to determine the results of the specified analysis.

# Arguments
- `sim::T where {T<:AbstractPDAnalysis}`: simulation job
"""
function submit end

function submit(sim::PDSingleBodyAnalysis)
    timingsummary = @timed begin
        log_header(sim.es)
        log_describe_sim(sim.name, sim.pc, sim.mat, sim.es)
        bodies = create_simmodel(sim.mat, sim.pc)
        log_describe_interactions(bodies, sim.es)
        for precrack in sim.precracks
            define_precrack!(bodies, precrack)
        end
        calc_damage!(bodies)
        if sim.td.Δt < 0.0 && sim.td.alg !== :dynrelax
            sim.td.Δt = calc_stable_timestep(bodies, sim.mat.rho, sim.mat.K, sim.mat.δ)
        end
        log_simsetup(sim)
        apply_ics!(bodies, sim.ics)
        if sim.es.exportflag
            export_results(bodies, sim.es.resultfile_prefix, 0, 0.0)
        end
        if sim.td.alg == :verlet
            velocity_verlet!(bodies, sim)
        elseif sim.td.alg == :dynrelax
            dynamic_relaxation_finite_difference!(bodies, sim)
        end
    end
    log_closesim(bodies, timingsummary, sim.es)
    return bodies
end

function log_header(es::ExportSettings)
    msg = "="^70 * "\n"
    msg *= string("PERIDYNAMIC SIMULATION ON ", nthreads(), " THREADS\n")
    msg *= "="^70 * "\n"
    print(msg)
    if es.exportflag
        open(es.logfile, "w") do io
            return write(io, msg)
        end
    end
    return nothing
end

function log_describe_sim(
    simname::String, pc::PointCloud, mat::AbstractPDMaterial, es::ExportSettings
)
    msg = "Peridynamic single body analysis: " * simname * "\nMaterial parameters:\n"
    msg *= @sprintf "  - Number of material points [-]:      %30g\n" pc.n_points
    msg *= describe_mat(mat)
    print(msg)
    if es.exportflag
        open(es.logfile, "a") do io
            return write(io, msg)
        end
    end
    return nothing
end

function describe_mat(mat::AbstractPDMaterial)
    msg = @sprintf "  - Material model:                     %30s\n" typeof(mat)
    msg *= @sprintf "  - Horizon δ [m]:                      %30g\n" mat.δ
    msg *= @sprintf "  - Density ρ [kg/m³]:                  %30g\n" mat.rho
    msg *= @sprintf "  - Young's modulus E [N/m²]:           %30g\n" mat.E
    msg *= @sprintf "  - Poisson ratio ν [-]:                %30g\n" mat.nu
    msg *= @sprintf "  - Critical bond stretch εc[-]:        %30g\n" mat.εc
    return msg
end

function log_describe_interactions(part::AbstractPDBody, es::ExportSettings)
    msg = "Interactions:\n"
    msg *= describe_interactions(part)
    msg *= @sprintf(
        "Total memory used by body [MB]:         %30g\n", Base.summarysize(part) * 1e-6
    )
    print(msg)
    if es.exportflag
        open(es.logfile, "a") do io
            return write(io, msg)
        end
    end
    return nothing
end

function describe_interactions(part::AbstractPDBody)
    msg = @sprintf "  - Number of bonds [-]:                %30d\n" part.n_bonds
    return msg
end

function log_simsetup(sim::AbstractPDAnalysis)
    msg = ""
    if sim.es.exportflag
        msg *= "Export setup:\n"
        msg *= @sprintf "  - Export frequency:                   %30g\n" sim.es.exportfreq
        if length(sim.es.resultfile_prefix) <= 48
            msg *= @sprintf "  - Export file name: %48s\n" sim.es.resultfile_prefix
        else
            msg *= @sprintf "  - Export file name:\n    %s\n" sim.es.resultfile_prefix
        end
    end
    msg *= "Time discretization:\n"
    msg *= @sprintf "  - Time step Δt [s]:                   %30g\n" sim.td.Δt
    msg *= @sprintf "  - Number of time steps [-]:           %30g\n" sim.td.n_timesteps
    msg *= @sprintf "  - Simulation time horizon [s]:        %30g\n" sim.td.n_timesteps *
        sim.td.Δt
    print(msg)
    if sim.es.exportflag
        open(sim.es.logfile, "a") do io
            return write(io, msg)
        end
    end
    return nothing
end

function log_closesim(part::AbstractPDBody, timingsummary::NamedTuple, es::ExportSettings)
    msg = @sprintf(
        "✓ Simulation completed after %g seconds\nResults:\n", timingsummary.time
    )
    msg *= @sprintf(
        "  - Max. abs. x-displacement [m]:       %30g\n",
        abs(maximum(part.displacement[1, :]))
    )
    msg *= @sprintf(
        "  - Max. abs. y-displacement [m]:       %30g\n",
        abs(maximum(part.displacement[2, :]))
    )
    msg *= @sprintf(
        "  - Max. abs. z-displacement [m]:       %30g\n",
        abs(maximum(part.displacement[3, :]))
    )
    msg *= @sprintf("  - Max. damage [-]:                    %30g\n", maximum(part.damage))
    print(msg)
    if es.exportflag
        open(es.logfile, "a") do io
            return write(io, msg)
        end
    end
    return nothing
end

function define_precrack!(body::AbstractPDBody, precrack::PreCrack)
    if body.unique_bonds # TODO: material dependent choicing with multiple dispatch
        @threads for _ in 1:(body.n_threads)
            tid = threadid()
            (@view body.n_active_family_members[:, tid]) .= 0
            for current_one_ni in body.owned_bonds[tid]
                i, j, _, _ = body.bond_data[current_one_ni]
                i_is_in_set_a = in(i, precrack.point_id_set_a)
                i_is_in_set_b = in(i, precrack.point_id_set_b)
                j_is_in_set_a = in(j, precrack.point_id_set_a)
                j_is_in_set_b = in(j, precrack.point_id_set_b)
                if i_is_in_set_a && j_is_in_set_b || i_is_in_set_b && j_is_in_set_a
                    body.bond_failure[current_one_ni] = 0
                end
                body.n_active_family_members[i, tid] += body.bond_failure[current_one_ni]
                body.n_active_family_members[j, tid] += body.bond_failure[current_one_ni]
            end
        end
    end
    return nothing
end

function calc_stable_timestep(body::AbstractPDBody, rho::Float64, K::Float64, δ::Float64)
    _Δt = zeros(Float64, body.n_threads)
    @threads for _ in 1:(body.n_threads)
        tid = threadid()
        timesteps = zeros(Float64, body.n_points)
        dtsum = zeros(Float64, (body.n_points, body.n_threads))
        for current_one_ni in body.owned_bonds[tid]
            (a, i, L, _) = body.bond_data[current_one_ni]
            dtsum[a, tid] += body.volumes[i] * 1 / L * 18 * K / (π * δ^4)
        end
        for a in body.owned_points[tid]
            dtsum[a, 1] = sum(@view dtsum[a, :])
            timesteps[a] = √(2 * rho / dtsum[a, 1])
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
function calc_stable_user_timestep(pc::PointCloud, mat::AbstractPDMaterial, Sf::Float64=0.7)
    owned_points = defaultdist(pc.n_points, nthreads())
    bond_data, _ = find_bonds(pc, mat.δ, owned_points)
    owned_bonds = defaultdist(length(bond_data), nthreads())
    _Δt = zeros(Float64, nthreads())
    @threads for _ in 1:nthreads()
        tid = threadid()
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

function calc_damage!(body::AbstractPDBody)
    @threads for _ in 1:(body.n_threads)
        for a in body.owned_points[threadid()]
            body.n_active_family_members[a, 1] = sum(
                @view body.n_active_family_members[a, :]
            )
            body.damage[a] =
                1 - body.n_active_family_members[a, 1] / body.n_family_members[a]
        end
    end
    return nothing
end

function apply_bcs!(body::AbstractPDBody, bcs::Vector{<:AbstractBC}, t::Float64)
    @sync for bc in bcs
        Threads.@spawn apply_boundarycondition!(body, bc, t)
    end
    return nothing
end

function apply_boundarycondition!(body::AbstractPDBody, bc::VelocityBC, t::Float64)
    value = bc.fun(t)
    dim = bc.dim
    @simd for i in bc.point_id_set
        body.velocity_half[dim, i] = value
    end
    return nothing
end

function apply_boundarycondition!(body::AbstractPDBody, bc::ForceDensityBC, t::Float64)
    value = bc.fun(t)
    dim = bc.dim
    @simd for i in bc.point_id_set
        body.b_ext[dim, i] = value
    end
    return nothing
end

function apply_boundarycondition!(body::AbstractPDBody, bc::PosDepVelBC, t::Float64)
    dim = bc.dim
    @simd for i in bc.point_id_set
        body.velocity_half[dim, i] = bc.fun(
            body.position[1, i], body.position[2, i], body.position[3, i], t
        )
    end
    return nothing
end

function apply_ics!(body::AbstractPDBody, ics::Vector{<:AbstractIC})
    @sync for ic in ics
        Threads.@spawn apply_initialcondition!(body, ic)
    end
    return nothing
end

function apply_initialcondition!(body::AbstractPDBody, ic::VelocityIC)
    dim = ic.dim
    @simd for i in ic.point_id_set
        body.velocity[dim, i] = ic.val
    end
    return nothing
end

function export_jld2(body::AbstractPDBody, expfile::String, timestep::Int, time::Float64)
    filename = expfile * "_t" * string(timestep) * ".jld2"
    save(
        filename,
        "position",
        body.position,
        "damage",
        body.damage,
        "b_int",
        (@view body.b_int[:, :, 1]),
        "displacement",
        body.displacement,
        "time",
        time,
    )
    return nothing
end

function export_vtk(body::AbstractPDBody, expfile::String, timestep::Int, time::Float64)
    filename = expfile * "_t" * string(timestep)
    cells = [MeshCell(VTKCellTypes.VTK_VERTEX, (j,)) for j in 1:(body.n_points)]
    vtkfile = vtk_grid(filename, body.position, cells)
    vtkfile["Damage", VTKPointData()] = body.damage
    vtkfile["Force density", VTKPointData()] = @views body.b_int[:, :, 1]
    vtkfile["Displacement", VTKPointData()] = body.displacement
    vtkfile["Velocity", VTKPointData()] = body.velocity
    vtkfile["Time", VTKFieldData()] = time
    vtk_save(vtkfile)
    return nothing
end

function export_results(body::AbstractPDBody, expfile::String, timestep::Int, time::Float64)
    export_vtk(body, expfile, timestep, time)
    export_jld2(body, expfile, timestep, time)
    return nothing
end

function velocity_verlet!(body::AbstractPDBody, sim::PDSingleBodyAnalysis)
    p = Progress(
        sim.td.n_timesteps; dt=1, desc="Time integration... ", barlen=30, color=:normal
    )
    Δt½ = 0.5 * sim.td.Δt
    for t in 1:(sim.td.n_timesteps)
        time = t * sim.td.Δt
        update_velhalf!(body, Δt½)
        apply_bcs!(body, sim.bcs, time)
        update_disp_and_position!(body, sim.td.Δt)
        compute_forcedensity!(body, sim.mat)
        calc_damage!(body)
        compute_equation_of_motion!(body, Δt½, sim.mat.rho)
        if mod(t, sim.es.exportfreq) == 0
            export_results(body, sim.es.resultfile_prefix, t, time)
        end
        next!(p)
    end
    finish!(p)
    return nothing
end

function dynamic_relaxation_finite_difference!(
    body::AbstractPDBody, sim::PDSingleBodyAnalysis
)
    damping_matrix = fill(18 * sim.mat.K * sim.td.Δt^2 / sim.mat.δ, (3, body.n_points))
    velocity_half_old = zeros(Float64, (3, body.n_points))
    b_int_old = zeros(Float64, (3, body.n_points))
    p = Progress(
        sim.td.n_timesteps; dt=1, desc="Time integration... ", barlen=30, color=:normal
    )
    for t in 1:(sim.td.n_timesteps)
        time = t * sim.td.Δt
        apply_bcs!(body, sim.bcs, time)
        compute_forcedensity!(body, sim.mat)
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
            export_results(body, sim.es.resultfile_prefix, t, time)
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
    for i in 1:(body.n_points)
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
    @threads for i in 1:(body.n_points)
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
    @threads for i in 1:(body.n_points)
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
    @threads for i in 1:(body.n_points)
        body.velocity_half[1, i] = body.velocity[1, i] + body.acceleration[1, i] * Δt½
        body.velocity_half[2, i] = body.velocity[2, i] + body.acceleration[2, i] * Δt½
        body.velocity_half[3, i] = body.velocity[3, i] + body.acceleration[3, i] * Δt½
    end
    return nothing
end

function update_disp_and_position!(body::AbstractPDBody, Δt::Float64)
    @threads for i in 1:(body.n_points)
        body.displacement[1, i] += body.velocity_half[1, i] * Δt
        body.displacement[2, i] += body.velocity_half[2, i] * Δt
        body.displacement[3, i] += body.velocity_half[3, i] * Δt
        body.position[1, i] += body.velocity_half[1, i] * Δt
        body.position[2, i] += body.velocity_half[2, i] * Δt
        body.position[3, i] += body.velocity_half[3, i] * Δt
    end
    return nothing
end

function compute_equation_of_motion!(body::AbstractPDBody, Δt½::Float64, rho::Float64)
    @threads for i in 1:(body.n_points)
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
