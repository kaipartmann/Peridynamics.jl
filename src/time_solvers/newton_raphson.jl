"""
    NewtonRaphson(; kwargs...)

$(internal_api_warning())
$(experimental_api_warning())

An implicit time integration solver for quasi-static peridynamic simulations using the
Jacobian-Free Newton-Krylov (JFNK) method.

The solver uses a matrix-free approach where the Jacobian-vector product is computed
via central finite differences: `J*v ≈ (F(u + ε*v) - F(u - ε*v)) / (2ε)`. This eliminates
the need to store the full Jacobian matrix, reducing memory from O(n²) to O(n), enabling
simulations with large numbers of points.

GMRES is used as the Krylov solver for the resulting linear system at each iteration.

# Keywords:
- `time::Real`: Total simulation time. Cannot be used together with `steps`.
- `steps::Int`: Number of time steps. Cannot be used together with `time`.
- `stepsize::Real`: Manual time step size (default: 1.0).
- `maxiter::Int`: Maximum number of iterations per step (default: 100).
- `tol::Real`: Tolerance for convergence (default: 1e-4).
- `perturbation::Real`: Perturbation size for finite difference approximation of the
    Jacobian-vector product (default: 1e-7 times the point spacing).
- `gmres_maxiter::Int`: Maximum number of iterations for GMRES (default: min(200, n_dof)).
- `gmres_reltol::Real`: Relative tolerance for GMRES (default: 1e-4).
- `gmres_abstol::Real`: Absolute tolerance for GMRES (default: 1e-8).

# Example:
```julia
solver = NewtonRaphson(steps=100, stepsize=1e-4, maxiter=50, tol=1e-6)
```
"""
mutable struct NewtonRaphson <: AbstractTimeSolver
    end_time::Float64
    n_steps::Int
    Δt::Float64
    maxiter::Int
    tol::Float64
    perturbation::Float64
    gmres_maxiter::Int
    gmres_reltol::Float64
    gmres_abstol::Float64

    function NewtonRaphson(; time::Real=-1, steps::Int=-1, stepsize::Real=1.0,
                            maxiter::Int=100, tol::Real=1e-4, perturbation::Real=-1,
                            gmres_maxiter::Int=-1, gmres_reltol::Real=1e-4,
                            gmres_abstol::Real=1e-8)
        if time > 0 && steps > 0
            msg = "specify either time or number of steps, not both!"
            throw(ArgumentError(msg))
        elseif time < 0 && steps < 0
            msg = "specify either time or number of steps!"
            throw(ArgumentError(msg))
        end
        if stepsize ≤ 0
            msg = "stepsize has to be larger than zero: stepsize > 0"
            throw(ArgumentError(msg))
        end

        # Calculate derived values
        if time > 0
            n_steps = Int(ceil(time / stepsize))
            end_time = time
        else
            n_steps = steps
            end_time = steps * stepsize
        end

        if maxiter ≤ 0
            msg = "maxiter has to be larger than zero: maxiter > 0"
            throw(ArgumentError(msg))
        end
        if tol ≤ 0
            msg = "tolerance has to be larger than zero: tol > 0"
            throw(ArgumentError(msg))
        end
        if gmres_reltol ≤ 0
            msg = "gmres_reltol has to be larger than zero: gmres_reltol > 0"
            throw(ArgumentError(msg))
        end
        if gmres_abstol ≤ 0
            msg = "gmres_abstol has to be larger than zero: gmres_abstol > 0"
            throw(ArgumentError(msg))
        end
        new(end_time, n_steps, stepsize, maxiter, tol, perturbation, gmres_maxiter,
            gmres_reltol, gmres_abstol)
    end
end

function Base.show(io::IO, @nospecialize(nr::NewtonRaphson))
    print(io, typeof(nr))
    fields = Vector{Symbol}()
    for field in fieldnames(typeof(nr))
        value = getfield(nr, field)
        if value > 0
            push!(fields, field)
        end
    end
    print(io, msg_fields_in_brackets(nr, Tuple(fields)))
    return nothing
end

function Base.show(io::IO, ::MIME"text/plain", @nospecialize(nr::NewtonRaphson))
    if get(io, :compact, false)
        show(io, nr)
    else
        println(io, typeof(nr), ":")
        fields = Vector{Symbol}()
        for field in fieldnames(typeof(nr))
            value = getfield(nr, field)
            if value > 0
                push!(fields, field)
            end
        end
        print(io, msg_fields(nr, Tuple(fields)))
    end
    return nothing
end

# Always use only one thread!
function chop_body_threads(body::AbstractBody, solver::NewtonRaphson,
                           point_decomp::PointDecomposition, param_spec::AbstractParamSpec)
    if point_decomp.n_chunks > 1
        msg = "NewtonRaphson solver currently not implemented for multithreading usage!\n"
        msg *= "For now, it will force `n_chunks = 1`\n"
        println() # to not interfere with logging
        @warn msg
    end
    _point_decomp = PointDecomposition(body, 1) # force a single chunk for now!
    ChunkType = body_chunk_type(body, solver, param_spec)
    chunks = _chop_body_threads(ChunkType, body, solver, _point_decomp, param_spec)
    return chunks
end

function init_time_solver!(nr::NewtonRaphson, dh::ThreadsBodyDataHandler)
    newton_raphson_check(nr)

    # Check that only a single chunk is created
    if dh.n_chunks != 1
        msg = "NewtonRaphson solver currently only implemented for single chunk!\n"
        msg *= "Number of chunks in data handler: $(dh.n_chunks)\n"
        throw(ArgumentError(msg))
    end
    chunk = dh.chunks[1] # extract for the checks

    # Boundary condition checks
    (; condhandler) = chunk
    n_sdbcs = length(condhandler.single_dim_bcs)
    n_pdsdbcs = length(condhandler.posdep_single_dim_bcs)
    n_dbcs = length(condhandler.data_bcs)
    if n_sdbcs + n_pdsdbcs + n_dbcs > 0
        msg = "NewtonRaphson solver does not support all boundary condition types!\n"
        msg *= "Specify only position-dependent boundary conditions!"
        throw(ArgumentError(msg))
    end
    if length(condhandler.pos_single_dim_bcs) == 0
        msg = "NewtonRaphson solver requires at least one boundary condition!\n"
        msg *= "Specify a boundary conditions with `displacement_bc!`, or \n"
        msg *= "`forcedensity_bc!` with a function signature `p -> ...`!"
        throw(ArgumentError(msg))
    end

    # Material limitations
    incompatible_mats = (:CRMaterial, :RKCRMaterial)
    if nameof(typeof(chunk.mat)) ∈ incompatible_mats
        msg = "$(nameof(typeof(chunk.mat))) not compatible with Newton-Raphson!\n"
        msg *= "The material uses stress rotation designed for dynamic simulations.\n"
        msg *= "For quasi-static Newton-Raphson simulations, consider using the material"
        msg *= " without stress rotation instead.\n"
        throw(ArgumentError(msg))
    end

    # No damage allowed in whole model
    (; mat, system, paramsetup) = chunk
    for i in each_point_idx(system)
        params = get_params(paramsetup, i)
        if has_fracture(mat, params)
            msg = "NewtonRaphson solver currently does not support damage or fracture!\n"
            msg *= "Fracture parameters found at point index $(i)!\n"
            throw(ArgumentError(msg))
        end
    end

    # Calculate perturbation size based on system characteristics
    if nr.perturbation ≤ 0
        # Use characteristic length scale based on volume
        Δx = (minimum_volume(dh))^(1/3) # minimum grid spacing
        nr.perturbation = Δx * 1e-7 # perturbation size
    end

    # Set GMRES max iterations based on system size
    if nr.gmres_maxiter ≤ 0
        nr.gmres_maxiter = min(200, get_n_dof(system))
    end

    return nothing
end

# Overload for MultibodySetups to throw an error
function init_time_solver!(::NewtonRaphson, ::AbstractDataHandler)
    msg = "NewtonRaphson solver only implemented for single body multithreading!\n"
    return throw(ArgumentError(msg))
end

function minimum_volume(dh::ThreadsBodyDataHandler)
    # min_vols = fill(Inf, dh.n_chunks)
    # @threads :static for chunk_id in eachindex(dh.chunks)
    #     chunk = dh.chunks[chunk_id]
    #     min_vols[chunk_id] = minimum(chunk.system.volume)
    # end
    # return minimum(min_vols)
    return minimum(dh.chunks[1].system.volume)
end

# TODO: Add MPI version later
# function minimum_volume(dh::MPIBodyDataHandler)
#     local_min_vol = minimum_volume(dh.chunk.system.volume)
#     global_min_vol = MPI.Allreduce(local_min_vol, MPI.MIN, mpi_comm())
#     return global_min_vol
# end

function newton_raphson_check(nr::NewtonRaphson)
    if nr.end_time < 0
        error("`end_time` of NewtonRaphson smaller than zero!\n")
    end
    if nr.n_steps < 0
        error("`n_steps` of NewtonRaphson smaller than zero!\n")
    end
    if nr.Δt < 0
        error("`Δt` of NewtonRaphson smaller than zero!\n")
    end
    if nr.maxiter < 0
        error("`maxiter` of NewtonRaphson smaller than zero!\n")
    end
    if nr.tol < 0
        error("`tol` of NewtonRaphson smaller than zero!\n")
    end
    return nothing
end

function solve!(dh::AbstractDataHandler, nr::NewtonRaphson, options::AbstractJobOptions)
    export_reference_results(dh, options)

    # Log solver info
    log_it(options, "NEWTON-RAPHSON ITERATIONS")
    print_log(stdout, "\n") # necessary for correct terminal output

    # Time stepping loop
    for n in 1:nr.n_steps
        newton_raphson_step!(dh, nr, options, n)
    end

    add_to_logfile(options, "\n") # necessary for correct logfile output

    return dh
end

# Newton-Raphson step function, for now only for a single chunk!
function newton_raphson_step!(dh::ThreadsBodyDataHandler,
                              nr::NewtonRaphson, options::AbstractJobOptions, n::Int)
    chunk = dh.chunks[1]
    (; n_steps, Δt, maxiter, tol) = nr
    (; system, storage, condhandler) = chunk
    (; free_dofs, constrained_dofs) = condhandler

    t = n * Δt
    β = n / n_steps # load step factor

    # Apply boundary conditions and update position
    apply_incr_boundary_conditions!(chunk, β)
    update_position!(storage, system, constrained_dofs)

    for iter in 1:maxiter
        # Calculate internal forces
        calc_force_density!(dh, t, Δt)
        calc_residual!(chunk)

        # Check convergence using combined residual from all chunks
        r = get_residual_norm(dh)
        msg = @sprintf("\r  step: %8d | iter: %8d | r: %16g", n, iter, r)
        log_it(options, msg)
        if r < tol
            print_log(stdout, " ✔\n")
            add_to_logfile(options, "\n  " * "-"^42 * "> converged ✔")
            break
        elseif iter == maxiter
            print_log(stdout, " ❌\n\n")
            add_to_logfile(options, "\n  " * "-"^35 * "> did not converge ❌")
            msg = "Newton-Raphson solver did not converge after max iterations!\n"
            msg *= "Consider increasing `maxiter` or `tol`.\n"
            throw(ErrorException(msg))
        end

        # Solve linear system using Jacobian-free Newton-Krylov method
        solve_linear_system!(chunk, nr, t, Δt)
    end

    @threads for chunk_id in eachindex(dh.chunks)
        export_results(dh, options, chunk_id, n, t)
    end

    return nothing
end

function update_position!(storage::AbstractStorage, system::AbstractSystem,
                          dofs::AbstractVector{Int})
    (; position, displacement) = storage
    for idx in dofs
        @inbounds position[idx] = system.position[idx] + displacement[idx]
    end
    return nothing
end

function calc_residual!(chunk::AbstractBodyChunk)
    (; system, storage, condhandler) = chunk
    (; residual, b_int, b_ext) = storage
    (; volume) = system

    # Convert 2D force matrices to 1D residual vector
    for (dof, _, i) in each_dof_idx(system)
        residual[dof] = (b_int[dof] + b_ext[dof]) * volume[i] # scaling with volume
    end

    # Set residual to zero for constrained DOFs
    for dof in condhandler.constrained_dofs
        residual[dof] = 0.0
    end

    return nothing
end

function get_residual_norm(dh::ThreadsBodyDataHandler)
    return get_residual_norm(dh.chunks[1])
end

function get_residual_norm(chunk::AbstractBodyChunk)
    return norm(chunk.storage.residual)
end

"""
    jacobian_vector_product!(Jv, chunk, solver, v, t, Δt)

Compute the Jacobian-vector product `J * v` using central finite differences.
The result is stored in `Jv`.

This is the core of the Jacobian-Free Newton-Krylov (JFNK) method.
Uses central difference approximation for O(ε²) accuracy:
    J * v ≈ (F(u + ε*v) - F(u - ε*v)) / (2ε)
"""
function jacobian_vector_product!(Jv::AbstractVector, chunk::AbstractBodyChunk,
                                  solver::NewtonRaphson, v::AbstractVector, t, Δt)
    (; storage, system, mat, paramsetup, condhandler) = chunk
    (; b_int, position, displacement) = storage
    (; b_int_copy, displacement_copy, temp_force_a) = storage
    (; volume) = system
    (; constrained_dofs) = condhandler
    ε = solver.perturbation

    # Store original displacement and force density
    displacement_copy .= displacement
    b_int_copy .= b_int

    # Apply positive perturbation: u + ε*v
    for (dof, _, _) in each_loc_dof_idx(system)
        @inbounds displacement[dof] = displacement_copy[dof] + ε * v[dof]
        @inbounds position[dof] = system.position[dof] + displacement[dof]
    end

    # Calculate forces at positively perturbed state
    b_int .= 0.0
    calc_force_density!(storage, system, mat, paramsetup, t, Δt)

    # Compute F(u + ε*v) * volume and store in temp_force_a
    for (dof, _, i) in each_dof_idx(system)
        @inbounds temp_force_a[dof] = b_int[dof] * volume[i]
    end

    # Apply negative perturbation: u - ε*v
    for (dof, _, _) in each_loc_dof_idx(system)
        @inbounds displacement[dof] = displacement_copy[dof] - ε * v[dof]
        @inbounds position[dof] = system.position[dof] + displacement[dof]
    end

    # Calculate forces at negatively perturbed state
    b_int .= 0.0
    calc_force_density!(storage, system, mat, paramsetup, t, Δt)

    # Compute Jacobian-vector product: Jv = (F(u + εv) - F(u - εv)) / (2ε)
    for (dof, _, i) in each_dof_idx(system)
        @inbounds Jv[dof] = (temp_force_a[dof] - b_int[dof] * volume[i]) / (2ε)
    end

    # Handle constrained DOFs: set Jv[dof] = v[dof] (identity row)
    for dof in constrained_dofs
        @inbounds Jv[dof] = v[dof]
    end

    # Restore original displacement, position, and force density
    for (dof, _, _) in each_loc_dof_idx(system)
        @inbounds displacement[dof] = displacement_copy[dof]
        @inbounds position[dof] = system.position[dof] + displacement[dof]
    end
    b_int .= b_int_copy

    return nothing
end

function solve_linear_system!(chunk::AbstractBodyChunk, solver::NewtonRaphson, t, Δt)
    (; storage, condhandler) = chunk
    (; free_dofs) = condhandler
    (; position, displacement, residual, Δu) = storage
    (; gmres_maxiter, gmres_reltol, gmres_abstol) = solver

    Δu .= 0.0

    # Create Jacobian-free linear operator using a closure
    n_dof = get_n_dof(chunk.system)
    J_linmap = LinearMap(n_dof; issymmetric=false, ismutating=true) do Jv, v
        jacobian_vector_product!(Jv, chunk, solver, v, t, Δt)
    end

    # Solve using GMRES with the matrix-free operator
    # Use larger restart value to improve convergence for ill-conditioned systems
    maxiter, reltol, abstol = gmres_maxiter, gmres_reltol, gmres_abstol
    restart = min(50, n_dof)  # Restart after 50 iterations
    gmres!(Δu, J_linmap, -residual; maxiter, reltol, abstol, restart)

    for dof in free_dofs
        @inbounds du = Δu[dof]
        @inbounds displacement[dof] += du
        @inbounds position[dof] += du
    end

    return nothing
end


# Required interface functions
function init_field_solver(::NewtonRaphson, system::AbstractSystem, ::Val{:position})
    return copy(system.position)
end

function init_field_solver(::NewtonRaphson, system::AbstractSystem, ::Val{:displacement})
    return zeros(3, get_n_loc_points(system))
end

function init_field_solver(::NewtonRaphson, system::AbstractSystem, ::Val{:b_int})
    return zeros(3, get_n_points(system))
end

function init_field_solver(::NewtonRaphson, system::AbstractSystem, ::Val{:b_ext})
    return zeros(3, get_n_points(system))
end

# Newton-Raphson specific fields
function init_field_solver(::NewtonRaphson, system::AbstractSystem,
                           ::Val{:displacement_copy})
    return zeros(3, get_n_loc_points(system))
end
function init_field_solver(::AbstractTimeSolver, system::AbstractSystem,
                           ::Val{:displacement_copy})
    return Array{Float64,2}(undef, 0, 0)
end

function init_field_solver(::NewtonRaphson, system::AbstractSystem,
                           ::Val{:b_int_copy})
    return zeros(3, get_n_points(system))
end
function init_field_solver(::AbstractTimeSolver, system::AbstractSystem,
                           ::Val{:b_int_copy})
    return Array{Float64,2}(undef, 0, 0)
end

function init_field_solver(::NewtonRaphson, system::AbstractSystem, ::Val{:residual})
    return zeros(get_n_dof(system))
end
function init_field_solver(::AbstractTimeSolver, system::AbstractSystem, ::Val{:residual})
    return Vector{Float64}()
end

# Note: jacobian field removed - using Jacobian-free Newton-Krylov method
function init_field_solver(::AbstractTimeSolver, system::AbstractSystem,
                           ::Val{:jacobian})
    return Array{Float64,2}(undef, 0, 0)
end

function init_field_solver(::NewtonRaphson, system::AbstractSystem, ::Val{:temp_force_a})
    return zeros(get_n_dof(system))
end
function init_field_solver(::AbstractTimeSolver, system::AbstractSystem,
                           ::Val{:temp_force_a})
    return Array{Float64,1}(undef, 0)
end

# Note: temp_force_b field removed - using forward difference in JFNK
function init_field_solver(::AbstractTimeSolver, system::AbstractSystem,
                           ::Val{:temp_force_b})
    return Array{Float64,1}(undef, 0)
end

function init_field_solver(::NewtonRaphson, system::AbstractSystem, ::Val{:Δu})
    return zeros(get_n_dof(system))
end
function init_field_solver(::AbstractTimeSolver, system::AbstractSystem, ::Val{:Δu})
    return Array{Float64,1}(undef, 0)
end

# Note: affected_points field removed - not needed for matrix-free approach
function init_field_solver(::AbstractTimeSolver, system::AbstractSystem,
                           ::Val{:affected_points})
    return Vector{Vector{Int}}()
end

function req_point_data_fields_timesolver(::Type{<:NewtonRaphson})
    # Jacobian-free Newton-Krylov method: no jacobian, temp_force_b, or affected_points
    fields = (:position, :displacement, :b_int, :b_ext,
              :residual, :displacement_copy, :temp_force_a, :b_int_copy, :Δu)
    return fields
end

function req_bond_data_fields_timesolver(::Type{<:NewtonRaphson})
    return ()
end

function req_data_fields_timesolver(::Type{<:NewtonRaphson})
    return ()
end

function log_timesolver(options::AbstractJobOptions, nr::NewtonRaphson)
    msg = "JACOBIAN-FREE NEWTON-KRYLOV SOLVER\n"
    msg *= msg_qty("number of time steps", nr.n_steps)
    msg *= msg_qty("time step size", nr.Δt)
    msg *= msg_qty("simulation time", nr.end_time)
    msg *= msg_qty("max iterations per step", nr.maxiter)
    msg *= msg_qty("convergence tolerance", nr.tol)
    msg *= msg_qty("perturbation size", nr.perturbation)
    msg *= msg_qty("GMRES max iterations", nr.gmres_maxiter)
    msg *= msg_qty("GMRES relative tolerance", nr.gmres_reltol)
    msg *= msg_qty("GMRES absolute tolerance", nr.gmres_abstol)
    log_it(options, msg)
    return nothing
end

# If this is uncommented, the solver is registered and all storages must support it!
# The internal storages in this package, already do, but the user-defined ones may not.
# Therefore, for now this is commented out.
# register_solver!(NewtonRaphson)
