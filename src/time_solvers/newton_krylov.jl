"""
    NewtonKrylov(; kwargs...)

An implicit time integration solver for quasi-static peridynamic simulations using the
Jacobian-Free Newton-Krylov (JFNK) method.

The solver uses a matrix-free approach where the Jacobian-vector product is computed
via central finite differences: `J*v ≈ (F(u + ε*v) - F(u - ε*v)) / (2ε)`. This eliminates
the need to store the full Jacobian matrix, reducing memory from O(n²) to O(n), enabling
simulations with large numbers of points.

The perturbation `ε` is computed adaptively using the standard JFNK heuristic for central
differences: `ε = scale * cbrt(ε_mach) * (||u|| + 1) / ||v||`, where `ε_mach` is machine
epsilon. This balances truncation error O(ε²) with roundoff error O(ε_mach/ε).

GMRES is used as the Krylov solver for the resulting linear system at each iteration,
with an optional diagonal preconditioner computed by probing the Jacobian.

# Keywords:
- `time::Real`: Total simulation time. Cannot be used together with `steps`.
- `steps::Int`: Number of time steps. Cannot be used together with `time`.
- `stepsize::Real`: Manual time step size (default: 1.0).
- `maxiter::Int`: Maximum number of iterations per step (default: 100).
- `tol::Real`: Tolerance for convergence (default: 1e-4).
- `perturbation_scale::Real`: Scale factor for the adaptive perturbation heuristic
    (default: 1.0). Increase if finite difference approximation is noisy.
- `gmres_maxiter::Int`: Maximum number of iterations for GMRES (default: 200).
- `gmres_reltol::Real`: Relative tolerance for GMRES (default: 1e-4).
- `gmres_abstol::Real`: Absolute tolerance for GMRES (default: 1e-8).
- `gmres_restart::Int`: Number of iterations between restarts for GMRES (default: 50).
- `preconditioner::Bool`: Whether to use a diagonal right preconditioner (default: true).
    The preconditioner is based on material stiffness scaling, computed once at initialization
    without any additional force evaluations.
- `damping::Real`: Damping coefficient for load-step regularization (default: 0.0).
    When damping > 0, adds a regularization term `damping * α * K * u` to the residual,
    where α linearly decays from 1 at step 1 to 0 at the final step. This helps stabilize
    convergence for materials with zero-energy modes (e.g., correspondence formulations).
    The damping is applied to the Jacobian as well, adding `damping * α * K` to the diagonal.
    Typical values are 0.01 to 0.1.

# Example:
```julia
solver = NewtonKrylov(steps=100, stepsize=1e-4, maxiter=50, tol=1e-6)
# With damping for correspondence materials:
solver = NewtonKrylov(steps=100, stepsize=1e-4, damping=0.05)
```
"""
mutable struct NewtonKrylov <: AbstractTimeSolver
    end_time::Float64
    n_steps::Int
    Δt::Float64
    maxiter::Int
    tol::Float64
    perturbation_scale::Float64
    gmres_maxiter::Int
    gmres_reltol::Float64
    gmres_abstol::Float64
    gmres_restart::Int
    use_preconditioner::Bool
    damping::Float64
    global_residual::Vector{Float64}
    global_Δu::Vector{Float64}
    global_precond::Vector{Float64}
    mpi_dof_counts::Vector{Int32}
    mpi_displs::Vector{Int32}
    # Pre-allocated buffers for norm computation (avoid allocations in hot path)
    u_norm_sq_buf::Vector{Float64}
    v_norm_sq_buf::Vector{Float64}
    # Damping stiffness scaling per DOF (computed at initialization)
    damping_diag::Vector{Float64}

    function NewtonKrylov(; time::Real=-1, steps::Int=-1, stepsize::Real=1.0,
                            maxiter::Int=100, tol::Real=1e-4, perturbation_scale::Real=1.0,
                            gmres_maxiter::Int=200, gmres_reltol::Real=1e-4,
                            gmres_abstol::Real=1e-8, gmres_restart::Int=50,
                            preconditioner::Bool=true, damping::Real=0.0)
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

        # Validate parameters
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
        if gmres_maxiter ≤ 0
            msg = "gmres_maxiter has to be larger than zero: gmres_maxiter > 0"
            throw(ArgumentError(msg))
        end
        if gmres_restart ≤ 0
            msg = "gmres_restart has to be larger than zero: gmres_restart > 0"
            throw(ArgumentError(msg))
        end

        # Validate perturbation_scale
        if perturbation_scale ≤ 0
            msg = "perturbation_scale has to be larger than zero: perturbation_scale > 0"
            throw(ArgumentError(msg))
        end

        # Validate damping
        if damping < 0
            msg = "damping must be non-negative: damping >= 0"
            throw(ArgumentError(msg))
        end

        # Initialize global GMRES buffers as empty (will be resized on first use)
        # MPI buffers also initialized empty
        global_residual = Vector{Float64}()
        global_Δu = Vector{Float64}()
        global_precond = Vector{Float64}()
        # MPI-specific cached data (lazily initialized)
        mpi_dof_counts = Vector{Int32}()
        mpi_displs = Vector{Int32}()
        # Norm computation buffers (lazily initialized)
        u_norm_sq_buf = Vector{Float64}()
        v_norm_sq_buf = Vector{Float64}()
        # Damping diagonal (lazily initialized)
        damping_diag = Vector{Float64}()

        new(end_time, n_steps, stepsize, maxiter, tol, perturbation_scale, gmres_maxiter,
            gmres_reltol, gmres_abstol, gmres_restart, preconditioner, damping,
            global_residual, global_Δu, global_precond, mpi_dof_counts, mpi_displs,
            u_norm_sq_buf, v_norm_sq_buf, damping_diag)
    end
end

function Base.show(io::IO, @nospecialize(nr::NewtonKrylov))
    print(io, typeof(nr))
    fields = Vector{Symbol}()
    excluded = (:global_residual, :global_Δu, :global_precond, :mpi_dof_counts, :mpi_displs,
                :u_norm_sq_buf, :v_norm_sq_buf, :damping_diag)
    for field in fieldnames(typeof(nr))
        field in excluded && continue
        value = getfield(nr, field)
        if value isa Number && value > 0
            push!(fields, field)
        elseif value isa Bool
            push!(fields, field)
        end
    end
    print(io, msg_fields_in_brackets(nr, Tuple(fields)))
    return nothing
end

function Base.show(io::IO, ::MIME"text/plain", @nospecialize(nr::NewtonKrylov))
    if get(io, :compact, false)
        show(io, nr)
    else
        println(io, typeof(nr), ":")
        fields = Vector{Symbol}()
        excluded = (:global_residual, :global_Δu, :global_precond, :mpi_dof_counts,
                    :mpi_displs, :u_norm_sq_buf, :v_norm_sq_buf, :damping_diag)
        for field in fieldnames(typeof(nr))
            field in excluded && continue
            value = getfield(nr, field)
            if value isa Number && value > 0
                push!(fields, field)
            elseif value isa Bool
                push!(fields, field)
            end
        end
        print(io, msg_fields(nr, Tuple(fields)))
    end
    return nothing
end

function init_time_solver!(nr::NewtonKrylov, dh::ThreadsBodyDataHandler)
    newton_krylov_check(nr)
    @threads :static for chunk in dh.chunks
        validate_chunk_for_newton_krylov!(chunk)
    end
    init_newton_buffers!(nr, dh)
    return nothing
end

# MPI support
function init_time_solver!(nr::NewtonKrylov, dh::MPIBodyDataHandler)
    newton_krylov_check(nr)
    validate_chunk_for_newton_krylov!(dh.chunk)
    init_newton_buffers!(nr, dh)
    return nothing
end

function init_newton_buffers!(nr::NewtonKrylov, dh::ThreadsBodyDataHandler)
    total_loc_dof = sum(get_n_loc_dof(c.system) for c in dh.chunks)
    resize!(nr.global_residual, total_loc_dof)
    resize!(nr.global_Δu, total_loc_dof)
    if nr.use_preconditioner
        resize!(nr.global_precond, total_loc_dof)
        # Compute material-based preconditioner at initialization (cheap, no force evals)
        compute_material_preconditioner!(dh, nr.global_precond)
    end
    # Pre-allocate buffers for norm computation to avoid allocations in hot path
    resize!(nr.u_norm_sq_buf, dh.n_chunks)
    resize!(nr.v_norm_sq_buf, dh.n_chunks)
    # Initialize damping diagonal if damping is enabled
    if nr.damping > 0
        resize!(nr.damping_diag, total_loc_dof)
        compute_damping_diagonal!(dh, nr.damping_diag)
    end
    return nothing
end

function init_newton_buffers!(nr::NewtonKrylov, dh::MPIBodyDataHandler)
    n_loc_dof = get_n_loc_dof(dh.chunk.system)
    n_ranks = mpi_nranks()

    # Initialize MPI-specific cached data
    nr.mpi_dof_counts = MPI.Allgather(Int32(n_loc_dof), mpi_comm())
    nr.mpi_displs = zeros(Int32, n_ranks)
    for i in 2:n_ranks
        nr.mpi_displs[i] = nr.mpi_displs[i-1] + nr.mpi_dof_counts[i-1]
    end

    total_dof = sum(nr.mpi_dof_counts)
    resize!(nr.global_residual, total_dof)
    resize!(nr.global_Δu, total_dof)
    if nr.use_preconditioner
        resize!(nr.global_precond, total_dof)
        # Compute material-based preconditioner at initialization (cheap, no force evals)
        compute_material_preconditioner_mpi!(dh, nr)
    end
    # Initialize damping diagonal if damping is enabled
    if nr.damping > 0
        resize!(nr.damping_diag, total_dof)
        compute_damping_diagonal_mpi!(dh, nr)
    end
    # MPI doesn't need norm buffers (uses reduction)
    return nothing
end

# Overload for MultibodySetups to throw an error
function init_time_solver!(::NewtonKrylov, ::AbstractDataHandler)
    msg = "NewtonKrylov solver only implemented for single body setups!\n"
    return throw(ArgumentError(msg))
end

function validate_chunk_for_newton_krylov!(chunk::AbstractBodyChunk)
    # Boundary condition checks
    (; condhandler) = chunk
    n_sdbcs = length(condhandler.single_dim_bcs)
    n_pdsdbcs = length(condhandler.posdep_single_dim_bcs)
    n_dbcs = length(condhandler.data_bcs)
    if n_sdbcs + n_pdsdbcs + n_dbcs > 0
        msg = "NewtonKrylov solver does not support all boundary condition types!\n"
        msg *= "Specify only position-dependent boundary conditions!"
        throw(ArgumentError(msg))
    end

    # Material limitations
    incompatible_mats = (:CRMaterial, :RKCRMaterial)
    if nameof(typeof(chunk.mat)) ∈ incompatible_mats
        msg = "$(nameof(typeof(chunk.mat))) not compatible with Newton-Krylov!\n"
        msg *= "The material uses stress rotation designed for dynamic simulations.\n"
        msg *= "For quasi-static Newton-Krylov simulations, consider using the material"
        msg *= " without stress rotation instead.\n"
        throw(ArgumentError(msg))
    end

    # No damage allowed
    (; mat, system, paramsetup) = chunk
    for i in each_point_idx(system)
        params = get_params(paramsetup, i)
        if has_fracture(mat, params)
            msg = "NewtonKrylov solver currently does not support damage or fracture!\n"
            msg *= "Fracture parameters found at point index $(i)!\n"
            throw(ArgumentError(msg))
        end
    end

    return nothing
end

function newton_krylov_check(nr::NewtonKrylov)
    if nr.end_time < 0
        error("`end_time` of NewtonKrylov smaller than zero!\n")
    end
    if nr.n_steps < 0
        error("`n_steps` of NewtonKrylov smaller than zero!\n")
    end
    if nr.Δt < 0
        error("`Δt` of NewtonKrylov smaller than zero!\n")
    end
    if nr.maxiter < 0
        error("`maxiter` of NewtonKrylov smaller than zero!\n")
    end
    if nr.tol < 0
        error("`tol` of NewtonKrylov smaller than zero!\n")
    end
    if nr.perturbation_scale ≤ 0
        error("`perturbation_scale` of NewtonKrylov must be larger than zero!\n")
    end
    return nothing
end

function solve!(dh::AbstractDataHandler, nr::NewtonKrylov, options::AbstractJobOptions)
    export_reference_results(dh, options)

    # Log solver info
    log_it(options, "NEWTON-KRYLOV ITERATIONS")
    print_log(stdout, "\n") # necessary for correct terminal output

    # Time stepping loop
    for n in 1:nr.n_steps
        newton_krylov_step!(dh, nr, options, n)
    end

    add_to_logfile(options, "\n") # necessary for correct logfile output

    return dh
end

# Newton-Krylov step function - distributed across all chunks
function newton_krylov_step!(dh::AbstractDataHandler, nr::NewtonKrylov,
                              options::AbstractJobOptions, n::Int)
    (; n_steps, Δt, maxiter, tol, damping) = nr

    t = n * Δt
    β = n / n_steps # load step factor

    # Compute damping decay factor: α = 1 at step 1, α = 0 at final step
    # This provides regularization early on that vanishes at the final solution
    if n_steps > 1
        α = (n_steps - n) / (n_steps - 1)
    else
        α = 0.0  # No damping for single-step problems
    end
    damping_factor = damping * α

    apply_newton_bcs_and_pos_update!(dh, β)

    for iter in 1:maxiter
        calc_force_density!(dh, t, Δt)
        calc_residual!(dh)

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
            msg = "Newton-Krylov solver did not converge after max iterations!\n"
            msg *= "Consider increasing `maxiter` or `tol`.\n"
            throw(ErrorException(msg))
        end

        solve_linear_system!(dh, nr, t, Δt, damping_factor)
    end

    newton_export_results(dh, options, n, t)

    return nothing
end

function apply_newton_bcs_and_pos_update!(dh::ThreadsBodyDataHandler, β)
    @threads :static for chunk_id in eachindex(dh.chunks)
        chunk = dh.chunks[chunk_id]
        apply_incr_boundary_conditions!(chunk, β)
        update_position!(chunk.storage, chunk.system, chunk.condhandler.constrained_dofs)
    end
    return nothing
end

function apply_newton_bcs_and_pos_update!(dh::MPIBodyDataHandler, β)
    (; chunk) = dh
    apply_incr_boundary_conditions!(chunk, β)
    update_position!(chunk.storage, chunk.system, chunk.condhandler.constrained_dofs)
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

function calc_residual!(dh::ThreadsBodyDataHandler)
    @threads :static for chunk_id in eachindex(dh.chunks)
        calc_residual!(dh.chunks[chunk_id])
    end
    return nothing
end

function calc_residual!(dh::MPIBodyDataHandler)
    calc_residual!(dh.chunk)
    return nothing
end

function calc_residual!(chunk::AbstractBodyChunk)
    (; system, storage, condhandler) = chunk
    (; residual, b_int, b_ext) = storage
    (; volume) = system

    # Convert 2D force matrices to 1D residual vector (local DOFs only)
    for (dof, _, i) in each_loc_dof_idx(system)
        residual[dof] = (b_int[dof] + b_ext[dof]) * volume[i]
    end

    # Set residual to zero for constrained DOFs
    for dof in condhandler.constrained_dofs
        @inbounds residual[dof] = 0.0
    end

    return nothing
end

function get_residual_norm(dh::ThreadsBodyDataHandler)
    # Compute global norm as sqrt(sum of squared local norms)
    local_norms_sq = zeros(dh.n_chunks)
    @threads :static for chunk_id in eachindex(dh.chunks)
        local_norms_sq[chunk_id] = get_local_residual_norm_sq(dh.chunks[chunk_id])
    end
    return sqrt(sum(local_norms_sq))
end

function get_residual_norm(dh::MPIBodyDataHandler)
    # Compute global norm using MPI reduction
    local_norm_sq = get_local_residual_norm_sq(dh.chunk)
    global_norm_sq = MPI.Allreduce(local_norm_sq, MPI.SUM, mpi_comm())
    return sqrt(global_norm_sq)
end

function get_local_residual_norm_sq(chunk::AbstractBodyChunk)
    residual = chunk.storage.residual
    s = 0.0
    @inbounds for i in each_loc_dof(chunk)
        r = residual[i]
        s += r * r
    end
    return s
end

# Diagonal preconditioner struct that supports ldiv! for IterativeSolvers.jl
struct DiagonalPreconditioner{T<:AbstractVector{Float64}}
    diag_inv::T
end

# ldiv!(y, P, x) computes y = P \ x = P⁻¹ * x
# For diagonal preconditioner P = diag(d), P⁻¹ = diag(1/d), so y[i] = x[i] / d[i]
# Since we store the inverse already: y[i] = x[i] * diag_inv[i]
function LinearAlgebra.ldiv!(y::AbstractVector, P::DiagonalPreconditioner,
                              x::AbstractVector)
    @inbounds for i in eachindex(x)
        y[i] = x[i] * P.diag_inv[i]
    end
    return y
end

# ldiv!(P, x) computes P \ x in-place of x
function LinearAlgebra.ldiv!(P::DiagonalPreconditioner, x::AbstractVector)
    @inbounds for i in eachindex(x)
        x[i] = x[i] * P.diag_inv[i]
    end
    return x
end

# P \ x - creates a new vector
function Base.:\(P::DiagonalPreconditioner, x::AbstractVector)
    y = similar(x)
    ldiv!(y, P, x)
    return y
end

function solve_linear_system!(dh::ThreadsBodyDataHandler, solver::NewtonKrylov, t, Δt,
                               damping_factor::Float64=0.0)
    (; global_residual, global_Δu, global_precond, damping_diag) = solver
    (; gmres_maxiter, gmres_reltol, gmres_abstol, gmres_restart, use_preconditioner) = solver

    # Reset Δu on all chunks
    @threads :static for chunk_id in eachindex(dh.chunks)
        dh.chunks[chunk_id].storage.Δu .= 0.0
    end

    # Calculate total local DOF (same across all Newton iterations)
    total_loc_dof = length(global_residual)

    # Preconditioner was already computed at initialization (material-based, very cheap)
    # No additional work needed per iteration

    # Create distributed Jacobian-free linear operator
    # This operator works on concatenated local DOF vectors from all chunks
    # Uses storage.v_temp and storage.Jv_temp from each chunk
    # If damping is enabled, adds damping_factor * D to the Jacobian
    J_linmap = LinearMap(total_loc_dof; issymmetric=false, ismutating=true) do Jv, v
        # Scatter global v to chunk storage.v_temp
        scatter_to_chunk!(dh, v)

        # Compute distributed JVP (uses storage.v_temp, writes to storage.Jv_temp)
        jacobian_vector_product!(dh, solver, t, Δt)

        # Gather storage.Jv_temp to global Jv
        gather_from_chunk!(Jv, dh)

        # Add damping term to Jacobian: Jv += damping_factor * D * v
        # where D is the diagonal stiffness matrix (same scaling as preconditioner)
        if damping_factor > 0
            @inbounds for i in eachindex(Jv)
                Jv[i] += damping_factor * damping_diag[i] * v[i]
            end
        end
    end

    # Gather residuals from all chunks into global residual vector
    gather_residuals!(global_residual, dh)

    # Add damping term to residual: r += damping_factor * D * u
    # This regularizes the system by penalizing large displacements early on
    if damping_factor > 0
        gather_displacements!(global_Δu, dh)  # temporarily store u in global_Δu
        @inbounds for i in eachindex(global_residual)
            global_residual[i] += damping_factor * damping_diag[i] * global_Δu[i]
        end
    end

    # Reset solution vector
    global_Δu .= 0.0

    # Solve using GMRES with the distributed matrix-free operator
    # Use right preconditioning if enabled: solve J*P⁻¹*y = -r, then x = P⁻¹*y
    if use_preconditioner
        # Create right preconditioner using our custom struct
        P_right = DiagonalPreconditioner(global_precond)

        # Solve with right preconditioning
        gmres!(global_Δu, J_linmap, -global_residual; Pr=P_right,
               maxiter=gmres_maxiter, reltol=gmres_reltol, abstol=gmres_abstol,
               restart=gmres_restart)
    else
        gmres!(global_Δu, J_linmap, -global_residual;
               maxiter=gmres_maxiter, reltol=gmres_reltol, abstol=gmres_abstol,
               restart=gmres_restart)
    end

    # Scatter global_Δu to chunks and update displacements
    scatter_to_chunk!(dh, global_Δu)
    @threads :static for chunk_id in eachindex(dh.chunks)
        chunk = dh.chunks[chunk_id]
        update_displacement_from_Δu!(chunk)
    end

    return nothing
end

# Compute material-based diagonal preconditioner at initialization.
# This is very cheap (no force evaluations) and provides a good scaling for GMRES.
# The preconditioner estimates the Jacobian diagonal using material stiffness:
# For peridynamics: d_i ~ bc * V_i * (number of neighbors) / δ²
# where bc is the bond constant and V_i is the volume associated with point i.
function compute_material_preconditioner!(dh::ThreadsBodyDataHandler,
                                           global_precond::Vector{Float64})
    offset = 0
    for chunk_id in eachindex(dh.chunks)
        chunk = dh.chunks[chunk_id]
        n_loc_dof = compute_chunk_precond!(chunk)
        # Gather to global precond vector
        (; storage, system) = chunk
        for i in each_loc_dof(system)
            @inbounds global_precond[offset + i] = storage.precond_diag[i]
        end
        offset += n_loc_dof
    end
    return nothing
end

function compute_chunk_precond!(chunk::AbstractBodyChunk)
    (; system, storage, paramsetup, condhandler) = chunk
    (; precond_diag) = storage
    (; volume) = system
    (; constrained_dofs) = condhandler

    # Compute diagonal scaling based on material stiffness
    for (dof, _, i) in each_loc_dof_idx(system)
        params = get_params(paramsetup, i)
        # Get stiffness scaling from material parameters
        # For most materials, K (bulk modulus) gives a good scaling
        # d ~ K * volume, but we need to account for the peridynamic kernel
        # A simple estimate: use K * volume / δ as diagonal scaling
        K = params.K
        δ = params.δ
        V_i = volume[i]

        # Stiffness scale: K * V / δ gives units of [force/length]
        # This captures the local stiffness contribution
        d = K * V_i / δ

        # Store inverse for preconditioning (avoid division by zero)
        if abs(d) < 1e-14
            d = 1e-14
        end
        @inbounds precond_diag[dof] = 1.0 / d
    end

    # For constrained DOFs, use identity (d=1, so inv(d)=1)
    for dof in constrained_dofs
        @inbounds precond_diag[dof] = 1.0
    end

    return get_n_loc_dof(system)
end

function compute_material_preconditioner_mpi!(dh::MPIBodyDataHandler, solver::NewtonKrylov)
    chunk = dh.chunk
    (; global_precond, mpi_dof_counts) = solver
    (; storage) = chunk

    # Compute local preconditioner
    compute_chunk_precond!(chunk)

    # Gather to all ranks
    MPI.Allgatherv!(storage.precond_diag, VBuffer(global_precond, mpi_dof_counts), mpi_comm())
    return nothing
end

# Compute damping diagonal based on material stiffness
# The damping term adds: damping * α * D * u to the residual
# where D is a diagonal matrix with stiffness scaling (same as preconditioner but not inverted)
# and α = (n_steps - n) / (n_steps - 1) is the decay factor (1 at step 1, 0 at final step)
function compute_damping_diagonal!(dh::ThreadsBodyDataHandler, damping_diag::Vector{Float64})
    offset = 0
    for chunk_id in eachindex(dh.chunks)
        chunk = dh.chunks[chunk_id]
        (; system, paramsetup, condhandler) = chunk
        (; volume) = system
        (; constrained_dofs) = condhandler

        # Compute stiffness-based damping diagonal
        for (dof, _, i) in each_loc_dof_idx(system)
            params = get_params(paramsetup, i)
            K = params.K
            δ = params.δ
            V_i = volume[i]

            # Same stiffness scaling as preconditioner: K * V / δ
            d = K * V_i / δ
            @inbounds damping_diag[offset + dof] = d
        end

        # Zero damping for constrained DOFs (they're fixed anyway)
        for dof in constrained_dofs
            @inbounds damping_diag[offset + dof] = 0.0
        end

        offset += get_n_loc_dof(system)
    end
    return nothing
end

function compute_damping_diagonal_mpi!(dh::MPIBodyDataHandler, solver::NewtonKrylov)
    chunk = dh.chunk
    (; damping_diag, mpi_dof_counts) = solver
    (; system, paramsetup, condhandler) = chunk
    (; volume) = system
    (; constrained_dofs) = condhandler
    n_loc_dof = get_n_loc_dof(system)

    # Use precond_diag as temporary buffer for local damping values
    local_damping = chunk.storage.precond_diag

    # Compute stiffness-based damping diagonal
    for (dof, _, i) in each_loc_dof_idx(system)
        params = get_params(paramsetup, i)
        K = params.K
        δ = params.δ
        V_i = volume[i]
        d = K * V_i / δ
        @inbounds local_damping[dof] = d
    end

    # Zero damping for constrained DOFs
    for dof in constrained_dofs
        @inbounds local_damping[dof] = 0.0
    end

    # Gather to all ranks
    MPI.Allgatherv!(local_damping, VBuffer(damping_diag, mpi_dof_counts), mpi_comm())
    return nothing
end

# New scatter function using storage.v_temp
function scatter_to_chunk!(dh::ThreadsBodyDataHandler, v::AbstractVector{Float64})
    offset = 0
    for chunk_id in eachindex(dh.chunks)
        (; storage, system) = dh.chunks[chunk_id]
        (; v_temp) = storage
        for i in each_loc_dof(system)
            @inbounds v_temp[i] = v[offset + i]
        end
        offset += get_n_loc_dof(system)
    end
    return nothing
end

# New gather function using storage.Jv_temp
function gather_from_chunk!(v::AbstractVector{Float64}, dh::ThreadsBodyDataHandler)
    offset = 0
    for chunk_id in eachindex(dh.chunks)
        (; storage, system) = dh.chunks[chunk_id]
        (; Jv_temp) = storage
        for i in each_loc_dof(system)
            @inbounds v[offset + i] = Jv_temp[i]
        end
        offset += get_n_loc_dof(system)
    end
    return nothing
end

function gather_residuals!(global_residual::Vector{Float64}, dh::ThreadsBodyDataHandler)
    offset = 0
    for chunk_id in eachindex(dh.chunks)
        (; storage, system) = dh.chunks[chunk_id]
        (; residual) = storage
        for i in each_loc_dof(system)
            @inbounds global_residual[offset + i] = residual[i]
        end
        offset += get_n_loc_dof(system)
    end
    return nothing
end

# Gather displacements from chunks into a global vector (for damping computation)
function gather_displacements!(global_u::Vector{Float64}, dh::ThreadsBodyDataHandler)
    offset = 0
    for chunk_id in eachindex(dh.chunks)
        (; storage, system) = dh.chunks[chunk_id]
        (; displacement) = storage
        for (dof, _, _) in each_loc_dof_idx(system)
            @inbounds global_u[offset + dof] = displacement[dof]
        end
        offset += get_n_loc_dof(system)
    end
    return nothing
end

# New JVP function using storage fields
function jacobian_vector_product!(dh::ThreadsBodyDataHandler, solver::NewtonKrylov, t, Δt)
    # Compute adaptive perturbation: ε = scale * cbrt(eps) * (||u|| + 1) / ||v||
    # For central differences, the optimal ε balances truncation error O(ε²) with
    # roundoff error O(eps_mach/ε), giving ε ~ cbrt(eps_mach)
    ε = compute_perturbation(dh, solver)

    # Step 1: Store original state and apply positive perturbation
    @threads :static for chunk_id in eachindex(dh.chunks)
        store_state_and_apply_perturbation!(dh.chunks[chunk_id], ε)
    end

    # Step 2: Calculate forces at positively perturbed state (with halo exchange)
    calc_force_density!(dh, t, Δt)

    # Step 3: Store positive forces and apply negative perturbation
    @threads :static for chunk_id in eachindex(dh.chunks)
        store_pos_force_and_apply_neg_perturbation!(dh.chunks[chunk_id], ε)
    end

    # Step 4: Calculate forces at negatively perturbed state (with halo exchange)
    calc_force_density!(dh, t, Δt)

    # Step 5: Compute Jv and restore original state
    @threads :static for chunk_id in eachindex(dh.chunks)
        compute_jvp_and_restore!(dh.chunks[chunk_id], ε)
    end

    return nothing
end

# Compute adaptive perturbation using the standard JFNK heuristic for central differences:
# ε = scale * cbrt(eps_mach) * (||u|| + 1) / ||v||
# where eps_mach is machine epsilon, u is current displacement, v is the direction vector
function compute_perturbation(dh::ThreadsBodyDataHandler, solver::NewtonKrylov)
    (; perturbation_scale, u_norm_sq_buf, v_norm_sq_buf) = solver

    # Compute ||u||² and ||v||² in parallel using pre-allocated buffers
    @threads :static for chunk_id in eachindex(dh.chunks)
        chunk = dh.chunks[chunk_id]
        u_sq, v_sq = compute_u_v_norm_sq(chunk)
        @inbounds u_norm_sq_buf[chunk_id] = u_sq
        @inbounds v_norm_sq_buf[chunk_id] = v_sq
    end

    # Sum across chunks (no allocation)
    u_norm_sq = 0.0
    v_norm_sq = 0.0
    @inbounds for i in eachindex(u_norm_sq_buf)
        u_norm_sq += u_norm_sq_buf[i]
        v_norm_sq += v_norm_sq_buf[i]
    end

    u_norm = sqrt(u_norm_sq)
    v_norm = sqrt(v_norm_sq)

    # Avoid division by zero
    if v_norm < eps(Float64)
        v_norm = 1.0
    end

    # cbrt(eps(Float64)) ≈ 6.055e-6 for Float64
    ε = perturbation_scale * cbrt(eps(Float64)) * (u_norm + 1.0) / v_norm
    return ε
end

function compute_perturbation(dh::MPIBodyDataHandler, solver::NewtonKrylov)
    (; perturbation_scale) = solver

    # Compute local norms
    u_sq, v_sq = compute_u_v_norm_sq(dh.chunk)

    # Global reduction
    u_norm_sq = MPI.Allreduce(u_sq, MPI.SUM, mpi_comm())
    v_norm_sq = MPI.Allreduce(v_sq, MPI.SUM, mpi_comm())

    u_norm = sqrt(u_norm_sq)
    v_norm = sqrt(v_norm_sq)

    # Avoid division by zero
    if v_norm < eps(Float64)
        v_norm = 1.0
    end

    ε = perturbation_scale * cbrt(eps(Float64)) * (u_norm + 1.0) / v_norm
    return ε
end

function compute_u_v_norm_sq(chunk::AbstractBodyChunk)
    (; storage, system) = chunk
    (; displacement, v_temp) = storage

    u_sq = 0.0
    v_sq = 0.0
    for (dof, _, _) in each_loc_dof_idx(system)
        @inbounds u_sq += displacement[dof]^2
        @inbounds v_sq += v_temp[dof]^2
    end
    return u_sq, v_sq
end

function store_state_and_apply_perturbation!(chunk::AbstractBodyChunk, ε::Float64)
    (; storage, system) = chunk
    (; b_int, position, displacement, b_int_copy, displacement_copy, v_temp) = storage

    # Store original state
    displacement_copy .= displacement
    b_int_copy .= b_int

    # Apply positive perturbation: u + ε*v (using v_temp)
    for (dof, _, _) in each_loc_dof_idx(system)
        @inbounds displacement[dof] = displacement_copy[dof] + ε * v_temp[dof]
        @inbounds position[dof] = system.position[dof] + displacement[dof]
    end

    # Reset forces for new calculation
    b_int .= 0.0

    return nothing
end

function store_pos_force_and_apply_neg_perturbation!(chunk::AbstractBodyChunk,
                                                              ε::Float64)
    (; storage, system) = chunk
    (; b_int, position, displacement, displacement_copy, temp_force, v_temp) = storage
    (; volume) = system

    # Store F(u + ε*v) * volume in temp_force
    for (dof, _, i) in each_loc_dof_idx(system)
        @inbounds temp_force[dof] = b_int[dof] * volume[i]
    end

    # Apply negative perturbation: u - ε*v (using v_temp)
    for (dof, _, _) in each_loc_dof_idx(system)
        @inbounds displacement[dof] = displacement_copy[dof] - ε * v_temp[dof]
        @inbounds position[dof] = system.position[dof] + displacement[dof]
    end

    # Reset forces for new calculation
    b_int .= 0.0

    return nothing
end

function compute_jvp_and_restore!(chunk::AbstractBodyChunk, ε::Float64)
    (; storage, system, condhandler) = chunk
    (; b_int, position, displacement, displacement_copy, b_int_copy, temp_force) = storage
    (; Jv_temp, v_temp) = storage
    (; volume) = system
    (; constrained_dofs) = condhandler

    # Compute Jacobian-vector product: Jv = (F(u + εv) - F(u - εv)) / (2ε)
    for (dof, _, i) in each_loc_dof_idx(system)
        @inbounds Jv_temp[dof] = (temp_force[dof] - b_int[dof] * volume[i]) / (2ε)
    end

    # Handle constrained DOFs: set Jv[dof] = v[dof] (identity row)
    for dof in constrained_dofs
        @inbounds Jv_temp[dof] = v_temp[dof]
    end

    # Restore original state
    for (dof, _, _) in each_loc_dof_idx(system)
        @inbounds displacement[dof] = displacement_copy[dof]
        @inbounds position[dof] = system.position[dof] + displacement[dof]
    end
    b_int .= b_int_copy

    return nothing
end

function update_displacement_from_Δu!(chunk::AbstractBodyChunk)
    (; storage, system, condhandler) = chunk
    (; position, displacement, v_temp) = storage
    (; free_dofs) = condhandler

    # v_temp now contains the Δu values
    for dof in free_dofs
        @inbounds du = v_temp[dof]
        @inbounds displacement[dof] += du
        @inbounds position[dof] += du
    end

    return nothing
end

# MPI version: Gather global vectors to all ranks, solve GMRES on full system, scatter back
# This ensures all ranks stay synchronized during GMRES iterations
function solve_linear_system!(dh::MPIBodyDataHandler, solver::NewtonKrylov, t, Δt,
                               damping_factor::Float64=0.0)
    chunk = dh.chunk
    (; storage, system) = chunk
    (; gmres_maxiter, gmres_reltol, gmres_abstol, gmres_restart, use_preconditioner) = solver
    (; global_residual, global_Δu, global_precond, mpi_dof_counts, mpi_displs) = solver
    (; damping_diag) = solver

    n_loc_dof = get_n_loc_dof(system)
    total_dof = length(global_residual)
    my_offset = mpi_displs[mpi_chunk_id()]

    # Preconditioner was already computed at initialization (material-based, very cheap)
    # No additional work needed per iteration

    # Use v_temp as local residual buffer (copy residual into it)
    v_temp = storage.v_temp
    for i in 1:n_loc_dof
        @inbounds v_temp[i] = storage.residual[i]
    end

    # Gather local residuals to all ranks
    MPI.Allgatherv!(v_temp, VBuffer(global_residual, mpi_dof_counts), mpi_comm())

    # Add damping term to residual: r += damping_factor * D * u
    if damping_factor > 0
        # First gather displacements
        (; displacement) = storage
        for (dof, _, _) in each_loc_dof_idx(system)
            @inbounds v_temp[dof] = displacement[dof]
        end
        MPI.Allgatherv!(v_temp, VBuffer(global_Δu, mpi_dof_counts), mpi_comm())

        # Add damping to residual
        @inbounds for i in eachindex(global_residual)
            global_residual[i] += damping_factor * damping_diag[i] * global_Δu[i]
        end
    end

    # Reset solution vector
    global_Δu .= 0.0
    storage.Δu .= 0.0

    # Get reference to Jv_temp for use in closure
    (; Jv_temp) = storage

    # Create global Jacobian-free linear operator
    # All ranks compute the same global JVP together
    # If damping is enabled, adds damping_factor * D to the Jacobian
    J_linmap = LinearMap(total_dof; issymmetric=false, ismutating=true) do Jv, v
        # Each rank extracts its local portion into v_temp
        for i in 1:n_loc_dof
            @inbounds v_temp[i] = v[my_offset + i]
        end

        # Compute JVP with halo exchange (all ranks synchronized)
        # Uses v_temp as input, Jv_temp as output
        jacobian_vector_product_mpi!(dh, solver, t, Δt)

        # Gather all local JVPs to global
        MPI.Allgatherv!(Jv_temp, VBuffer(Jv, mpi_dof_counts), mpi_comm())

        # Add damping term to Jacobian: Jv += damping_factor * D * v
        if damping_factor > 0
            @inbounds for i in eachindex(Jv)
                Jv[i] += damping_factor * damping_diag[i] * v[i]
            end
        end
    end

    # Solve using GMRES with right preconditioning if enabled
    if use_preconditioner
        P_right = DiagonalPreconditioner(global_precond)

        gmres!(global_Δu, J_linmap, -global_residual; Pr=P_right,
               maxiter=gmres_maxiter, reltol=gmres_reltol, abstol=gmres_abstol,
               restart=gmres_restart)
    else
        gmres!(global_Δu, J_linmap, -global_residual;
               maxiter=gmres_maxiter, reltol=gmres_reltol, abstol=gmres_abstol,
               restart=gmres_restart)
    end

    # Extract local solution into v_temp and update displacements
    for i in 1:n_loc_dof
        @inbounds v_temp[i] = global_Δu[my_offset + i]
    end
    update_displacement_from_Δu!(chunk)

    return nothing
end

# MPI version of JVP using storage fields
function jacobian_vector_product_mpi!(dh::MPIBodyDataHandler, solver::NewtonKrylov, t, Δt)
    chunk = dh.chunk

    # Compute adaptive perturbation
    ε = compute_perturbation(dh, solver)

    # Step 1: Store original state and apply positive perturbation
    store_state_and_apply_perturbation!(chunk, ε)

    # Step 2: Calculate forces at positively perturbed state (with MPI halo exchange)
    calc_force_density!(dh, t, Δt)

    # Step 3: Store positive forces and apply negative perturbation
    store_pos_force_and_apply_neg_perturbation!(chunk, ε)

    # Step 4: Calculate forces at negatively perturbed state (with MPI halo exchange)
    calc_force_density!(dh, t, Δt)

    # Step 5: Compute Jv and restore original state
    compute_jvp_and_restore!(chunk, ε)

    return nothing
end

function newton_export_results(dh::ThreadsBodyDataHandler, options::JobOptions, n, t)
    @threads :static for chunk_id in eachindex(dh.chunks)
        export_results(dh, options, chunk_id, n, t)
    end
    return nothing
end

function newton_export_results(dh::MPIBodyDataHandler, options::JobOptions, n, t)
    export_results(dh, options, n, t)
    return nothing
end

# Required interface functions
function init_field_solver(::NewtonKrylov, system::AbstractSystem, ::Val{:position})
    return copy(system.position)
end

function init_field_solver(::NewtonKrylov, system::AbstractSystem, ::Val{:displacement})
    return zeros(3, get_n_loc_points(system))
end

function init_field_solver(::NewtonKrylov, system::AbstractSystem, ::Val{:b_int})
    return zeros(3, get_n_points(system))
end

function init_field_solver(::NewtonKrylov, system::AbstractSystem, ::Val{:b_ext})
    return zeros(3, get_n_points(system))
end

# Newton-Krylov specific fields
function init_field_solver(::NewtonKrylov, system::AbstractSystem,
                           ::Val{:displacement_copy})
    return zeros(3, get_n_loc_points(system))
end
function init_field_solver(::AbstractTimeSolver, system::AbstractSystem,
                           ::Val{:displacement_copy})
    return Array{Float64,2}(undef, 0, 0)
end

function init_field_solver(::NewtonKrylov, system::AbstractSystem,
                           ::Val{:b_int_copy})
    return zeros(3, get_n_points(system))
end
function init_field_solver(::AbstractTimeSolver, system::AbstractSystem,
                           ::Val{:b_int_copy})
    return Array{Float64,2}(undef, 0, 0)
end

function init_field_solver(::NewtonKrylov, system::AbstractSystem, ::Val{:residual})
    return zeros(get_n_loc_dof(system))
end
function init_field_solver(::AbstractTimeSolver, system::AbstractSystem, ::Val{:residual})
    return Vector{Float64}()
end

function init_field_solver(::NewtonKrylov, system::AbstractSystem, ::Val{:temp_force})
    return zeros(get_n_loc_dof(system))
end
function init_field_solver(::AbstractTimeSolver, system::AbstractSystem,
                           ::Val{:temp_force})
    return Array{Float64,1}(undef, 0)
end

function init_field_solver(::NewtonKrylov, system::AbstractSystem, ::Val{:Δu})
    return zeros(get_n_loc_dof(system))
end
function init_field_solver(::AbstractTimeSolver, system::AbstractSystem, ::Val{:Δu})
    return Array{Float64,1}(undef, 0)
end

# GMRES temporary buffers for Jacobian-free Newton-Krylov
function init_field_solver(::NewtonKrylov, system::AbstractSystem, ::Val{:v_temp})
    return zeros(get_n_loc_dof(system))
end
function init_field_solver(::AbstractTimeSolver, system::AbstractSystem, ::Val{:v_temp})
    return Vector{Float64}()
end

function init_field_solver(::NewtonKrylov, system::AbstractSystem, ::Val{:Jv_temp})
    return zeros(get_n_loc_dof(system))
end
function init_field_solver(::AbstractTimeSolver, system::AbstractSystem, ::Val{:Jv_temp})
    return Vector{Float64}()
end

# Diagonal preconditioner storage
function init_field_solver(::NewtonKrylov, system::AbstractSystem, ::Val{:precond_diag})
    return zeros(get_n_loc_dof(system))
end
function init_field_solver(::AbstractTimeSolver, system::AbstractSystem,
                           ::Val{:precond_diag})
    return Vector{Float64}()
end

function req_point_data_fields_timesolver(::Type{<:NewtonKrylov})
    # Jacobian-free Newton-Krylov method: no jacobian, temp_force_b, or affected_points
    # Added v_temp, Jv_temp for GMRES temporary buffers
    # Added precond_diag for diagonal preconditioner
    fields = (:position, :displacement, :b_int, :b_ext, :residual, :displacement_copy,
              :temp_force, :b_int_copy, :Δu, :v_temp, :Jv_temp, :precond_diag)
    return fields
end

function req_bond_data_fields_timesolver(::Type{<:NewtonKrylov})
    return ()
end

function req_data_fields_timesolver(::Type{<:NewtonKrylov})
    return ()
end

function log_timesolver(options::AbstractJobOptions, nr::NewtonKrylov)
    msg = "JACOBIAN-FREE NEWTON-KRYLOV SOLVER\n"
    msg *= msg_qty("number of time steps", nr.n_steps)
    msg *= msg_qty("time step size", nr.Δt)
    msg *= msg_qty("simulation time", nr.end_time)
    msg *= msg_qty("max iterations per step", nr.maxiter)
    msg *= msg_qty("convergence tolerance", nr.tol)
    msg *= msg_qty("perturbation scale", nr.perturbation_scale)
    msg *= msg_qty("GMRES max iterations", nr.gmres_maxiter)
    msg *= msg_qty("GMRES relative tolerance", nr.gmres_reltol)
    msg *= msg_qty("GMRES absolute tolerance", nr.gmres_abstol)
    msg *= msg_qty("using preconditioner", nr.use_preconditioner)
    if nr.damping > 0
        msg *= msg_qty("load-step damping", nr.damping)
    end
    log_it(options, msg)
    return nothing
end

# If this is uncommented, the solver is registered and all storages must support it!
# The internal storages in this package, already do, but the user-defined ones may not.
# Therefore, for now this is commented out.
# register_solver!(NewtonKrylov)
