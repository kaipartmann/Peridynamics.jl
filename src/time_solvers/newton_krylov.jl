"""
    NewtonKrylov(; kwargs...)

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
solver = NewtonKrylov(steps=100, stepsize=1e-4, maxiter=50, tol=1e-6)
```
"""
mutable struct NewtonKrylov <: AbstractTimeSolver
    end_time::Float64
    n_steps::Int
    Δt::Float64
    maxiter::Int
    tol::Float64
    perturbation::Float64
    gmres_maxiter::Int
    gmres_reltol::Float64
    gmres_abstol::Float64
    gmres_restart::Int
    global_residual::Vector{Float64}
    global_Δu::Vector{Float64}
    mpi_dof_counts::Vector{Int32}
    mpi_displs::Vector{Int32}

    function NewtonKrylov(; time::Real=-1, steps::Int=-1, stepsize::Real=1.0,
                            maxiter::Int=100, tol::Real=1e-4, perturbation::Real=-1,
                            gmres_maxiter::Int=200, gmres_reltol::Real=1e-4,
                            gmres_abstol::Real=1e-8, gmres_restart::Int=50)
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

        # Initialize global GMRES buffers as empty (will be resized on first use)
        # MPI buffers also initialized empty
        global_residual = Vector{Float64}()
        global_Δu = Vector{Float64}()
        # MPI-specific cached data (lazily initialized)
        mpi_dof_counts = Vector{Int32}()
        mpi_displs = Vector{Int32}()

        new(end_time, n_steps, stepsize, maxiter, tol, perturbation, gmres_maxiter,
            gmres_reltol, gmres_abstol, gmres_restart, global_residual, global_Δu,
            mpi_dof_counts, mpi_displs)
    end
end

function Base.show(io::IO, @nospecialize(nr::NewtonKrylov))
    print(io, typeof(nr))
    fields = Vector{Symbol}()
    excluded = (:global_residual, :global_Δu, :mpi_dof_counts, :mpi_displs)
    for field in fieldnames(typeof(nr))
        field in excluded && continue
        value = getfield(nr, field)
        if value > 0
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
        excluded = (:global_residual, :global_Δu, :mpi_dof_counts, :mpi_displs)
        for field in fieldnames(typeof(nr))
            field in excluded && continue
            value = getfield(nr, field)
            if value > 0
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
    calculate_perturbation!(nr, dh)
    init_newton_buffers!(nr, dh)
    return nothing
end

# MPI support
function init_time_solver!(nr::NewtonKrylov, dh::MPIBodyDataHandler)
    newton_krylov_check(nr)
    validate_chunk_for_newton_krylov!(dh.chunk)
    calculate_perturbation!(nr, dh)
    init_newton_buffers!(nr, dh)
    return nothing
end

function init_newton_buffers!(nr::NewtonKrylov, dh::ThreadsBodyDataHandler)
    total_loc_dof = sum(get_n_loc_dof(c.system) for c in dh.chunks)
    resize!(nr.global_residual, total_loc_dof)
    resize!(nr.global_Δu, total_loc_dof)
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

function calculate_perturbation!(nr::NewtonKrylov, dh::AbstractDataHandler)
    if nr.perturbation ≤ 0
        Δx = (minimum_volume(dh))^(1/3)
        nr.perturbation = Δx * 1e-5
    end
    return nothing
end

function minimum_volume(dh::ThreadsBodyDataHandler)
    min_vols = fill(Inf, dh.n_chunks)
    @threads :static for chunk_id in eachindex(dh.chunks)
        chunk = dh.chunks[chunk_id]
        vol = chunk.system.volume
        if !isempty(vol)
            min_vols[chunk_id] = minimum(vol)
        end
    end
    return minimum(min_vols)
end

function minimum_volume(dh::MPIBodyDataHandler)
    local_min_vol = minimum(dh.chunk.system.volume)
    global_min_vol = MPI.Allreduce(local_min_vol, MPI.MIN, mpi_comm())
    return global_min_vol
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
    (; n_steps, Δt, maxiter, tol) = nr

    t = n * Δt
    β = n / n_steps # load step factor

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

        solve_linear_system!(dh, nr, t, Δt)
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

function solve_linear_system!(dh::ThreadsBodyDataHandler, solver::NewtonKrylov, t, Δt)
    (; global_residual, global_Δu) = solver
    (; gmres_maxiter, gmres_reltol, gmres_abstol, gmres_restart) = solver

    # Reset Δu on all chunks
    @threads :static for chunk_id in eachindex(dh.chunks)
        dh.chunks[chunk_id].storage.Δu .= 0.0
    end

    # Calculate total local DOF (same across all Newton iterations)
    total_loc_dof = length(global_residual)

    # Create distributed Jacobian-free linear operator
    # This operator works on concatenated local DOF vectors from all chunks
    # Uses storage.v_temp and storage.Jv_temp from each chunk
    J_linmap = LinearMap(total_loc_dof; issymmetric=false, ismutating=true) do Jv, v
        # Scatter global v to chunk storage.v_temp
        scatter_to_chunk!(dh, v)

        # Compute distributed JVP (uses storage.v_temp, writes to storage.Jv_temp)
        jacobian_vector_product!(dh, solver, t, Δt)

        # Gather storage.Jv_temp to global Jv
        gather_from_chunk!(Jv, dh)
    end

    # Gather residuals from all chunks into global residual vector
    gather_residuals!(global_residual, dh)

    # Reset solution vector
    global_Δu .= 0.0

    # Solve using GMRES with the distributed matrix-free operator
    gmres!(global_Δu, J_linmap, -global_residual;
           maxiter=gmres_maxiter, reltol=gmres_reltol, abstol=gmres_abstol,
           restart=gmres_restart)

    # Scatter global_Δu to chunks and update displacements
    scatter_to_chunk!(dh, global_Δu)
    @threads :static for chunk_id in eachindex(dh.chunks)
        chunk = dh.chunks[chunk_id]
        update_displacement_from_Δu!(chunk)
    end

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

# New JVP function using storage fields
function jacobian_vector_product!(dh::ThreadsBodyDataHandler, solver::NewtonKrylov, t, Δt)
    ε = solver.perturbation

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
function solve_linear_system!(dh::MPIBodyDataHandler, solver::NewtonKrylov, t, Δt)
    chunk = dh.chunk
    (; storage, system) = chunk
    (; gmres_maxiter, gmres_reltol, gmres_abstol, gmres_restart) = solver
    (; global_residual, global_Δu, mpi_dof_counts, mpi_displs) = solver

    n_loc_dof = get_n_loc_dof(system)
    total_dof = length(global_residual)
    my_offset = mpi_displs[mpi_chunk_id()]

    # Use v_temp as local residual buffer (copy residual into it)
    v_temp = storage.v_temp
    for i in 1:n_loc_dof
        @inbounds v_temp[i] = storage.residual[i]
    end

    # Gather local residuals to all ranks
    MPI.Allgatherv!(v_temp, VBuffer(global_residual, mpi_dof_counts), mpi_comm())

    # Reset solution vector
    global_Δu .= 0.0
    storage.Δu .= 0.0

    # Get reference to Jv_temp for use in closure
    (; Jv_temp) = storage

    # Create global Jacobian-free linear operator
    # All ranks compute the same global JVP together
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
    end

    # All ranks solve the same global GMRES problem
    gmres!(global_Δu, J_linmap, -global_residual;
           maxiter=gmres_maxiter, reltol=gmres_reltol, abstol=gmres_abstol,
           restart=gmres_restart)

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
    ε = solver.perturbation

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

function req_point_data_fields_timesolver(::Type{<:NewtonKrylov})
    # Jacobian-free Newton-Krylov method: no jacobian, temp_force_b, or affected_points
    # Added v_temp, Jv_temp for GMRES temporary buffers
    fields = (:position, :displacement, :b_int, :b_ext, :residual, :displacement_copy,
              :temp_force, :b_int_copy, :Δu, :v_temp, :Jv_temp)
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
# register_solver!(NewtonKrylov)
