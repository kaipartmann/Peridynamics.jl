"""
    NOSBMaterial <: AbstractPDMaterial
"""
struct NOSBMaterial <: AbstractPDMaterial
    δ::Float64
    rho::Float64
    E::Float64
    nu::Float64
    G::Float64
    K::Float64
    bc::Float64
    Gc::Float64
    εc::Float64
    hs::Float64
end


function NOSBMaterial(; horizon::Real, rho::Real, E::Real, nu::Real, Gc::Real = -1,
    epsilon_c::Real = -1, stabilization::Real=0.01)
    G = E / (2 * (1 + nu))
    K = E / (3 * (1 - 2 * nu))
    bc = 18 * K / (π * horizon^4)
    if (Gc !== -1) && (epsilon_c == -1)
        epsilon_c = sqrt(5.0 * Gc / (9.0 * K * horizon))
    elseif (Gc == -1) && (epsilon_c !== -1)
        Gc = 9.0 / 5.0 * K * horizon * epsilon_c^2
    elseif (Gc !== -1) && (epsilon_c !== -1)
        msg = "Duplicate definition: define either Gc or epsilon_c, not both!"
        throw(ArgumentError(msg))
    elseif (Gc == -1) && (epsilon_c == -1)
        throw(ArgumentError("Either Gc or epsilon_c have to be defined!"))
    end
    return NOSBMaterial(horizon, rho, E, nu, G, K, bc, Gc, epsilon_c, stabilization)
end

function Base.show(io::IO, ::MIME"text/plain", mat::NOSBMaterial)
    print(io, typeof(mat), ":")
    for field in fieldnames(typeof(mat))
        msg = "\n  " * rpad(string(field) * ":", 5) * string(getfield(mat, field))
        print(io, msg)
    end
    return nothing
end

"""
    NOSBBody
"""
struct NOSBBody <: AbstractPDBody
    n_points::Int
    n_bonds::Int
    n_threads::Int
    unique_bonds::Bool
    owned_points::Vector{UnitRange{Int}}
    owned_bonds::Vector{UnitRange{Int}}
    single_tids::Vector{Tuple{Int, Int}}
    multi_tids::Vector{Tuple{Int, Vector{Int}}}
    bond_data::Vector{Tuple{Int, Int, Float64, Bool}}
    volume::Vector{Float64}
    cells::Vector{MeshCell{VTKCellType, Tuple{Int64}}}
    n_family_members::Vector{Int}
    n_active_family_members::Matrix{Int}
    hood_range::Vector{UnitRange{Int}}
    refposition::Matrix{Float64}
    position::Matrix{Float64}
    displacement::Matrix{Float64}
    velocity::Matrix{Float64}
    velocity_half::Matrix{Float64}
    acceleration::Matrix{Float64}
    b_int::Array{Float64, 3}
    b_ext::Matrix{Float64}
    damage::Vector{Float64}
    bond_failure::Vector{Int}
    sigma::Matrix{Float64}
end

function NOSBBody(mat::PDMaterial{NOSBMaterial}, pc::PointCloud)
    BLAS.set_num_threads(1)
    n_threads = nthreads()
    n_points = pc.n_points
    @assert n_points>=n_threads "n_points < n_threads"
    owned_points = defaultdist(n_points, n_threads)
    volume = pc.volume
    cells = get_cells(n_points)
    refposition = copy(pc.position)
    position = copy(pc.position)
    displacement = zeros(Float64, (3, n_points))
    velocity = zeros(Float64, (3, n_points))
    velocity_half = zeros(Float64, (3, n_points))
    acceleration = zeros(Float64, (3, n_points))
    forcedensity_ext = zeros(Float64, (3, n_points))
    sigma = zeros(9, n_points)
    damage = zeros(Int, n_points)
    forcedensity_int = zeros(Float64, (3, n_points, n_threads))
    unique_bonds = false
    bond_data, n_family_members = find_bonds(pc, mat, owned_points)
    n_bonds = length(bond_data)
    owned_bonds = defaultdist(n_bonds, n_threads)
    hood_range = find_hood_range(n_family_members, pc.n_points)
    bond_failure = ones(Int, n_bonds)
    n_active_family_members = zeros(Int, (n_points, n_threads))
    _sum_tids = zeros(Bool, (n_points, n_threads))
    _sum_tids .= false
    @threads for tid in 1:n_threads
        # for current_bond in owned_bonds[tid]
        #     (i, j, _, _) = bond_data[current_bond]
        #     n_active_family_members[i, tid] += 1
        #     n_active_family_members[j, tid] += 1
        #     _sum_tids[i, tid] = true
        #     _sum_tids[j, tid] = true
        # end
        for i in owned_points[tid]
            for current_bond in hood_range[i]
                _, j, _, _ = bond_data[current_bond]
                n_active_family_members[i, 1] += 1
                _sum_tids[i, 1] = true
                _sum_tids[j, tid] = true
            end
        end
    end
    sum_tids = [findall(row) for row in eachrow(_sum_tids)]
    single_tids, multi_tids = find_tids(sum_tids)
    return NOSBBody(n_points, n_bonds, n_threads, unique_bonds, owned_points,
                    owned_bonds, single_tids, multi_tids, bond_data, volume, cells,
                    n_family_members, n_active_family_members, hood_range, refposition, position,
                    displacement, velocity, velocity_half, acceleration, forcedensity_int,
                    forcedensity_ext, damage, bond_failure, sigma)
end

function Base.show(io::IO, ::MIME"text/plain", body::NOSBBody)
    println(io, body.n_points, "-points ", typeof(body), ":")
    print(io, "  ", body.n_bonds, " bonds")
    return nothing
end

init_body(mat::PDMaterial{NOSBMaterial}, pc::PointCloud) = NOSBBody(mat, pc)

function compute_forcedensity!(body::NOSBBody, mat::PDMaterial{NOSBMaterial})
    body.b_int .= 0.0
    body.n_active_family_members .= 0

    @inbounds @threads for tid in 1:body.n_threads

        # initializations
        ΔXij = @MVector zeros(3)
        Δxij = @MVector zeros(3)
        K = @MMatrix zeros(3,3)
        _F = @MMatrix zeros(3,3)

        # loop over all material points
        for i in body.owned_points[tid]
            # calculate the deformation gradient F
            # K and _F need to start from zero for every point
            K .= 0.0 # symmetric shape tensor
            _F .= 0.0 # part of deformation gradient (calculated inside hood-loop)
            for current_bond in body.hood_range[i]
                _, j, L, _ = body.bond_data[current_bond]
                ΔXij .= @views body.refposition[:, j] - body.refposition[:, i]
                Δxij .= @views body.position[:, j] - body.position[:, i]
                Vj = body.volume[j]
                ωij = (1 + mat[i].δ / L) * body.bond_failure[current_bond]
                K .+= ωij .* ΔXij * ΔXij' * Vj
                _F .+= ωij .* Δxij * ΔXij' * Vj
            end
            Kinv = inv(K)
            F = _F * Kinv
            J = det(F)

            # calculate second piola kirchhoff stress
            C = F' * F
            Cinv = inv(C)
            μ = mat[i].G
            κ = mat[i].K
            S = μ .* (I - 1/3 .* tr(C) .* Cinv) .* J^(-2/3) .+ κ/4 .* (J^2 - J^(-2)) .* Cinv
            P = F * S
            sigma = 1/J * P * F'
            body.sigma[:,i] .= sigma[:]

            # calculate PK⁻¹
            PK⁻¹ = P * Kinv

            # calculate force densities
            for current_bond in body.hood_range[i]
                _, j, L, failureflag = body.bond_data[current_bond]

                ΔXij .= @views body.refposition[:, j] - body.refposition[:, i]
                Δxij .= @views body.position[:, j] - body.position[:, i]
                l = sqrt(Δxij[1] * Δxij[1] + Δxij[2] * Δxij[2] + Δxij[3] * Δxij[3])
                ε = (l - L) / L

                # stabilization
                T = mat[i].hs * mat[i].bc .* ε .* Δxij ./ l

                # failure mechanism
                if ε > mat[i].εc
                    if failureflag
                        body.bond_failure[current_bond] = 0
                    end
                end
                body.n_active_family_members[i, 1] += body.bond_failure[current_bond]

                # update of force density
                ωij = (1 + mat[i].δ / L) * body.bond_failure[current_bond]
                tij = ωij * PK⁻¹ * ΔXij + T
                body.b_int[:, i, 1] .+= tij * body.volume[j]
                body.b_int[:, j, tid] .-= tij * body.volume[i]
            end
        end
    end

    return nothing
end

function export_vtk(body::NOSBBody, expfile::String, timestep::Int, time::Float64)
    filename = @sprintf("%s_t%04d", expfile, timestep)
    vtkfile = vtk_grid(filename, body.position, body.cells)
    vtkfile["Damage", VTKPointData()] = body.damage
    vtkfile["ForceDensity", VTKPointData()] = @views body.b_int[:, :, 1]
    vtkfile["Displacement", VTKPointData()] = body.displacement
    vtkfile["Velocity", VTKPointData()] = body.velocity
    vtkfile["TraceOfStress", VTKPointData()] = @views body.sigma[1,:] + body.sigma[5,:] + body.sigma[9,:]
    vtkfile["Time", VTKFieldData()] = time
    vtk_save(vtkfile)
    return nothing
end
