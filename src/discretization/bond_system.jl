struct Bond
    neighbor::Int
    length::Float64
    fail_permit::Bool
end

struct BondSystem{Correction<:AbstractCorrection} <: AbstractSystem
    position::Matrix{Float64}
    volume::Vector{Float64}
    bonds::Vector{Bond}
    n_neighbors::Vector{Int}
    bond_ids::Vector{UnitRange{Int}}
    correction::Correction
end

function BondSystem(body::AbstractBody, pd::PointDecomposition, chunk_id::Int)
    check_bond_system_compat(body.mat)
    bonds, n_neighbors = find_bonds(body, pd.decomp[chunk_id])
    bond_ids = find_bond_ids(n_neighbors)
    ch = get_chunk_handler(bonds, pd, chunk_id)
    localize!(bonds, ch.localizer)
    position, volume = get_pos_and_vol_chunk(body, ch.point_ids)
    correction = get_correction(body.mat, ch.n_loc_points, length(ch.point_ids),
                                length(bonds))
    bs = BondSystem(position, volume, bonds, n_neighbors, bond_ids, correction)
    return bs, ch
end

function get_system(body::AbstractBody{Material}, pd::PointDecomposition,
                    chunk_id::Int) where {Material<:AbstractBondSystemMaterial}
    return BondSystem(body, pd, chunk_id)
end

@inline function system_type(mat::AbstractBondSystemMaterial)
    return BondSystem{correction_type(mat)}
end

function check_bond_system_compat(::M) where {M<:AbstractMaterial}
    msg = "body with material `$(M)` incompatible to `BondSystem`!\n"
    msg *= "The material has to be a subtype of `AbstractBondSystemMaterial`!\n"
    return throw(ArgumentError(msg))
end

function check_bond_system_compat(::AbstractBondSystemMaterial)
    return nothing
end

function find_bonds(body::AbstractBody, loc_points::UnitRange{Int})
    δmax = maximum_horizon(body)
    nhs = GridNeighborhoodSearch{3}(search_radius=δmax, n_points=body.n_points,
                                    threaded_update=false)
    initialize_grid!(nhs, body.position)
    bonds = Vector{Bond}()
    sizehint!(bonds, body.n_points * 300)
    n_neighbors = zeros(Int, length(loc_points))
    for (li, i) in enumerate(loc_points)
        n_neighbors[li] = find_bonds!(bonds, nhs, body.position, body.fail_permit,
                                      get_point_param(body, :δ, i), i)
    end
    filter_bonds!(bonds, n_neighbors, loc_points, body)
    return bonds, n_neighbors
end

function find_bonds!(bonds::Vector{Bond}, nhs::PointNeighbors.GridNeighborhoodSearch,
                     position::Matrix{Float64}, fail_permit::Vector{Bool}, δ::Float64,
                     point_id::Int)
    n_bonds_pre = length(bonds)
    foreach_neighbor(position, position, nhs, point_id; search_radius=δ) do i, j, _, L
        if i != j
            push!(bonds, Bond(j, L, fail_permit[i] & fail_permit[j]))
        end
    end
    n_neighbors = length(bonds) - n_bonds_pre
    return n_neighbors
end

function filter_bonds!(bonds::Vector{Bond}, n_neighbors::Vector{Int},
                       loc_points::UnitRange{Int}, body::AbstractBody)
    for crack in body.point_sets_precracks
        filter_bonds_by_crack!(bonds, n_neighbors, loc_points, crack, body)
    end
    return nothing
end

function filter_bonds_by_crack!(bonds::Vector{Bond}, n_neighbors::Vector{Int},
                                loc_points::UnitRange{Int}, crack::PointSetsPreCrack,
                                body::AbstractBody)
    filter_bonds(crack) || return nothing
    set_a, set_b = body.point_sets[crack.set_a], body.point_sets[crack.set_b]
    bond_ids = find_bond_ids(n_neighbors)
    bonds_to_delete = fill(false, length(bonds))
    for (loc_point_id, point_id) in enumerate(loc_points)
        for bond_id in bond_ids[loc_point_id]
            bond = bonds[bond_id]
            neighbor_id = bond.neighbor
            point_in_a = in(point_id, set_a)
            point_in_b = in(point_id, set_b)
            neigh_in_a = in(neighbor_id, set_a)
            neigh_in_b = in(neighbor_id, set_b)
            if (point_in_a && neigh_in_b) || (point_in_b && neigh_in_a)
                bonds_to_delete[bond_id] = true
                n_neighbors[loc_point_id] -= 1
            end
        end
    end
    deleteat!(bonds, bonds_to_delete)
    return nothing
end

function find_bond_ids(n_neighbors::Vector{Int})
    bond_ids = fill(0:0, length(n_neighbors))
    bonds_start, bonds_end = 1, 0
    for i in eachindex(n_neighbors)
        bonds_end = bonds_start + n_neighbors[i] - 1
        bond_ids[i] = bonds_start:bonds_end
        bonds_start += n_neighbors[i]
    end
    return bond_ids
end

function get_pos_and_vol_chunk(body::AbstractBody, point_ids::AbstractVector{<:Integer})
    position = body.position[:, point_ids]
    volume = body.volume[point_ids]
    return position, volume
end

function get_chunk_handler(bonds::Vector{Bond}, pd::PointDecomposition, chunk_id::Int)
    loc_points = pd.decomp[chunk_id]
    n_loc_points = length(loc_points)
    halo_points = find_halo_points(bonds, loc_points)
    hidxs_by_src = sort_halo_by_src!(halo_points, pd.point_src, length(loc_points))
    point_ids = vcat(loc_points, halo_points)
    localizer = find_localizer(point_ids)
    return ChunkHandler(n_loc_points, point_ids, loc_points, halo_points, hidxs_by_src,
                        localizer)
end

function find_halo_points(bonds::Vector{Bond}, loc_points::UnitRange{Int})
    halo_points = Vector{Int}()
    for bond in bonds
        j = bond.neighbor
        if !in(j, loc_points) && !in(j, halo_points)
            push!(halo_points, j)
        end
    end
    return halo_points
end

@inline each_bond_idx(bd::BondSystem, point_id::Int) = bd.bond_ids[point_id]

function localize!(bonds::Vector{Bond}, localizer::Dict{Int,Int})
    for i in eachindex(bonds)
        bond = bonds[i]
        bonds[i] = Bond(localizer[bond.neighbor], bond.length, bond.fail_permit)
    end
    return nothing
end

function break_bonds!(s::AbstractStorage, system::BondSystem, ch::ChunkHandler,
                      set_a::Vector{Int}, set_b::Vector{Int})
    s.n_active_bonds .= 0
    for point_id in each_point_idx(ch)
        for bond_id in each_bond_idx(system, point_id)
            bond = system.bonds[bond_id]
            neighbor_id = bond.neighbor
            point_in_a = in(point_id, set_a)
            point_in_b = in(point_id, set_b)
            neigh_in_a = in(neighbor_id, set_a)
            neigh_in_b = in(neighbor_id, set_b)
            if (point_in_a && neigh_in_b) || (point_in_b && neigh_in_a)
                s.bond_active[bond_id] = false
            end
            s.n_active_bonds[point_id] += s.bond_active[bond_id]
        end
    end
    return nothing
end

function calc_timestep_point(bd::BondSystem, params::AbstractPointParameters, point_id::Int)
    dtsum = 0.0
    for bond_id in each_bond_idx(bd, point_id)
        bond = bd.bonds[bond_id]
        dtsum += bd.volume[bond.neighbor] * params.bc / bond.length
    end
    return sqrt(2 * params.rho / dtsum)
end

function calc_force_density!(chunk::AbstractBodyChunk{S,M}) where {S<:BondSystem,M}
    (; system, mat, paramsetup, storage) = chunk
    storage.b_int .= 0
    storage.n_active_bonds .= 0
    for point_id in each_point_idx(chunk)
        params = get_params(paramsetup, point_id)
        force_density_point!(storage, system, mat, params, point_id)
    end
    return nothing
end

@inline function calc_damage!(chunk::AbstractBodyChunk{S,M}) where {S<:BondSystem,M}
    (; n_neighbors) = chunk.system
    (; n_active_bonds, damage) = chunk.storage
    for point_id in each_point_idx(chunk)
        @inbounds damage[point_id] = 1 - n_active_bonds[point_id] / n_neighbors[point_id]
    end
    return nothing
end

function log_system(::Type{B}, options::AbstractJobOptions,
                    dh::AbstractDataHandler) where {B<:BondSystem}
    n_bonds = calc_n_bonds(dh)
    msg = "BOND SYSTEM\n"
    msg *= msg_qty("number of bonds", n_bonds)
    log_it(options, msg)
    return nothing
end

function calc_n_bonds(dh::AbstractThreadsBodyDataHandler)
    n_bonds = 0
    for chunk in dh.chunks
        n_bonds += length(chunk.system.bonds)
    end
    return n_bonds
end

function calc_n_bonds(dh::AbstractMPIBodyDataHandler)
    n_bonds = MPI.Reduce(length(dh.chunk.system.bonds), MPI.SUM, mpi_comm())
    return n_bonds
end
