
struct BondAssociatedSystem <: AbstractSystem
    position::Matrix{Float64}
    volume::Vector{Float64}
    bonds::Vector{Bond}
    n_neighbors::Vector{Int}
    bond_ids::Vector{UnitRange{Int}}
    intersection_bond_ids::Vector{Vector{Int}}
    volume_factor::Vector{Float64}
end

function BondAssociatedSystem(body::AbstractBody, pd::PointDecomposition, chunk_id::Int)
    check_bond_associated_system_compat(body.mat)
    loc_points = pd.decomp[chunk_id]
    bonds, n_neighbors = find_bonds(body, loc_points)
    bond_ids = find_bond_ids(n_neighbors)
    intersection_bond_ids = find_intersection_bond_ids(body, loc_points, bonds, bond_ids)
    ch = get_chunk_handler(bonds, pd, chunk_id)
    localize!(bonds, ch.localizer)
    position, volume = get_pos_and_vol_chunk(body, ch.point_ids)
    volume_factor = find_volume_factor(bonds, bond_ids, volume, intersection_bond_ids)
    bs = BondAssociatedSystem(position, volume, bonds, n_neighbors, bond_ids,
                              intersection_bond_ids, volume_factor)
    return bs, ch
end

function get_system(body::AbstractBody{Material}, pd::PointDecomposition,
                    chunk_id::Int) where {Material<:AbstractBondAssociatedSystemMaterial}
    return BondAssociatedSystem(body, pd, chunk_id)
end

@inline function system_type(::AbstractBondAssociatedSystemMaterial)
    return BondAssociatedSystem
end

function check_bond_associated_system_compat(::M) where {M<:AbstractMaterial}
    msg = "body with material `$(M)` incompatible to `BondAssociatedSystem`!\n"
    msg *= "The material has to be a subtype of `AbstractBondAssociatedSystemMaterial`!\n"
    return throw(ArgumentError(msg))
end

function check_bond_associated_system_compat(::AbstractBondAssociatedSystemMaterial)
    return nothing
end

function find_intersection_bond_ids(body, loc_points, bonds, bond_ids)
    intersection_bond_ids = Vector{Vector{Int}}(undef, length(bonds))
    for (li, i) in enumerate(loc_points)
        δ = get_point_param(body, :δ, i)
        δ² = δ * δ
        bond_ids_of_i = bond_ids[li]
        for bond_id in bond_ids_of_i
            bond = bonds[bond_id]
            j = bond.neighbor
            Xj = get_coordinates(body, j)
            intersecting_bonds = Vector{Int}()
            for (ibond_id, bond_id) in enumerate(bond_ids_of_i)
                bond = bonds[bond_id]
                jj = bond.neighbor
                if jj != j
                    Xjj = get_coordinates(body, jj)
                    ΔX = Xj - Xjj
                    L² = dot(ΔX, ΔX)
                    if L² < δ²
                        push!(intersecting_bonds, ibond_id)
                    end
                end
            end
            intersection_bond_ids[bond_id] = intersecting_bonds
        end
    end
    return intersection_bond_ids
end

function find_volume_factor(bonds, bond_ids, volume, intersection_bond_ids)
    volume_Hx = zeros(length(bond_ids))
    volume_hx = zeros(length(bonds))
    for i in eachindex(bond_ids)
        bond_ids_of_i = bond_ids[i]
        vol_Hx = 0.0
        for bond_id in bond_ids_of_i
            intersection_bond_ids_of_bond = intersection_bond_ids[bond_id]
            bond_j = bonds[bond_id]
            j = bond_j.neighbor
            vol_Hx += volume[j]
            vol_hx = 0.0
            for i_bond_id in bond_ids_of_i[intersection_bond_ids_of_bond]
                bond_jj = bonds[i_bond_id]
                jj = bond_jj.neighbor
                vol_hx += volume[jj]
            end
            volume_hx[bond_id] = vol_hx
        end
        volume_Hx[i] = vol_Hx
    end
    volume_factor = zeros(length(bonds))
    for i in eachindex(bond_ids)
        bond_ids_of_i = bond_ids[i]
        for bond_id in bond_ids_of_i
            bond_j = bonds[bond_id]
            j = bond_j.neighbor
            volume_factor[bond_id] = volume_hx[bond_id] / volume_Hx[i]
        end
    end
    return volume_factor
end

@inline function each_bond_idx(system::BondAssociatedSystem, point_id::Int)
    return system.bond_ids[point_id]
end

@inline function each_intersecting_bond_idx(system::BondAssociatedSystem, point_id::Int,
                                            bond_id::Int)
    return view(each_bond_idx(system, point_id), system.intersection_bond_ids[bond_id])
end
