struct ChunkLocalHandler{D<:AbstractDiscretization,S<:AbstractStorage}
    point_ids::Vector{Int}
    loc_points::UnitRange{Int}
    halo_points::Vector{Int}
    localizer::Dict{Int,Int}
    d::D
    s::S
end

function chunk_local_handlers(pc::PointCloud, mat::AbstractMaterialConfig, n_chunks::Int)
    point_decomp = defaultdist(pc.n_points, n_chunks)
    vclh = Vector{ChunkLocalHandler}(undef, n_chunks)

    @threads for cid in eachindex(point_decomp)
        loc_points = point_decomp[cid]
        vclh[cid] = ChunkLocalHandler(pc, mat, loc_points)
    end

    return vclh
end

#TODO: dispatch on time solver for the storage!
function ChunkLocalHandler(pc::PointCloud, mat::AbstractMaterialConfig,
                           loc_points::UnitRange{Int})
    bonds, n_neighbors = find_bonds(pc, mat, loc_points)
    bond_range = find_bond_range(n_neighbors)
    halo_points = find_halo_points(bonds, loc_points)
    point_ids = vcat(loc_points, halo_points)
    localizer = find_localizer(point_ids)
    localize!(bonds, localizer)
    position, volume = get_pos_and_vol(pc, point_ids)
    pbd = PointBondDiscretization(position, volume, bonds, n_neighbors, bond_range)
    s = init_storage(mat, pbd, loc_points, halo_points) #TODO: dispatch on time solver!
    return ChunkLocalHandler(point_ids, loc_points, halo_points, localizer, pbd, s)
end
