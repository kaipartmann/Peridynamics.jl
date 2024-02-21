struct ChunkHandler
    point_ids::Vector{Int}
    loc_points::UnitRange{Int}
    halo_points::Vector{Int}
    halo_by_src::Dict{Int,UnitRange{Int}}
    localizer::Dict{Int,Int}
end

function ChunkHandler(bonds::Vector{Bond}, pd::PointDecomposition, chunk_id::Int)
    loc_points = pd.decomp[chunk_id]
    halo_points = find_halo_points(bonds, loc_points)
    halo_by_src = sort_halo_by_src!(halo_points, pd.point_src, length(loc_points))
    point_ids = vcat(loc_points, halo_points)
    localizer = find_localizer(point_ids)
    return ChunkHandler(point_ids, loc_points, halo_points, halo_by_src, localizer)
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

function sort_halo_by_src!(halo_points::Vector{Int}, point_src::Dict{Int,Int},
                           n_loc_points::Int)
    halo_sources = [point_src[i] for i in halo_points]
    idx_sorted = sortperm(halo_sources)
    halo_sources .= halo_sources[idx_sorted]
    halo_points .= halo_points[idx_sorted]

    halo_by_src = Dict{Int,UnitRange{Int}}()
    unique_sources = unique(halo_sources)

    for source in unique_sources
        idxs = findall(x -> x == source, halo_sources)
        idx_begin = first(idxs) + n_loc_points
        idx_end = last(idxs) + n_loc_points
        halo_by_src[source] = idx_begin:idx_end
    end

    return halo_by_src
end

function find_localizer(point_ids::Vector{Int})
    localizer = Dict{Int,Int}()
    for (li, i) in enumerate(point_ids)
        localizer[i] = li
    end
    return localizer
end

function localize!(point_ids::Vector{Int}, localizer::Dict{Int,Int})
    for i in eachindex(point_ids)
        point_ids[i] = localizer[point_ids[i]]
    end
    return nothing
end

function localize(point_ids::Vector{Int}, ch::ChunkHandler)
    is_loc_point = zeros(Bool, length(point_ids))
    for i in eachindex(point_ids)
        glob_index = point_ids[i]
        if in(glob_index, ch.loc_points)
            is_loc_point[i] = true
        end
    end
    loc_point_ids = point_ids[is_loc_point]
    localize!(loc_point_ids, ch.localizer)
    return loc_point_ids
end

function localize!(bonds::Vector{Bond}, localizer::Dict{Int,Int})
    for i in eachindex(bonds)
        bond = bonds[i]
        bonds[i] = Bond(localizer[bond.neighbor], bond.length, bond.fail_permit)
    end
    return nothing
end

function localized_point_sets(point_sets::Dict{Symbol,Vector{Int}}, ch::ChunkHandler)
    loc_point_sets = Dict{Symbol,Vector{Int}}()
    for (name, ids) in point_sets
        loc_point_sets[name] = Vector{Int}()
        for id in ids
            if id in ch.loc_points
                push!(loc_point_sets[name], ch.localizer[id])
            end
        end
    end
    return loc_point_sets
end

@inline function each_point_idx(ch::ChunkHandler)
    return eachindex(ch.loc_points)
end
