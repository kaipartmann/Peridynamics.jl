struct ChunkHandler <: AbstractChunkHandler
    n_loc_points::Int
    point_ids::Vector{Int}
    loc_points::UnitRange{Int}
    halo_points::Vector{Int}
    hidxs_by_src::Dict{Int,UnitRange{Int}}
    localizer::Dict{Int,Int}
end

for __field in fieldnames(ChunkHandler)
    local __funcname = Symbol("get_$(__field)")
    local __accessor_func = quote
        @inline function $(__funcname)(chunk_handler::ChunkHandler)
            return chunk_handler.$(__field)
        end
    end
    eval(__accessor_func)
end

function sort_halo_by_src!(halo_points::Vector{Int}, point_src::Dict{Int,Int},
                           n_loc_points::Int)
    halo_sources = [point_src[i] for i in halo_points]
    idx_sorted = sortperm(halo_sources)
    halo_sources .= halo_sources[idx_sorted]
    halo_points .= halo_points[idx_sorted]
    hidxs_by_src = get_hidxs_by_src(halo_sources, n_loc_points)
    return hidxs_by_src
end

function get_hidxs_by_src(halo_sources::Vector{Int}, n_loc_points::Int)
    @assert sort(halo_sources) == halo_sources
    hidxs_by_src = Dict{Int,UnitRange{Int}}()
    unique_sources = unique(halo_sources)

    for source in unique_sources
        idxs = findall(x -> x == source, halo_sources)
        idx_begin = first(idxs) + n_loc_points
        idx_end = last(idxs) + n_loc_points
        hidxs_by_src[source] = idx_begin:idx_end
    end

    return hidxs_by_src
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

@inline function each_point_idx(chunk_handler::ChunkHandler)
    return eachindex(chunk_handler.loc_points)
end

@inline function each_point_idx_pair(chunk_handler::ChunkHandler)
    return enumerate(chunk_handler.loc_points)
end

@inline function get_loc_view(a::Matrix{T}, chunk_handler::ChunkHandler) where {T}
    return view(a, :, 1:chunk_handler.n_loc_points)
end

@inline function get_loc_view(a::Vector{T}, chunk_handler::ChunkHandler) where {T}
    return view(a, 1:chunk_handler.n_loc_points)
end

@inline function get_n_points(chunk_handler::ChunkHandler)
    return length(chunk_handler.point_ids)
end
