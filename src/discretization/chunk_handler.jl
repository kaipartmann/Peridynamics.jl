"""
    ChunkHandler

$(internal_api_warning())

A type to handle a body chunk and its communication to other chunks.

# Fields

- `n_loc_points::Int`: Number of local points that belong to the body chunk.
- `point_ids::Vector{Int}`: Indices of all local and halo points of the chunk.
- `loc_points::UnitRange{Int}`: Indices of local points of the chunk.
- `halo_points::Vector{Int}`: Indices of halo points of the chunk.
- `hidxs_by_src::Dict{Int,UnitRange{Int}}`: Dict specifying the indices of halo Points
    depending on the body chunk they belong to. So `body_chunk => indices`, with indices
    being the indices of the halo points in `point_ids`.
- `localizer::Dict{Int,Int}`: Localizes global indices to local indices in this chunk.
"""
struct ChunkHandler <: AbstractChunkHandler
    n_loc_points::Int
    point_ids::Vector{Int}
    loc_points::UnitRange{Int}
    halo_points::Vector{Int}
    hidxs_by_src::Dict{Int,UnitRange{Int}}
    localizer::Dict{Int,Int}
end

function ChunkHandler(pd::PointDecomposition, halo_points::Vector{Int}, chunk_id::Int)
    loc_points = pd.decomp[chunk_id]
    n_loc_points = length(loc_points)
    hidxs_by_src = sort_halo_by_src!(halo_points, pd.point_src, length(loc_points))
    point_ids = vcat(loc_points, halo_points)
    localizer = find_localizer(point_ids)
    chunk_handler = ChunkHandler(n_loc_points, point_ids, loc_points, halo_points,
                                 hidxs_by_src, localizer)
    return chunk_handler
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

"""
    DeviceChunkHandler

$(internal_api_warning())

What is left of a [`ChunkHandler`](@ref) once a chunk is moved to another array backend with
`Adapt.adapt`: the two point counts, which is everything a kernel reads. The point ids, the
halo bookkeeping and the localizer stay on the host, because the halo exchange, the point
sets and the export are host code and read them there.

# Fields

- `n_loc_points::Int`: Number of local points that belong to the body chunk.
- `n_points::Int`: Number of local and halo points of the body chunk.
"""
struct DeviceChunkHandler <: AbstractChunkHandler
    n_loc_points::Int
    n_points::Int
end

#=
Only a target that really moves arrays gets a `DeviceChunkHandler`. Adapting to something
that leaves an array alone, e.g. `Adapt.adapt(Array, chunk)` on the host, has to return the
chunk handler it was given, so that the whole chunk comes back `===` to itself and nothing
of the halo bookkeeping is lost by an adapt that moves nothing.
=#
function Adapt.adapt_structure(to, ch::ChunkHandler)
    Adapt.adapt(to, ch.point_ids) === ch.point_ids && return ch
    return DeviceChunkHandler(ch.n_loc_points, length(ch.point_ids))
end

@inline function each_point_idx(chunk_handler::ChunkHandler)
    return eachindex(chunk_handler.loc_points)
end

@inline function each_point_idx_pair(chunk_handler::ChunkHandler)
    return enumerate(chunk_handler.loc_points)
end

@inline function get_loc_view(a::AbstractMatrix, chunk_handler::AbstractChunkHandler)
    return view(a, :, 1:get_n_loc_points(chunk_handler))
end

@inline function get_loc_view(a::AbstractVector, chunk_handler::AbstractChunkHandler)
    return view(a, 1:get_n_loc_points(chunk_handler))
end

@inline function get_n_points(chunk_handler::ChunkHandler)
    return length(chunk_handler.point_ids)
end

# the whole shape layer a device chunk answers, see `AbstractSystem` in `core/systems.jl`
@inline get_n_loc_points(ch::DeviceChunkHandler) = ch.n_loc_points
@inline get_n_points(ch::DeviceChunkHandler) = ch.n_points
@inline each_point_idx(ch::DeviceChunkHandler) = Base.OneTo(ch.n_loc_points)
