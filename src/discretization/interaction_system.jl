struct TwoNeighborInteraction
    oni_j::Int
    oni_k::Int
    surface::Float64
end

struct ThreeNeighborInteraction
    oni_j::Int
    oni_k::Int
    oni_l::Int
    volume::Float64
end

struct InteractionSystem <: AbstractSystem
    position::Matrix{Float64}
    one_nis::Vector{Bond}
    two_nis::Vector{TwoNeighborInteraction}
    three_nis::Vector{ThreeNeighborInteraction}
    volume::Vector{Float64}
    volume_one_nis::Vector{Float64}
    volume_two_nis::Vector{Float64}
    volume_three_nis::Vector{Float64}
    n_one_nis::Vector{Int}
    n_two_nis::Vector{Int}
    n_three_nis::Vector{Int}
    bond_ids::Vector{UnitRange{Int}}
end

function InteractionSystem(body::AbstractBody, pd::PointDecomposition, chunk_id::Int)
    check_interaction_system_compat(body.mat)
    loc_points = pd.decomp[chunk_id]
    bonds, n_one_nis = find_bonds(body, loc_points)
    bond_ids = find_bond_ids(n_one_nis)
    volume_one_nis = zeros(length(n_one_nis))
    if has_two_nis(body)
        two_nis, n_two_nis = find_two_nis(body, loc_points, bonds, bond_ids)
        volume_two_nis = zeros(length(n_two_nis))
    else
        two_nis = Vector{TwoNeighborInteraction}()
        n_two_nis = Vector{Float64}()
        volume_two_nis = Vector{Float64}()
    end
    if has_three_nis(body)
        three_nis, n_three_nis = find_three_nis(body, loc_points, bonds, bond_ids)
        volume_three_nis = zeros(length(n_three_nis))
    else
        three_nis = Vector{ThreeNeighborInteraction}()
        n_three_nis = Vector{Float64}()
        volume_three_nis = Vector{Float64}()
    end
    ch = get_chunk_handler(bonds, pd, chunk_id)
    localize!(bonds, ch.localizer)
    position, volume = get_pos_and_vol_chunk(body, ch.point_ids)
    is = InteractionSystem(position, bonds, two_nis, three_nis, volume, volume_one_nis,
                           volume_two_nis, volume_three_nis, n_one_nis, n_two_nis,
                           n_three_nis, bond_ids)
    return is, ch
end

function get_system(body::AbstractBody{Material}, pd::PointDecomposition,
                    chunk_id::Int) where {Material<:AbstractInteractionSystemMaterial}
    return InteractionSystem(body, pd, chunk_id)
end

@inline function system_type(::AbstractInteractionSystemMaterial)
    return InteractionSystem
end

function check_interaction_system_compat(::M) where {M<:AbstractMaterial}
    msg = "body with material `$(M)` incompatible to `InteractionSystem`!\n"
    msg *= "The material has to be a subtype of `AbstractInteractionSystemMaterial`!\n"
    return throw(ArgumentError(msg))
end

function check_interaction_system_compat(::AbstractInteractionSystemMaterial)
    return nothing
end

function has_two_nis(body::AbstractBody)
    for params in body.point_params
        get_c2(params) ≈ 0 || return true
    end
    return false
end

@inline function has_two_nis(body::AbstractBody, point_id::Int)
    get_c2(body.point_params[body.params_map[point_id]]) ≈ 0 || return true
    return false
end

function has_three_nis(body::AbstractBody)
    for params in body.point_params
        get_c3(params) ≈ 0 || return true
    end
    return false
end

@inline function has_three_nis(body::AbstractBody, point_id::Int)
    get_c3(body.point_params[body.params_map[point_id]]) ≈ 0 || return true
    return false
end

function find_two_nis(body, loc_points, bonds, bond_ids)
    two_nis = Vector{TwoNeighborInteraction}()
    sizehint!(two_nis, n_points(body) * 1000)
    n_two_nis = zeros(length(loc_points))
    position = body.position
    points_with_two_nis = filter(x -> has_two_nis(body, x), loc_points)
    for i in points_with_two_nis
        num = 0
        δ = get_point_param(body, :δ, i)
        jk_seen = Set{Tuple{Int,Int}}()
        for oni_j in bond_ids[i], oni_k in bond_ids[i]
            j, k = bonds[oni_j].neighbor, bonds[oni_k].neighbor
            if k !== j && !in((j, k), jk_seen)
                Ξijx = position[1, j] - position[1, i]
                Ξijy = position[2, j] - position[2, i]
                Ξijz = position[3, j] - position[3, i]
                Ξikx = position[1, k] - position[1, i]
                Ξiky = position[2, k] - position[2, i]
                Ξikz = position[3, k] - position[3, i]
                Ξjkx = position[1, k] - position[1, j]
                Ξjky = position[2, k] - position[2, j]
                Ξjkz = position[3, k] - position[3, j]
                Ξjk = sqrt(Ξjkx * Ξjkx + Ξjky * Ξjky + Ξjkz * Ξjkz)
                surface = surf_two_neigh(Ξijx, Ξijy, Ξijz, Ξikx, Ξiky, Ξikz)
                if surface > eps() && Ξjk <= δ
                    num += 1
                    push!(two_nis, TwoNeighborInteraction(oni_j, oni_k, surface))
                    push!(jk_seen, (k, j))
                end
            end
        end
        n_two_nis[i] = num
    end
    return two_nis, n_two_nis
end

@inline function surf_two_neigh(ξijx, ξijy, ξijz, ξikx, ξiky, ξikz)
    return sqrt((ξijy * ξikz - ξijz * ξiky)^2 +
                (ξijz * ξikx - ξijx * ξikz)^2 +
                (ξijx * ξiky - ξijy * ξikx)^2)
end

function find_three_nis(body, loc_points, bonds, bond_ids)
    three_nis = Vector{ThreeNeighborInteraction}()
    sizehint!(three_nis, n_points(body) * 1000)
    n_three_nis = zeros(length(loc_points))
    position = body.position
    points_with_three_nis = filter(x -> has_three_nis(body, x), loc_points)
    for i in points_with_three_nis
        num = 0
        δ = get_point_param(body, :δ, i)
        jkl_seen = Set{Tuple{Int,Int,Int}}()
        for oni_j in bond_ids[i], oni_k in bond_ids[i], oni_l in bond_ids[i]
            j = bonds[oni_j].neighbor
            k = bonds[oni_k].neighbor
            l = bonds[oni_l].neighbor
            if k !== j && l !== j && l !== k && !in((j, k, l), jkl_seen)
                Ξijx = position[1, j] - position[1, i]
                Ξijy = position[2, j] - position[2, i]
                Ξijz = position[3, j] - position[3, i]
                Ξikx = position[1, k] - position[1, i]
                Ξiky = position[2, k] - position[2, i]
                Ξikz = position[3, k] - position[3, i]
                Ξilx = position[1, l] - position[1, i]
                Ξily = position[2, l] - position[2, i]
                Ξilz = position[3, l] - position[3, i]
                Ξjkx = position[1, k] - position[1, j]
                Ξjky = position[2, k] - position[2, j]
                Ξjkz = position[3, k] - position[3, j]
                Ξjlx = position[1, l] - position[1, j]
                Ξjly = position[2, l] - position[2, j]
                Ξjlz = position[3, l] - position[3, j]
                Ξlkx = position[1, k] - position[1, l]
                Ξlky = position[2, k] - position[2, l]
                Ξlkz = position[3, k] - position[3, l]
                _Ξjk = sqrt(Ξjkx * Ξjkx + Ξjky * Ξjky + Ξjkz * Ξjkz)
                _Ξjl = sqrt(Ξjlx * Ξjlx + Ξjly * Ξjly + Ξjlz * Ξjlz)
                _Ξlk = sqrt(Ξlkx * Ξlkx + Ξlky * Ξlky + Ξlkz * Ξlkz)
                Aijkx = Ξijy * Ξikz - Ξijz * Ξiky
                Aijky = Ξijz * Ξikx - Ξijx * Ξikz
                Aijkz = Ξijx * Ξiky - Ξijy * Ξikx
                volume = abs(Aijkx * Ξilx + Aijky * Ξily + Aijkz * Ξilz)
                if volume > eps() && _Ξjk <= δ && _Ξjl <= δ && _Ξlk <= δ
                    num += 1
                    tni = ThreeNeighborInteraction(oni_j, oni_k, oni_l, volume)
                    push!(three_nis, tni)
                    push!(jkl_seen, (l, j, k))
                    push!(jkl_seen, (k, l, j))
                    push!(jkl_seen, (j, l, k))
                    push!(jkl_seen, (l, k, j))
                    push!(jkl_seen, (k, j, l))
                end
            end
        end
        n_three_nis[i] = num
    end
    return three_nis, n_three_nis
end

function initialize!(chunk::AbstractBodyChunk{InteractionSystem})
    update_volumes!(chunk.system)
    return nothing
end

function update_volumes!(system::InteractionSystem)
    return nothing
end
