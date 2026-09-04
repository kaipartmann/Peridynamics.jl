"""
    TwoNeighborInteraction

$(internal_api_warning())

Type for two-neighbor interactions.

# Fields

- `oni_j::Int`: One-neighbor interaction of considered point with point j.
- `oni_k::Int`: One-neighbor interaction of considered point with point k.
- `surface::Float64`: Surface spread by this two-neighbor interaction.
"""
struct TwoNeighborInteraction
    oni_j::Int
    oni_k::Int
    surface::Float64
end

"""
    ThreeNeighborInteraction

$(internal_api_warning())

Type for three-neighbor interactions.

# Fields

- `oni_j::Int`: One-neighbor interaction of considered point with point j.
- `oni_k::Int`: One-neighbor interaction of considered point with point k.
- `oni_l::Int`: One-neighbor interaction of considered point with point l.
- `volume::Float64`: Volume spread by this three-neighbor interaction.
"""
struct ThreeNeighborInteraction
    oni_j::Int
    oni_k::Int
    oni_l::Int
    volume::Float64
end

"""
    InteractionSystem

$(extension_api_note())

A peridynamic system type that is mainly designed for continuum-kinematics-inspired
peridynamics [Javili2019](@cite).

Its one-neighbor interactions are the bonds of the shared bond API, so `neighbor`,
`bond_length`, `fail_permit`, `n_neighbors` and `bond_ids` describe them and
[`each_bond_idx`](@ref), [`get_neighbor`](@ref) and [`reference_bond_length`](@ref) read
them. `each_one_ni_idx` and `get_n_one_nis` are the names of the same two functions in the
vocabulary of this system.

`two_nis` and `three_nis` hold the two- and three-neighbor interactions, `n_two_nis` and
`n_three_nis` their number per point, `two_ni_idxs` and `three_ni_idxs` the range of each
point in those vectors, and the three `volume_*_nis` their effective volumes.

# Fields

$(block_table(InteractionSystem))

"""
@system struct InteractionSystem <: AbstractSystem
    position::PointVector{Float64}
    volume::PointScalar
    neighbor::BondScalar{Int}
    bond_length::BondScalar
    fail_permit::BondScalar{Bool}
    n_neighbors::PointScalar{Int}
    bond_ids::PointScalar{UnitRange{Int}}
    two_nis::Vector{TwoNeighborInteraction}
    three_nis::Vector{ThreeNeighborInteraction}
    volume_one_nis::PointScalar
    volume_two_nis::PointScalar
    volume_three_nis::PointScalar
    n_two_nis::PointScalar{Int}
    n_three_nis::PointScalar{Int}
    two_ni_idxs::PointScalar{UnitRange{Int}}
    three_ni_idxs::PointScalar{UnitRange{Int}}
end

function InteractionSystem(body::AbstractBody, pd::PointDecomposition, chunk_id::Int)
    check_system_compat(InteractionSystem, body.mat)
    loc_points = pd.decomp[chunk_id]
    neighbor, bond_length, fail_permit, n_neighbors, bond_ids, chunk_handler = get_bond_data(body, pd, chunk_id)
    position, volume = get_pos_and_vol_chunk(body, chunk_handler.point_ids)
    N, FT = size(position, 1), default_float_type()
    sizes = SystemSizes{N,FT}(chunk_handler, length(neighbor))
    volume_one_nis = alloc_field(PointScalar(), sizes, LocalPoints())
    if has_two_nis(body)
        two_nis, n_two_nis, two_ni_idxs = find_two_nis(body, loc_points, neighbor, bond_ids)
        volume_two_nis = alloc_field(PointScalar(), sizes, LocalPoints())
    else
        two_nis = Vector{TwoNeighborInteraction}()
        n_two_nis = Vector{Int}()
        volume_two_nis = Vector{FT}()
        two_ni_idxs = Vector{UnitRange{Int}}()
    end
    if has_three_nis(body)
        three_nis, n_three_nis, three_ni_idxs = find_three_nis(body, loc_points, neighbor,
                                                               bond_ids)
        volume_three_nis = alloc_field(PointScalar(), sizes, LocalPoints())
    else
        three_nis = Vector{ThreeNeighborInteraction}()
        n_three_nis = Vector{Int}()
        volume_three_nis = Vector{FT}()
        three_ni_idxs = Vector{UnitRange{Int}}()
    end
    system = InteractionSystem{N,FT}(position, volume, neighbor, bond_length, fail_permit,
                                     n_neighbors, bond_ids, two_nis, three_nis,
                                     volume_one_nis, volume_two_nis, volume_three_nis,
                                     n_two_nis, n_three_nis, two_ni_idxs, three_ni_idxs,
                                     chunk_handler)
    return system
end

function get_system(body::AbstractBody{Material}, pd::PointDecomposition,
                    chunk_id::Int) where {Material<:AbstractInteractionSystemMaterial}
    return InteractionSystem(body, pd, chunk_id)
end

@inline function system_type(::AbstractInteractionSystemMaterial,
                             ::Type{FT}=default_float_type(),
                             ::Val{N}=Val(3)) where {FT,N}
    return host_system_type(InteractionSystem, Val(N), FT)
end

function check_system_compat(::Type{S},
                             ::M) where {S<:InteractionSystem,M<:AbstractMaterial}
    msg = "body with material `$(M)` incompatible to `InteractionSystem`!\n"
    msg *= "The material has to be a subtype of `AbstractInteractionSystemMaterial`!\n"
    return throw(ArgumentError(msg))
end

function check_system_compat(::Type{<:InteractionSystem},
                             ::AbstractInteractionSystemMaterial)
    return nothing
end

function get_c2(params::AbstractPointParameters)
    hasproperty(params, :C2) || return 0.0
    return params.C2
end

function get_c3(params::AbstractPointParameters)
    hasproperty(params, :C3) || return 0.0
    return params.C3
end

function has_two_nis(body::AbstractBody)
    for params in body.point_params
        get_c2(params) ≈ 0 || return true
    end
    return false
end

@inline function has_two_nis(chunk::AbstractBodyChunk{<:InteractionSystem})
    return has_two_nis(chunk.paramsetup)
end

@inline function has_two_nis(param_setup::AbstractParameterHandler)
    for params in param_setup.parameters
        get_c2(params) ≈ 0 || return true
    end
    return false
end

@inline function has_two_nis(params::AbstractPointParameters)
    get_c2(params) ≈ 0 || return true
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

@inline function has_three_nis(chunk::AbstractBodyChunk{<:InteractionSystem})
    return has_three_nis(chunk.paramsetup)
end

@inline function has_three_nis(param_setup::AbstractParameterHandler)
    for params in param_setup.parameters
        get_c3(params) ≈ 0 || return true
    end
    return false
end

@inline function has_three_nis(params::AbstractPointParameters)
    get_c3(params) ≈ 0 || return true
    return false
end

@inline function has_three_nis(body::AbstractBody, point_id::Int)
    get_c3(body.point_params[body.params_map[point_id]]) ≈ 0 || return true
    return false
end

function find_two_nis(body, loc_points, neighbor, bond_ids)
    two_nis = Vector{TwoNeighborInteraction}()
    sizehint!(two_nis, n_points(body) * 1000)
    n_two_nis = zeros(Int, length(loc_points))
    two_ni_idxs = fill(0:-1, length(loc_points))
    two_ni_idx_start, two_ni_idx_end = 1, 0
    position = body.position
    for (li, i) in enumerate(loc_points)
        num = 0
        δ = get_point_param(body, :δ, i)
        jk_seen = Set{Tuple{Int,Int}}()
        for oni_j in bond_ids[li], oni_k in bond_ids[li]
            j, k = neighbor[oni_j], neighbor[oni_k]
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
        n_two_nis[li] = num
        two_ni_idx_end = two_ni_idx_start + num - 1
        two_ni_idxs[li] = two_ni_idx_start:two_ni_idx_end
        two_ni_idx_start = two_ni_idx_end + 1
    end
    return two_nis, n_two_nis, two_ni_idxs
end

@inline function surf_two_neigh(ξijx, ξijy, ξijz, ξikx, ξiky, ξikz)
    return sqrt((ξijy * ξikz - ξijz * ξiky)^2 +
                (ξijz * ξikx - ξijx * ξikz)^2 +
                (ξijx * ξiky - ξijy * ξikx)^2)
end

function find_three_nis(body, loc_points, neighbor, bond_ids)
    three_nis = Vector{ThreeNeighborInteraction}()
    sizehint!(three_nis, n_points(body) * 1000)
    n_three_nis = zeros(Int, length(loc_points))
    three_ni_idxs = fill(0:-1, length(loc_points))
    three_ni_idx_start, three_ni_idx_end = 1, 0
    position = body.position
    for (li, i) in enumerate(loc_points)
        num = 0
        δ = get_point_param(body, :δ, i)
        jkl_seen = Set{Tuple{Int,Int,Int}}()
        for oni_j in bond_ids[li], oni_k in bond_ids[li], oni_l in bond_ids[li]
            j = neighbor[oni_j]
            k = neighbor[oni_k]
            l = neighbor[oni_l]
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
        n_three_nis[li] = num
        three_ni_idx_end = three_ni_idx_start + num - 1
        three_ni_idxs[li] = three_ni_idx_start:three_ni_idx_end
        three_ni_idx_start = three_ni_idx_end + 1
    end
    return three_nis, n_three_nis, three_ni_idxs
end

# `each_one_ni_idx` and `get_n_one_nis` are the interaction-system names of the shared bond
# accessors `each_bond_idx` and `get_n_bonds`, kept as aliases because the one-neighbor
# interactions of this system are addressed as one-neighbor interactions in its own API and
# as bonds by the shared kinematics of `current_bond_length`, `bond_stretch` and
# `update_bond_lengths!`. No material of this system inherits `BondLengthCache`, so those
# shared methods always compute the distance here, exactly as the dedicated methods used to.
@inline each_one_ni_idx(is::InteractionSystem, point_id::Int) = each_bond_idx(is, point_id)
@inline get_n_one_nis(is::InteractionSystem) = get_n_bonds(is)

@inline each_two_ni_idx(is::InteractionSystem, point_id::Int) = is.two_ni_idxs[point_id]
@inline each_three_ni_idx(is::InteractionSystem, point_id::Int) = is.three_ni_idxs[point_id]

function initialize!(chunk::AbstractBodyChunk{<:InteractionSystem})
    update_volumes!(chunk)
    return nothing
end

function update_volumes!(chunk::AbstractBodyChunk{<:InteractionSystem})
    volume_hood = get_neighborhood_volume(chunk)
    update_volume_one_nis!(chunk.system, volume_hood)
    has_two_nis(chunk) && update_volume_two_nis!(chunk.system, volume_hood)
    has_three_nis(chunk) && update_volume_three_nis!(chunk.system, volume_hood)
    return nothing
end

@inline function get_neighborhood_volume(chunk::AbstractBodyChunk{<:InteractionSystem})
    system = chunk.system
    δ = [get_params(chunk, i).δ for i in each_point_idx(chunk)]
    full_volume_hoods = 4 / 3 * π .* δ .^ 3
    discrete_volume_hoods = zeros(get_n_loc_points(chunk))
    for i in each_point_idx(chunk)
        volume_hood_point = system.volume[i]
        for bond_id in each_one_ni_idx(system, i)
            j = get_neighbor(system, bond_id)
            volume_hood_point += system.volume[j]
        end
        discrete_volume_hoods[i] = volume_hood_point
    end
    β = discrete_volume_hoods ./ full_volume_hoods
    volume_hood = full_volume_hoods .* β
    return volume_hood
end

function update_volume_one_nis!(system, volume_hood)
    (; volume_one_nis, n_neighbors) = system
    for (i, n) in enumerate(n_neighbors)
        if n > 0
            volume_one_nis[i] = volume_hood[i] / n
        end
    end
    return nothing
end

function update_volume_two_nis!(system, volume_hood)
    (; volume_two_nis, n_two_nis) = system
    for (i, n) in enumerate(n_two_nis)
        if n > 0
            volume_two_nis[i] = volume_hood[i] / n
        end
    end
    return nothing
end

function update_volume_three_nis!(system, volume_hood)
    (; volume_three_nis, n_three_nis) = system
    for (i, n) in enumerate(n_three_nis)
        if n > 0
            volume_three_nis[i] = volume_hood[i] / n
        end
    end
    return nothing
end

function calc_timestep_point(system::InteractionSystem, params::AbstractPointParameters,
                             point_id::Int)
    dtsum = 0.0
    for bond_id in each_one_ni_idx(system, point_id)
        j = get_neighbor(system, bond_id)
        L = reference_bond_length(system, bond_id)
        dtsum += system.volume[j] * params.C1 / L
    end
    return sqrt(2 * params.rho / dtsum)
end

function calc_force_density!(chunk::AbstractBodyChunk{<:InteractionSystem}, t, Δt)
    (; system, mat, paramsetup, storage) = chunk
    calc_force_density!(storage, system, mat, paramsetup, t, Δt)
    return nothing
end

function calc_force_density!(storage::AbstractStorage, system::InteractionSystem,
                             mat::AbstractInteractionSystemMaterial,
                             paramsetup::AbstractParameterSetup, t, Δt)
    (; dmgmodel) = mat
    storage.b_int .= 0
    for i in each_point_idx(system)
        calc_failure!(storage, system, mat, dmgmodel, paramsetup, t, Δt, i)
        calc_damage!(storage, system, mat, dmgmodel, paramsetup, i)
        force_density_point!(storage, system, mat, paramsetup, t, Δt, i)
    end
    return nothing
end

function calc_damage!(chunk::AbstractBodyChunk{<:InteractionSystem})
    (; system, mat, paramsetup, storage) = chunk
    (; dmgmodel) = mat
    for point_id in each_point_idx(chunk)
        calc_damage!(storage, system, mat, dmgmodel, paramsetup, point_id)
    end
    return nothing
end

function log_msg_interaction_system(n_one_nis::Int, n_two_nis::Int, n_three_nis::Int)
    msg = msg_qty("number of one-neighbor-interactions", n_one_nis)
    msg *= msg_qty("number of two-neighbor-interactions", n_two_nis)
    msg *= msg_qty("number of three-neighbor-interactions", n_three_nis)
    return msg
end

function log_system(::Type{I}, options::AbstractJobOptions,
                    dh::AbstractDataHandler) where {I<:InteractionSystem}
    n_one_nis, n_two_nis, n_three_nis = calc_n_interactions(dh)
    msg = "INTERACTION SYSTEM"
    body_name = string(get_body_name(dh))
    isempty(body_name) || (msg *= " `" * body_name * "`")
    msg *= "\n"
    msg *= msg_qty("number of one-neighbor-interactions", n_one_nis)
    msg *= msg_qty("number of two-neighbor-interactions", n_two_nis)
    msg *= msg_qty("number of three-neighbor-interactions", n_three_nis)
    log_it(options, msg)
    return nothing
end

function calc_n_interactions(dh::AbstractThreadsBodyDataHandler)
    n_one_nis = 0
    n_two_nis = 0
    n_three_nis = 0
    for chunk in dh.chunks
        (; two_nis, three_nis) = chunk.system
        n_one_nis += get_n_one_nis(chunk.system)
        n_two_nis += length(two_nis)
        n_three_nis += length(three_nis)
    end
    return n_one_nis, n_two_nis, n_three_nis
end

function calc_n_interactions(dh::AbstractMPIBodyDataHandler)
    n_one_nis = MPI.Reduce(get_n_one_nis(dh.chunk.system), MPI.SUM, mpi_comm())
    n_two_nis = MPI.Reduce(length(dh.chunk.system.two_nis), MPI.SUM, mpi_comm())
    n_three_nis = MPI.Reduce(length(dh.chunk.system.three_nis), MPI.SUM, mpi_comm())
    return n_one_nis, n_two_nis, n_three_nis
end

# the interaction count is system knowledge, so this initial value serves the field inside
# a damage state and as a flat storage field alike; `one_ni_active` and `damage` need no
# hook, their initial values follow from the declarations of `InteractionFracFields`
function init_field_system(system::InteractionSystem, ::Val{:n_active_one_nis})
    return copy(system.n_neighbors)
end

function required_point_parameters(::Type{<:AbstractInteractionSystemMaterial})
    return (:δ, :rho, elasticity_parameters()..., :C1, :C2, :C3)
end

function get_interaction_parameters(mat::AbstractInteractionSystemMaterial, params;
                                    C1=nothing, C2=nothing, C3=nothing)
    (; δ, μ, λ) = params
    _C1::Float64 = isnothing(C1) ? 0.0 : float(C1)
    _C2::Float64 = isnothing(C2) ? 0.0 : float(C2)
    _C3::Float64 = isnothing(C3) ? 0.0 : float(C3)

    if _C1 ≈ 0 && _C2 ≈ 0 && _C3 ≈ 0
        _C1 = 30 / π * μ / δ^4
        _C2 = 0.0
        _C3 = 32 / π^4 * (λ - μ) / δ^12
    else
        msg = "interaction parameters for $(typeof(mat)) specified manually!\n"
        msg *= "Be careful when adjusting these parameters to avoid unexpected outcomes!"
        @mpiroot @warn msg
    end

    return (; C1=_C1, C2=_C2, C3=_C3)
end

function log_material_property(::Val{:dmgmodel}, mat::AbstractInteractionSystemMaterial;
                               indentation::Int=2)
    return log_dmgmodel(mat.dmgmodel; indentation)
end

"""
    InteractionParameters

$(extension_api_note())

Parameter block of the three material constants of the continuum-kinematics-inspired
formulation. They are resolved together, because they are either all derived from the
elastic parameters or all specified by hand. See [`@params_fields`](@ref).

$(block_table(InteractionParameters))
"""
@params_fields InteractionParameters begin
    @derived (; C1, C2, C3) = get_interaction_parameters(mat, (; δ, μ, λ); C1, C2, C3)
    @log "parameter one-neighbor interactions" C1
    @log "parameter two-neighbor interactions" C2
    @log "parameter three-neighbor interactions" C3
end
