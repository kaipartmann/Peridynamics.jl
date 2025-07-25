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

$(internal_api_warning())

A peridynamic system type that is mainly designed for continuum-kinematics-inspired
peridynamics [Javili2019](@cite).

# Fields

- `position::Matrix{Float64}`: Positions of all points of the system.
- `one_nis::Vector{Bond}`: Vector containing all one-neighbor interactions (bonds) of the
    system.
- `two_nis::Vector{TwoNeighborInteraction}`: Vector containing all two-neighbor
    interactions of the system.
- `three_nis::Vector{ThreeNeighborInteraction}`: Vector containing all three-neighbor
    interactions of the system.
- `volume::Vector{Float64}`: Volumes of the points of the system.
- `volume_one_nis::Vector{Float64}`: Effective volumes of one-neighbor interactions.
- `volume_two_nis::Vector{Float64}`: Effective volumes of two-neighbor interactions.
- `volume_three_nis::Vector{Float64}`: Effective volumes of three-neighbor interactions.
- `n_one_nis::Vector{Int}`: Number of one-neighbor interactions for each point of the
    system.
- `n_two_nis::Vector{Int}`: Number of two-neighbor interactions for each point of the
    system.
- `n_three_nis::Vector{Int}`: Number of three-neighbor interactions for each point of the
    system.
- `one_ni_idxs::Vector{UnitRange{Int}}`: Range of the one-neighbor interactions vector
    containing interactions of considered point.
- `two_ni_idxs::Vector{UnitRange{Int}}`: Range of the two-neighbor interactions vector
    containing interactions of considered point.
- `three_ni_idxs::Vector{UnitRange{Int}}`: Range of the three-neighbor interactions vector
    containing interactions of considered point.
- `chunk_handler::ChunkHandler`: Type to handle the chunks for the simulation.
    See [`ChunkHandler`](@ref).
"""
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
    one_ni_idxs::Vector{UnitRange{Int}}
    two_ni_idxs::Vector{UnitRange{Int}}
    three_ni_idxs::Vector{UnitRange{Int}}
    chunk_handler::ChunkHandler
end

function InteractionSystem(body::AbstractBody, pd::PointDecomposition, chunk_id::Int)
    check_interaction_system_compat(body.mat)
    loc_points = pd.decomp[chunk_id]
    bonds, n_one_nis = find_bonds(body, loc_points)
    one_ni_idxs = find_bond_ids(n_one_nis)
    volume_one_nis = zeros(length(n_one_nis))
    if has_two_nis(body)
        two_nis, n_two_nis, two_ni_idxs = find_two_nis(body, loc_points, bonds, one_ni_idxs)
        volume_two_nis = zeros(length(n_two_nis))
    else
        two_nis = Vector{TwoNeighborInteraction}()
        n_two_nis = Vector{Float64}()
        volume_two_nis = Vector{Float64}()
        two_ni_idxs = Vector{UnitRange{Int}}()
    end
    if has_three_nis(body)
        three_nis, n_three_nis, three_ni_idxs = find_three_nis(body, loc_points, bonds,
                                                               one_ni_idxs)
        volume_three_nis = zeros(length(n_three_nis))
    else
        three_nis = Vector{ThreeNeighborInteraction}()
        n_three_nis = Vector{Float64}()
        volume_three_nis = Vector{Float64}()
        three_ni_idxs = Vector{UnitRange{Int}}()
    end
    chunk_handler = get_chunk_handler(bonds, pd, chunk_id)
    localize!(bonds, chunk_handler.localizer)
    position, volume = get_pos_and_vol_chunk(body, chunk_handler.point_ids)
    system = InteractionSystem(position, bonds, two_nis, three_nis, volume, volume_one_nis,
                               volume_two_nis, volume_three_nis, n_one_nis, n_two_nis,
                               n_three_nis, one_ni_idxs, two_ni_idxs, three_ni_idxs,
                               chunk_handler)
    return system
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

@inline function has_two_nis(chunk::AbstractBodyChunk{InteractionSystem})
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

@inline function has_three_nis(chunk::AbstractBodyChunk{InteractionSystem})
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

function find_two_nis(body, loc_points, bonds, bond_ids)
    two_nis = Vector{TwoNeighborInteraction}()
    sizehint!(two_nis, n_points(body) * 1000)
    n_two_nis = zeros(length(loc_points))
    two_ni_idxs = fill(0:-1, length(loc_points))
    two_ni_idx_start, two_ni_idx_end = 1, 0
    position = body.position
    for (li, i) in enumerate(loc_points)
        num = 0
        δ = get_point_param(body, :δ, i)
        jk_seen = Set{Tuple{Int,Int}}()
        for oni_j in bond_ids[li], oni_k in bond_ids[li]
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

function find_three_nis(body, loc_points, bonds, bond_ids)
    three_nis = Vector{ThreeNeighborInteraction}()
    sizehint!(three_nis, n_points(body) * 1000)
    n_three_nis = zeros(length(loc_points))
    three_ni_idxs = fill(0:-1, length(loc_points))
    three_ni_idx_start, three_ni_idx_end = 1, 0
    position = body.position
    for (li, i) in enumerate(loc_points)
        num = 0
        δ = get_point_param(body, :δ, i)
        jkl_seen = Set{Tuple{Int,Int,Int}}()
        for oni_j in bond_ids[li], oni_k in bond_ids[li], oni_l in bond_ids[li]
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
        n_three_nis[li] = num
        three_ni_idx_end = three_ni_idx_start + num - 1
        three_ni_idxs[li] = three_ni_idx_start:three_ni_idx_end
        three_ni_idx_start = three_ni_idx_end + 1
    end
    return three_nis, n_three_nis, three_ni_idxs
end

@inline each_one_ni_idx(is::InteractionSystem, point_id::Int) = is.one_ni_idxs[point_id]
@inline each_bond_idx(is::InteractionSystem, point_id::Int) = each_one_ni_idx(is, point_id)

@inline each_two_ni_idx(is::InteractionSystem, point_id::Int) = is.two_ni_idxs[point_id]
@inline each_three_ni_idx(is::InteractionSystem, point_id::Int) = is.three_ni_idxs[point_id]

function initialize!(chunk::AbstractBodyChunk{InteractionSystem})
    update_volumes!(chunk)
    return nothing
end

function update_volumes!(chunk::AbstractBodyChunk{InteractionSystem})
    volume_hood = get_neighborhood_volume(chunk)
    update_volume_one_nis!(chunk.system, volume_hood)
    has_two_nis(chunk) && update_volume_two_nis!(chunk.system, volume_hood)
    has_three_nis(chunk) && update_volume_three_nis!(chunk.system, volume_hood)
    return nothing
end

@inline function get_neighborhood_volume(chunk::AbstractBodyChunk{InteractionSystem})
    system = chunk.system
    δ = [get_params(chunk, i).δ for i in each_point_idx(chunk)]
    full_volume_hoods = 4 / 3 * π .* δ .^ 3
    discrete_volume_hoods = zeros(n_loc_points(chunk))
    for i in each_point_idx(chunk)
        volume_hood_point = system.volume[i]
        for bond_id in each_one_ni_idx(system, i)
            one_ni = system.one_nis[bond_id]
            j = one_ni.neighbor
            volume_hood_point += system.volume[j]
        end
        discrete_volume_hoods[i] = volume_hood_point
    end
    β = discrete_volume_hoods ./ full_volume_hoods
    volume_hood = full_volume_hoods .* β
    return volume_hood
end

function update_volume_one_nis!(system, volume_hood)
    (; volume_one_nis, n_one_nis) = system
    for (i, n) in enumerate(n_one_nis)
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

function break_bonds!(storage::AbstractStorage, system::InteractionSystem,
                      set_a::Vector{Int}, set_b::Vector{Int})
    storage.n_active_one_nis .= 0
    for point_id in each_point_idx(system)
        for bond_id in each_one_ni_idx(system, point_id)
            bond = system.one_nis[bond_id]
            neighbor_id = bond.neighbor
            point_in_a = in(point_id, set_a)
            point_in_b = in(point_id, set_b)
            neigh_in_a = in(neighbor_id, set_a)
            neigh_in_b = in(neighbor_id, set_b)
            if (point_in_a && neigh_in_b) || (point_in_b && neigh_in_a)
                storage.one_ni_active[bond_id] = false
            end
            storage.n_active_one_nis[point_id] += storage.one_ni_active[bond_id]
        end
    end
    return nothing
end

function calc_timestep_point(system::InteractionSystem, params::AbstractPointParameters,
                             point_id::Int)
    dtsum = 0.0
    for bond_id in each_one_ni_idx(system, point_id)
        one_ni = system.one_nis[bond_id]
        dtsum += system.volume[one_ni.neighbor] * params.C1 / one_ni.length
    end
    return sqrt(2 * params.rho / dtsum)
end

function calc_force_density!(chunk::AbstractBodyChunk{<:InteractionSystem}, t, Δt)
    (; system, mat, paramsetup, storage) = chunk
    (; dmgmodel) = mat
    storage.b_int .= 0
    storage.n_active_one_nis .= 0
    for point_id in each_point_idx(chunk)
        calc_failure!(storage, system, mat, dmgmodel, paramsetup, point_id)
        calc_damage!(storage, system, mat, dmgmodel, paramsetup, point_id)
        force_density_point!(storage, system, mat, paramsetup, t, Δt, point_id)
    end
    nancheck(chunk, t, Δt)
    return nothing
end

function calc_failure!(storage::AbstractStorage, system::InteractionSystem,
                       mat::AbstractInteractionSystemMaterial, dmgmodel::CriticalStretch,
                       paramsetup::AbstractParameterSetup, i)
    (; εc) = get_params(paramsetup, i)
    (; position, n_active_one_nis, one_ni_active) = storage
    (; one_nis) = system
    for bond_id in each_one_ni_idx(system, i)
        one_ni = one_nis[bond_id]
        j, L = one_ni.neighbor, one_ni.length
        Δxij = get_vector_diff(position, i, j)
        l = norm(Δxij)
        ε = (l - L) / L
        if ε > εc && one_ni.fail_permit
            one_ni_active[bond_id] = false
        end
        n_active_one_nis[i] += one_ni_active[bond_id]
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

function calc_damage!(storage::AbstractStorage, system::InteractionSystem,
                      mat::AbstractInteractionSystemMaterial, dmgmodel::AbstractDamageModel,
                      paramsetup::AbstractParameterSetup, i)
    @inbounds storage.damage[i] = 1 - storage.n_active_one_nis[i] / system.n_one_nis[i]
end

@inline function one_ni_failure(storage::AbstractStorage, one_ni_id::Int)
    return storage.one_ni_active[one_ni_id]
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
        (; one_nis, two_nis, three_nis) = chunk.system
        n_one_nis += length(one_nis)
        n_two_nis += length(two_nis)
        n_three_nis += length(three_nis)
    end
    return n_one_nis, n_two_nis, n_three_nis
end

function calc_n_interactions(dh::AbstractMPIBodyDataHandler)
    n_one_nis = MPI.Reduce(length(dh.chunk.system.one_nis), MPI.SUM, mpi_comm())
    n_two_nis = MPI.Reduce(length(dh.chunk.system.two_nis), MPI.SUM, mpi_comm())
    n_three_nis = MPI.Reduce(length(dh.chunk.system.three_nis), MPI.SUM, mpi_comm())
    return n_one_nis, n_two_nis, n_three_nis
end

function init_field_system(system::InteractionSystem, ::Val{:one_ni_active})
    return ones(Bool, get_n_one_nis(system))
end

function init_field_system(system::InteractionSystem, ::Val{:n_active_one_nis})
    return copy(system.n_one_nis)
end

function init_field_system(system::InteractionSystem, ::Val{:damage})
    return zeros(get_n_loc_points(system))
end

function req_point_data_fields_fracture(::Type{<:AbstractInteractionSystemMaterial})
    return (:damage, :n_active_one_nis)
end

function req_bond_data_fields_fracture(::Type{<:AbstractInteractionSystemMaterial})
    return (:one_ni_active,)
end

function req_data_fields_fracture(::Type{<:AbstractInteractionSystemMaterial})
    return ()
end

function required_point_parameters(::Type{<:AbstractInteractionSystemMaterial})
    return (:δ, :rho, elasticity_parameters()..., :C1, :C2, :C3)
end

function get_required_point_parameters(mat::AbstractInteractionSystemMaterial,
                                       p::Dict{Symbol,Any})
    params = (; get_horizon(p)..., get_density(p)..., get_elastic_params(p)...)
    params = (; params..., get_interaction_parameters(mat, p, params)...)
    return params
end

function get_interaction_parameters(mat::AbstractInteractionSystemMaterial,
                                    p::Dict{Symbol,Any}, params)
    (; δ, μ, λ) = params
    if haskey(p, :C1)
        C1::Float64 = float(p[:C1])
    else
        C1 = 0.0
    end

    if haskey(p, :C2)
        C2::Float64 = float(p[:C2])
    else
        C2 = 0.0
    end

    if haskey(p, :C3)
        C3::Float64 = float(p[:C3])
    else
        C3 = 0.0
    end

    if C1 ≈ 0 && C2 ≈ 0 && C3 ≈ 0
        C1 = 30 / π * μ / δ^4
        C2 = 0.0
        C3 = 32 / π^4 * (λ - μ) / δ^12
    else
        msg = "interaction parameters for $(typeof(mat)) specified manually!\n"
        msg *= "Be careful when adjusting these parameters to avoid unexpected outcomes!"
        @mpiroot @warn msg
    end

    return (; C1, C2, C3)
end

function allowed_material_kwargs(::AbstractInteractionSystemMaterial)
    kwargs = (discretization_kwargs()..., elasticity_kwargs()..., fracture_kwargs()...,
              :C1, :C2, :C3)
    return kwargs
end

@inline get_n_one_nis(system::InteractionSystem) = length(system.one_nis)

function log_material_property(::Val{:dmgmodel}, mat::AbstractInteractionSystemMaterial;
                               indentation::Int=2)
    msg = msg_qty("damage model type", typeof(mat.dmgmodel); indentation)
    return msg
end

function log_param_property(::Val{:C1}, param; indentation)
    return msg_qty("parameter one-neighbor interactions", param.C1; indentation)
end

function log_param_property(::Val{:C2}, param; indentation)
    return msg_qty("parameter two-neighbor interactions", param.C2; indentation)
end

function log_param_property(::Val{:C3}, param; indentation)
    return msg_qty("parameter three-neighbor interactions", param.C3; indentation)
end
