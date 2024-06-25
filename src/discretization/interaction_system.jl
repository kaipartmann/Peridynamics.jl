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
    bonds, n_one_nis = find_bonds(body, pd.decomp[chunk_id])
    bond_ids = find_bond_ids(n_one_nis)
    if has_two_nis(body)
        two_nis, n_two_nis = find_two_nis(body.position, bonds, bond_ids)
    else
        two_nis = Vector{TwoNeighborInteraction}()
        n_two_nis = Vector{Float64}()
    end
    if has_two_nis(body)
        three_nis, n_three_nis = find_three_nis(body.position, bonds, bond_ids)
    else
        three_nis = Vector{ThreeNeighborInteraction}()
        n_three_nis = Vector{Float64}()
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

function has_three_nis(body::AbstractBody)
    for params in body.point_params
        get_c3(params) ≈ 0 || return true
    end
    return false
end
