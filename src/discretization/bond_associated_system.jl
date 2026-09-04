"""
    BondAssociatedSystem

$(internal_api_warning())

A [`BondSystem`](@ref) that also knows, for every bond, which bonds of the same point lie
within the bond horizon `δb` of its neighbor. Those are the bond-associated family of the
bond, and they are stored flat: `intersection_bond_ids` holds the bond indices of all
families one after the other, and `intersection_ids[bond_id]` is the range of that vector
which belongs to `bond_id`, see `each_intersecting_bond_idx`.

# Fields

$(block_table(BondAssociatedSystem))

`ba_hood_volume` is the volume of the bond-associated family of a bond and `ba_volume_sum`
the sum of those volumes over all bonds of a point.
"""
@system struct BondAssociatedSystem <: AbstractBondSystem
    position::PointVector{Float64}
    volume::PointScalar
    neighbor::BondScalar{Int}
    bond_length::BondScalar
    fail_permit::BondScalar{Bool}
    n_neighbors::PointScalar{Int}
    bond_ids::PointScalar{UnitRange{Int}}
    intersection_bond_ids::Vector{Int}
    intersection_ids::BondScalar{UnitRange{Int}}
    ba_volume_sum::PointScalar
    ba_hood_volume::BondScalar
    kernels::BondScalar
end

function BondAssociatedSystem(body::AbstractBody, pd::PointDecomposition, chunk_id::Int)
    check_system_compat(BondAssociatedSystem, body.mat)
    neighbor, bond_length, fail_permit, n_neighbors, bond_ids, chunk_handler = get_bond_data(body, pd, chunk_id)
    position, volume = get_pos_and_vol_chunk(body, chunk_handler.point_ids)
    N, FT = size(position, 1), default_float_type()
    sizes = SystemSizes{N,FT}(chunk_handler, length(neighbor))
    intersection_bond_ids, intersection_ids = find_intersection_bond_ids(body, position,
                                                                        chunk_handler.loc_points,
                                                                        neighbor, bond_ids)
    ba_volume_sum = alloc_field(PointScalar(), sizes, HaloPoints())
    ba_hood_volume = alloc_field(BondScalar(), sizes, LocalPoints())
    kernels = find_kernels(body, chunk_handler, neighbor, bond_length, bond_ids)
    bas = BondAssociatedSystem{N,FT}(position, volume, neighbor, bond_length, fail_permit,
                                     n_neighbors, bond_ids, intersection_bond_ids,
                                     intersection_ids, ba_volume_sum, ba_hood_volume, kernels,
                                     chunk_handler)
    return bas
end

function get_system(body::AbstractBody{Material}, pd::PointDecomposition,
                    chunk_id::Int) where {Material<:AbstractBondAssociatedSystemMaterial}
    return BondAssociatedSystem(body, pd, chunk_id)
end

@inline function system_type(::AbstractBondAssociatedSystemMaterial,
                             ::Type{FT}=default_float_type(),
                             ::Val{N}=Val(3)) where {FT,N}
    return host_system_type(BondAssociatedSystem, Val(N), FT)
end

function check_system_compat(::Type{S},
                             ::M) where {S<:BondAssociatedSystem,M<:AbstractMaterial}
    msg = "body with material `$(M)` incompatible to `BondAssociatedSystem`!\n"
    msg *= "The material has to be a subtype of `AbstractBondAssociatedSystemMaterial`!\n"
    return throw(ArgumentError(msg))
end

function check_system_compat(::Type{<:BondAssociatedSystem},
                             ::AbstractBondAssociatedSystemMaterial)
    return nothing
end

# `position` has to be the chunk-local matrix and not `body.position`, because
# `get_bond_data` already localized the neighbors. Only for the first chunk both are the
# same, all others would get the bond-associated families of the wrong points.
function find_intersection_bond_ids(body, position, loc_points, neighbor, bond_ids)
    intersection_bond_ids = Vector{Int}()
    sizehint!(intersection_bond_ids, 10 * length(neighbor))
    intersection_ids = fill(1:0, length(neighbor))
    for (li, i) in enumerate(loc_points)
        δb = get_point_param(body, :δb, i)
        δb² = δb * δb
        bond_ids_of_i = bond_ids[li]
        for bond_id in bond_ids_of_i
            j = neighbor[bond_id]
            # no system in scope here, and this system is 3D until 2D physics lands
            Xj = get_vector(position, j, Val(3))
            first_id = length(intersection_bond_ids) + 1
            for other_bond_id in bond_ids_of_i
                jj = neighbor[other_bond_id]
                Xjj = get_vector(position, jj, Val(3))
                ΔX = Xj - Xjj
                L² = dot(ΔX, ΔX)
                if L² < δb²
                    push!(intersection_bond_ids, other_bond_id)
                end
            end
            intersection_ids[bond_id] = first_id:length(intersection_bond_ids)
        end
    end
    return intersection_bond_ids, intersection_ids
end

"""
    each_intersecting_bond_idx(system, point_id, bond_id)

$(internal_api_warning())

The bond indices of the bond-associated family of bond `bond_id` of point `point_id`, i.e.
the bonds of the point whose neighbor lies within the bond horizon `δb` of the neighbor of
`bond_id`. These are indices of the chunk, so they address every bond field of the system
and of the storage.
"""
@inline function each_intersecting_bond_idx(system::BondAssociatedSystem, point_id::Int,
                                            bond_id::Int)
    return view(system.intersection_bond_ids, system.intersection_ids[bond_id])
end

# `ba_hood_volume`: volume of the bond-associated family of a bond.
# `ba_volume_sum`: sum of those volumes over all bonds of a point. Note that this is much
# larger than the volume of the family of the point, because the families overlap.
function calc_ba_volumes!(chunk::AbstractBodyChunk{<:BondAssociatedSystem})
    (; system) = chunk
    (; volume, ba_volume_sum, ba_hood_volume) = system

    for i in each_point_idx(chunk)
        _volume_sum = 0.0
        for bond_idx in each_bond_idx(system, i)
            _ba_hood_volume = 0.0
            for i_bond_idx in each_intersecting_bond_idx(system, i, bond_idx)
                jj = get_neighbor(system, i_bond_idx)
                _ba_hood_volume += volume[jj]
            end
            ba_hood_volume[bond_idx] = _ba_hood_volume
            _volume_sum += _ba_hood_volume
        end
        ba_volume_sum[i] = _volume_sum
    end

    return nothing
end

@inline get_ba_volume_sum(chunk::AbstractBodyChunk) = chunk.system.ba_volume_sum

function initialize!(dh::AbstractThreadsBodyDataHandler{<:BondAssociatedSystem},
                     solver::AbstractTimeSolver)
    @threads :static for chunk in dh.chunks
        calc_ba_volumes!(chunk)
    end
    @threads :static for chunk_id in eachindex(dh.chunks)
        exchange_loc_to_halo!(get_ba_volume_sum, dh, chunk_id)
    end
    calc_force_density!(dh, 0.0, solver.Δt)
    return nothing
end

function initialize!(dh::AbstractMPIBodyDataHandler{<:BondAssociatedSystem},
                     solver::AbstractTimeSolver)
    calc_ba_volumes!(dh.chunk)
    exchange_loc_to_halo!(get_ba_volume_sum, dh)
    calc_force_density!(dh, 0.0, solver.Δt)
    return nothing
end

# Share of the strain energy of a point carried by one of its bonds. These shares have to
# sum to one over the bonds of a point, otherwise the stiffness of the material is scaled by
# whatever they add up to.
@inline function volume_fraction_factor(system::BondAssociatedSystem, point_idx::Int,
                                        bond_idx::Int)
    return system.ba_hood_volume[bond_idx] / system.ba_volume_sum[point_idx]
end

function required_point_parameters(::Type{<:AbstractBondAssociatedSystemMaterial})
    return (:δ, :δb, :rho, elasticity_parameters()...)
end

function get_bond_horizon(δ::Float64; bond_horizon=nothing)
    δb::Float64 = isnothing(bond_horizon) ? δ : float(bond_horizon)
    if δb ≤ 0
        throw(ArgumentError("`bond_horizon` should be larger than zero!\n"))
    end
    if δb < δ
        @warn "a small bond horizon < δ will possibly lead to numerical instabilities!"
    end
    return (; δb)
end

"""
    BondHorizonParameters

$(extension_api_note())

Parameter block of the bond horizon `δb` of a bond-associated material. It defaults to the
horizon `δ`, so the block has to follow the one that provides it. See
[`@params_fields`](@ref).

$(block_table(BondHorizonParameters))
"""
@params_fields BondHorizonParameters begin
    @derived (; δb) = get_bond_horizon(δ; bond_horizon)
    @log "bond horizon" δb
end
