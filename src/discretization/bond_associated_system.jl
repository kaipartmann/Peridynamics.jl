
struct BondAssociatedSystem <: AbstractBondSystem
    position::Matrix{Float64}
    volume::Vector{Float64}
    bonds::Vector{Bond}
    n_neighbors::Vector{Int}
    bond_ids::Vector{UnitRange{Int}}
    intersection_bond_ids::Vector{Vector{Int}}
    ba_volume_sum::Vector{Float64}
    ba_hood_volume::Vector{Float64}
    kernels::Vector{Float64}
    chunk_handler::ChunkHandler
end

function BondAssociatedSystem(body::AbstractBody, pd::PointDecomposition, chunk_id::Int)
    check_bond_associated_system_compat(body.mat)
    bonds, n_neighbors, bond_ids, chunk_handler = get_bond_data(body, pd, chunk_id)
    position, volume = get_pos_and_vol_chunk(body, chunk_handler.point_ids)
    intersection_bond_ids = find_intersection_bond_ids(body, position,
                                                       chunk_handler.loc_points, bonds,
                                                       bond_ids)
    ba_volume_sum = zeros(get_n_points(chunk_handler))
    ba_hood_volume = zeros(length(bonds))
    kernels = find_kernels(body, chunk_handler, bonds, bond_ids)
    bas = BondAssociatedSystem(position, volume, bonds, n_neighbors, bond_ids,
                               intersection_bond_ids, ba_volume_sum, ba_hood_volume,
                               kernels, chunk_handler)
    return bas
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

# `position` has to be the chunk-local matrix and not `body.position`, because
# `get_bond_data` already localized `bond.neighbor`. Only for the first chunk both are the
# same, all others would get the bond-associated families of the wrong points.
function find_intersection_bond_ids(body, position, loc_points, bonds, bond_ids)
    intersection_bond_ids = Vector{Vector{Int}}(undef, length(bonds))
    for (li, i) in enumerate(loc_points)
        δb = get_point_param(body, :δb, i)
        δb² = δb * δb
        bond_ids_of_i = bond_ids[li]
        for bond_id in bond_ids_of_i
            bond = bonds[bond_id]
            j = bond.neighbor
            Xj = get_vector(position, j)
            intersecting_bonds = Vector{Int}()
            for (ibond_id, bond_id) in enumerate(bond_ids_of_i)
                bond = bonds[bond_id]
                jj = bond.neighbor
                Xjj = get_vector(position, jj)
                ΔX = Xj - Xjj
                L² = dot(ΔX, ΔX)
                if L² < δb²
                    push!(intersecting_bonds, ibond_id)
                end
            end
            intersection_bond_ids[bond_id] = intersecting_bonds
        end
    end
    return intersection_bond_ids
end

@inline function each_intersecting_bond_idx(system::BondAssociatedSystem, point_id::Int,
                                            bond_id::Int)
    return view(each_bond_idx(system, point_id), system.intersection_bond_ids[bond_id])
end

# `ba_hood_volume`: volume of the bond-associated family of a bond.
# `ba_volume_sum`: sum of those volumes over all bonds of a point. Note that this is much
# larger than the volume of the family of the point, because the families overlap.
function calc_ba_volumes!(chunk::AbstractBodyChunk{<:BondAssociatedSystem})
    (; system) = chunk
    (; volume, bonds, ba_volume_sum, ba_hood_volume) = system

    for i in each_point_idx(chunk)
        _volume_sum = 0.0
        for bond_idx in each_bond_idx(system, i)
            _ba_hood_volume = 0.0
            for i_bond_idx in each_intersecting_bond_idx(system, i, bond_idx)
                i_bond = bonds[i_bond_idx]
                jj = i_bond.neighbor
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

function req_point_data_fields_fracture(::Type{<:AbstractBondAssociatedSystemMaterial})
    return (:damage, :n_active_bonds)
end

function req_bond_data_fields_fracture(::Type{<:AbstractBondAssociatedSystemMaterial})
    return (:bond_active,)
end

function req_data_fields_fracture(::Type{<:AbstractBondAssociatedSystemMaterial})
    return ()
end

function required_point_parameters(::Type{<:AbstractBondAssociatedSystemMaterial})
    return (:δ, :δb, :rho, elasticity_parameters()...)
end

function get_required_point_parameters(::AbstractBondAssociatedSystemMaterial,
                                       p::Dict{Symbol,Any})
    δ_params = get_horizon(p)
    δb_params = get_bond_horizon(p, δ_params.δ)
    return (; δ_params..., δb_params..., get_density(p)..., get_elastic_params(p)...)
end

function get_bond_horizon(p::Dict{Symbol,Any}, δ::Float64)
    δb::Float64 = float(get(p, :bond_horizon, δ))
    if δb ≤ 0
        throw(ArgumentError("`bond_horizon` should be larger than zero!\n"))
    end
    if δb < δ
        @warn "a small bond horizon < δ will possibly lead to numerical instabilities!"
    end
    return (; δb)
end

function allowed_material_kwargs(::AbstractBondAssociatedSystemMaterial)
    return (discretization_kwargs()..., elasticity_kwargs()..., fracture_kwargs()...,
            :bond_horizon)
end

function log_param_property(::Val{:δb}, param; indentation)
    return msg_qty("bond horizon", param.δb; indentation)
end
