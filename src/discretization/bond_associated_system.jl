
struct BondAssociatedSystem <: AbstractBondSystem
    position::Matrix{Float64}
    volume::Vector{Float64}
    bonds::Vector{Bond}
    n_neighbors::Vector{Int}
    bond_ids::Vector{UnitRange{Int}}
    intersection_bond_ids::Vector{Vector{Int}}
    hood_volume::Vector{Float64}
    ba_hood_volume::Vector{Float64}
    kernels::Vector{Float64}
    chunk_handler::ChunkHandler
end

function BondAssociatedSystem(body::AbstractBody, pd::PointDecomposition, chunk_id::Int)
    check_bond_associated_system_compat(body.mat)
    loc_points = pd.decomp[chunk_id]
    bonds, n_neighbors = find_bonds(body, loc_points)
    bond_ids = find_bond_ids(n_neighbors)
    intersection_bond_ids = find_intersection_bond_ids(body, loc_points, bonds, bond_ids)
    chunk_handler = get_chunk_handler(bonds, pd, chunk_id)
    localize!(bonds, chunk_handler.localizer)
    position, volume = get_pos_and_vol_chunk(body, chunk_handler.point_ids)
    hood_volume = zeros(get_n_points(chunk_handler))
    ba_hood_volume = zeros(length(bonds))
    kernels = find_kernels(body, chunk_handler, bonds, bond_ids)
    bas = BondAssociatedSystem(position, volume, bonds, n_neighbors, bond_ids,
                               intersection_bond_ids, hood_volume, ba_hood_volume, kernels,
                               chunk_handler)
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

function find_intersection_bond_ids(body, loc_points, bonds, bond_ids)
    intersection_bond_ids = Vector{Vector{Int}}(undef, length(bonds))
    for (li, i) in enumerate(loc_points)
        δb = get_point_param(body, :δb, i)
        δb² = δb * δb
        bond_ids_of_i = bond_ids[li]
        for bond_id in bond_ids_of_i
            bond = bonds[bond_id]
            j = bond.neighbor
            Xj = get_vector(body.position, j)
            intersecting_bonds = Vector{Int}()
            for (ibond_id, bond_id) in enumerate(bond_ids_of_i)
                bond = bonds[bond_id]
                jj = bond.neighbor
                Xjj = get_vector(body.position, jj)
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

function calc_hood_volumes!(chunk::AbstractBodyChunk{<:BondAssociatedSystem})
    (; system) = chunk
    (; volume, bonds, hood_volume, ba_hood_volume) = system

    for i in each_point_idx(chunk)
        _hood_volume = volume[i]
        for bond_idx in each_bond_idx(system, i)
            bond = bonds[bond_idx]
            j = bond.neighbor
            _hood_volume += volume[j]
            _ba_hood_volume = 0.0
            for i_bond_idx in each_intersecting_bond_idx(system, i, bond_idx)
                i_bond = bonds[i_bond_idx]
                jj = i_bond.neighbor
                _ba_hood_volume += volume[jj]
            end
            ba_hood_volume[bond_idx] = _ba_hood_volume
        end
        hood_volume[i] = _hood_volume
    end

    return nothing
end

@inline get_hood_volume(chunk::AbstractBodyChunk) = chunk.system.hood_volume

function initialize!(dh::AbstractThreadsBodyDataHandler{<:BondAssociatedSystem},
                     ::AbstractTimeSolver)
    @threads :static for chunk in dh.chunks
        calc_hood_volumes!(chunk)
    end
    @threads :static for chunk_id in eachindex(dh.chunks)
        exchange_loc_to_halo!(get_hood_volume, dh, chunk_id)
    end
    return nothing
end

function initialize!(dh::AbstractMPIBodyDataHandler{<:BondAssociatedSystem},
                     ::AbstractTimeSolver)
    calc_hood_volumes!(dh.chunk)
    exchange_loc_to_halo!(get_hood_volume, dh)
    return nothing
end

@inline function volume_fraction_factor(system::BondAssociatedSystem, point_idx::Int,
                                        bond_idx::Int)
    return system.ba_hood_volume[bond_idx] / system.hood_volume[point_idx]
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
