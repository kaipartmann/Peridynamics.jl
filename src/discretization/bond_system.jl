struct Bond
    neighbor::Int
    length::Float64
    fail_permit::Bool
end

struct BondSystem{Correction<:AbstractCorrection} <: AbstractBondSystem
    position::Matrix{Float64}
    volume::Vector{Float64}
    bonds::Vector{Bond}
    n_neighbors::Vector{Int}
    bond_ids::Vector{UnitRange{Int}}
    kernels::Vector{Float64}
    correction::Correction
    chunk_handler::ChunkHandler
end

function BondSystem(body::AbstractBody, pd::PointDecomposition, chunk_id::Int)
    check_bond_system_compat(body.mat)
    bonds, n_neighbors = find_bonds(body, pd.decomp[chunk_id])
    bond_ids = find_bond_ids(n_neighbors)
    chunk_handler = get_chunk_handler(bonds, pd, chunk_id)
    localize!(bonds, chunk_handler.localizer)
    position, volume = get_pos_and_vol_chunk(body, chunk_handler.point_ids)
    correction = get_correction(body.mat, chunk_handler.n_loc_points,
                                length(chunk_handler.point_ids), length(bonds))
    kernels = find_kernels(body, chunk_handler, bonds, bond_ids)
    system = BondSystem(position, volume, bonds, n_neighbors, bond_ids, kernels, correction,
                        chunk_handler)
    return system
end

function get_system(body::AbstractBody{Material}, pd::PointDecomposition,
                    chunk_id::Int) where {Material<:AbstractBondSystemMaterial}
    return BondSystem(body, pd, chunk_id)
end

@inline function system_type(mat::AbstractBondSystemMaterial)
    return BondSystem{correction_type(mat)}
end

function check_bond_system_compat(::M) where {M<:AbstractMaterial}
    msg = "body with material `$(M)` incompatible to `BondSystem`!\n"
    msg *= "The material has to be a subtype of `AbstractBondSystemMaterial`!\n"
    return throw(ArgumentError(msg))
end

function check_bond_system_compat(::AbstractBondSystemMaterial)
    return nothing
end

function find_bonds(body::AbstractBody, loc_points::UnitRange{Int})
    δmax = maximum_horizon(body)
    nhs = GridNeighborhoodSearch{3}(search_radius=δmax, n_points=body.n_points)
    initialize_grid!(nhs, body.position)
    bonds = Vector{Bond}()
    sizehint!(bonds, body.n_points * 300)
    n_neighbors = zeros(Int, length(loc_points))
    for (li, i) in enumerate(loc_points)
        n_neighbors[li] = find_bonds!(bonds, nhs, body.position, body.fail_permit,
                                      get_point_param(body, :δ, i), i)
    end
    filter_bonds!(bonds, n_neighbors, loc_points, body)
    return bonds, n_neighbors
end

function find_bonds!(bonds::Vector{Bond}, nhs::PointNeighbors.GridNeighborhoodSearch,
                     position::Matrix{Float64}, fail_permit::Vector{Bool}, δ::Float64,
                     point_id::Int)
    n_bonds_pre = length(bonds)
    foreach_neighbor(position, position, nhs, point_id; search_radius=δ) do i, j, _, L
        if i != j
            check_point_duplicates(L, i, j)
            push!(bonds, Bond(j, L, fail_permit[i] & fail_permit[j]))
        end
    end
    n_neighbors = length(bonds) - n_bonds_pre
    return n_neighbors
end

@inline function check_point_duplicates(L::Float64, i::Int, j::Int)
    if L == 0
        msg = "point duplicate found!\n"
        msg *= "Point #$(i) has a duplicate #$(j) which will lead to `NaN`s!\n"
        error(msg)
    end
    return nothing
end

function filter_bonds!(bonds::Vector{Bond}, n_neighbors::Vector{Int},
                       loc_points::UnitRange{Int}, body::AbstractBody)
    for crack in body.point_sets_precracks
        filter_bonds_by_crack!(bonds, n_neighbors, loc_points, crack, body)
    end
    return nothing
end

function filter_bonds_by_crack!(bonds::Vector{Bond}, n_neighbors::Vector{Int},
                                loc_points::UnitRange{Int}, crack::PointSetsPreCrack,
                                body::AbstractBody)
    filter_bonds(crack) || return nothing
    set_a, set_b = body.point_sets[crack.set_a], body.point_sets[crack.set_b]
    bond_ids = find_bond_ids(n_neighbors)
    bonds_to_delete = fill(false, length(bonds))
    for (loc_point_id, point_id) in enumerate(loc_points)
        for bond_id in bond_ids[loc_point_id]
            bond = bonds[bond_id]
            neighbor_id = bond.neighbor
            point_in_a = in(point_id, set_a)
            point_in_b = in(point_id, set_b)
            neigh_in_a = in(neighbor_id, set_a)
            neigh_in_b = in(neighbor_id, set_b)
            if (point_in_a && neigh_in_b) || (point_in_b && neigh_in_a)
                bonds_to_delete[bond_id] = true
                n_neighbors[loc_point_id] -= 1
            end
        end
    end
    deleteat!(bonds, bonds_to_delete)
    return nothing
end

function find_bond_ids(n_neighbors::Vector{Int})
    bond_ids = fill(0:0, length(n_neighbors))
    bonds_start, bonds_end = 1, 0
    for i in eachindex(n_neighbors)
        bonds_end = bonds_start + n_neighbors[i] - 1
        bond_ids[i] = bonds_start:bonds_end
        bonds_start += n_neighbors[i]
    end
    return bond_ids
end

function get_pos_and_vol_chunk(body::AbstractBody, point_ids::AbstractVector{<:Integer})
    position = body.position[:, point_ids]
    volume = body.volume[point_ids]
    return position, volume
end

function get_chunk_handler(bonds::Vector{Bond}, pd::PointDecomposition, chunk_id::Int)
    loc_points = pd.decomp[chunk_id]
    n_loc_points = length(loc_points)
    halo_points = find_halo_points(bonds, loc_points)
    hidxs_by_src = sort_halo_by_src!(halo_points, pd.point_src, length(loc_points))
    point_ids = vcat(loc_points, halo_points)
    localizer = find_localizer(point_ids)
    return ChunkHandler(n_loc_points, point_ids, loc_points, halo_points, hidxs_by_src,
                        localizer)
end

function find_kernels(body::AbstractBody, chunk_handler::ChunkHandler, bonds::Vector{Bond},
                      bond_ids::Vector{UnitRange{Int}})
    hasproperty(body.mat, :kernel) || return Vector{Float64}()
    kernels = zeros(length(bonds))
    for i in each_point_idx(chunk_handler)
        params = get_point_param(body, i)
        for bond_id in bond_ids[i]
            bond = bonds[bond_id]
            kernels[bond_id] = get_kernel(body.mat, params, bond.length)
        end
    end
    return kernels
end

function get_kernel(mat::AbstractMaterial, params::AbstractPointParameters, L)
    ω = mat.kernel(params.δ, L)
    return ω
end

@inline function kernel(system::AbstractBondSystem, bond_id::Int)
    return system.kernels[bond_id]
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

@inline each_bond_idx(system::AbstractBondSystem, point_id::Int) = system.bond_ids[point_id]

function localize!(bonds::Vector{Bond}, localizer::Dict{Int,Int})
    for i in eachindex(bonds)
        bond = bonds[i]
        bonds[i] = Bond(localizer[bond.neighbor], bond.length, bond.fail_permit)
    end
    return nothing
end

function break_bonds!(storage::AbstractStorage, system::AbstractBondSystem,
                      set_a::Vector{Int}, set_b::Vector{Int})
    storage.n_active_bonds .= 0
    for point_id in each_point_idx(system)
        for bond_id in each_bond_idx(system, point_id)
            bond = system.bonds[bond_id]
            neighbor_id = bond.neighbor
            point_in_a = in(point_id, set_a)
            point_in_b = in(point_id, set_b)
            neigh_in_a = in(neighbor_id, set_a)
            neigh_in_b = in(neighbor_id, set_b)
            if (point_in_a && neigh_in_b) || (point_in_b && neigh_in_a)
                storage.bond_active[bond_id] = false
            end
            storage.n_active_bonds[point_id] += storage.bond_active[bond_id]
        end
    end
    return nothing
end

function calc_timestep_point(system::AbstractBondSystem, params::AbstractPointParameters,
                             point_id::Int)
    dtsum = 0.0
    for bond_id in each_bond_idx(system, point_id)
        bond = system.bonds[bond_id]
        dtsum += system.volume[bond.neighbor] * params.bc / bond.length
    end
    return sqrt(2 * params.rho / dtsum)
end

function calc_force_density!(chunk::AbstractBodyChunk{<:AbstractBondSystem}, t, Δt)
    (; system, mat, paramsetup, storage) = chunk
    (; dmgmodel) = mat
    storage.b_int .= 0
    storage.n_active_bonds .= 0
    for point_id in each_point_idx(chunk)
        calc_failure!(storage, system, mat, dmgmodel, paramsetup, point_id)
        calc_damage!(storage, system, mat, dmgmodel, paramsetup, point_id)
        force_density_point!(storage, system, mat, paramsetup, t, Δt, point_id)
    end
    nancheck(chunk, t)
    return nothing
end

function calc_failure!(storage::AbstractStorage, system::AbstractBondSystem,
                       mat::AbstractMaterial, dmgmodel::CriticalStretch,
                       paramsetup::AbstractParameterSetup, i)
    (; εc) = get_params(paramsetup, i)
    (; position, n_active_bonds, bond_active) = storage
    (; bonds) = system
    for bond_id in each_bond_idx(system, i)
        bond = bonds[bond_id]
        j, L = bond.neighbor, bond.length
        Δxij = get_vector_diff(position, i, j)
        l = norm(Δxij)
        ε = (l - L) / L
        if ε > εc && bond.fail_permit
            bond_active[bond_id] = false
        end
        n_active_bonds[i] += bond_active[bond_id]
    end
    return nothing
end

function calc_damage!(chunk::AbstractBodyChunk{<:AbstractBondSystem})
    (; system, mat, paramsetup, storage) = chunk
    (; dmgmodel) = mat
    for point_id in each_point_idx(chunk)
        calc_damage!(storage, system, mat, dmgmodel, paramsetup, point_id)
    end
    return nothing
end

function calc_damage!(storage::AbstractStorage, system::AbstractBondSystem,
                      mat::AbstractMaterial, dmgmodel::AbstractDamageModel,
                      paramsetup::AbstractParameterSetup, i)
    @inbounds storage.damage[i] = 1 - storage.n_active_bonds[i] / system.n_neighbors[i]
    return nothing
end

function log_system(::Type{System}, options::AbstractJobOptions,
                    dh::AbstractDataHandler) where {System<:AbstractBondSystem}
    n_bonds = calc_n_bonds(dh)
    msg = "BOND SYSTEM"
    body_name = string(get_body_name(dh))
    isempty(body_name) || (msg *= " `" * body_name * "`")
    msg *= "\n"
    msg *= msg_qty("number of bonds", n_bonds)
    log_it(options, msg)
    return nothing
end

function calc_n_bonds(dh::AbstractThreadsBodyDataHandler)
    n_bonds = 0
    for chunk in dh.chunks
        n_bonds += length(chunk.system.bonds)
    end
    return n_bonds
end

function calc_n_bonds(dh::AbstractMPIBodyDataHandler)
    n_bonds = MPI.Reduce(length(dh.chunk.system.bonds), MPI.SUM, mpi_comm())
    return n_bonds
end

function init_field_system(system::AbstractBondSystem, ::Val{:bond_active})
    return ones(Bool, get_n_bonds(system))
end

function init_field_system(system::AbstractBondSystem, ::Val{:n_active_bonds})
    return copy(system.n_neighbors)
end

function init_field_system(system::AbstractBondSystem, ::Val{:damage})
    return zeros(get_n_loc_points(system))
end

function req_point_data_fields_fracture(::Type{<:AbstractBondSystemMaterial})
    return (:damage, :n_active_bonds)
end

function req_bond_data_fields_fracture(::Type{<:AbstractBondSystemMaterial})
    return (:bond_active,)
end

function req_data_fields_fracture(::Type{<:AbstractBondSystemMaterial})
    return ()
end

function required_point_parameters(::Type{<:AbstractBondSystemMaterial})
    return (:δ, :rho, elasticity_parameters()...)
end

function get_required_point_parameters(::AbstractBondSystemMaterial, p::Dict{Symbol,Any})
    return (; get_horizon(p)..., get_density(p)..., get_elastic_params(p)...)
end

function allowed_material_kwargs(::AbstractBondSystemMaterial)
    return (discretization_kwargs()..., elasticity_kwargs()..., fracture_kwargs()...)
end

@inline get_n_bonds(system::AbstractBondSystem) = length(system.bonds)

function log_material(mat::M; indentation::Int=2) where {M<:AbstractBondSystemMaterial}
    msg = msg_qty("material type", nameof(M); indentation)
    msg *= msg_qty("correction type", correction_type(mat); indentation)
    for prop in fieldnames(M)
        msg *= log_material_property(Val(prop), mat; indentation)
    end
    return msg
end

function log_material_property(prop::Val{S}, mat::AbstractBondSystemMaterial;
                               indentation::Int=2) where {S}
    msg = msg_qty(string(prop), getfield(mat, S); indentation)
    return msg
end

function log_material_property(::Val{:dmgmodel}, mat::AbstractBondSystemMaterial;
                               indentation::Int=2)
    msg = msg_qty("damage model type", typeof(mat.dmgmodel); indentation)
    return msg
end

function log_material_property(::Val{:kernel}, mat::AbstractBondSystemMaterial;
                               indentation::Int=2)
    msg = msg_qty("kernel function", mat.kernel; indentation)
    return msg
end

function log_material(mat::M; indentation::Int=2) where {M<:AbstractCorrespondenceMaterial}
    msg = msg_qty("material type", nameof(M); indentation)
    for prop in fieldnames(M)
        msg *= log_material_property(Val(prop), mat; indentation)
    end
    return msg
end
