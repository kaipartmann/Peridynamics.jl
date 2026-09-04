"""
    BondSystem{Correction}

$(extension_api_note())

A type for a system for all peridynamic formulations that work with just bonds of two
points.

# Type Parameters

- `Correction<:AbstractCorrection`: Applied surface correction.
- `N`: Number of spatial dimensions, see [`get_n_dim`](@ref).
- `FT`: Float type of the simulation data, see [`float_type`](@ref).

Every array field is a type parameter as well, so that a whole system can be moved to another
array backend with `Adapt.adapt`, see [`@system`](@ref) and `host_system_type`. Dispatch on
this system therefore has to be written `BondSystem{<:EnergySurfaceCorrection}` and never
`BondSystem{EnergySurfaceCorrection}`.

See [`AbstractBondSystem`](@ref) for what a bond is on this system and how a kernel reads
one.

# Fields

$(block_table(BondSystem))

`correction` is the applied surface correction and `kernels` holds the value of the influence
function of every bond, see [`kernel`](@ref).
"""
@system struct BondSystem{Correction<:AbstractCorrection} <: AbstractBondSystem
    position::PointVector{Float64}
    volume::PointScalar
    neighbor::BondScalar{Int}
    bond_length::BondScalar
    fail_permit::BondScalar{Bool}
    n_neighbors::PointScalar{Int}
    bond_ids::PointScalar{UnitRange{Int}}
    kernels::BondScalar
    correction::Correction
end

function BondSystem(body::AbstractBody, pd::PointDecomposition, chunk_id::Int)
    check_system_compat(BondSystem, body.mat)
    neighbor, bond_length, fail_permit, n_neighbors, bond_ids, chunk_handler = get_bond_data(body, pd, chunk_id)
    position, volume = get_pos_and_vol_chunk(body, chunk_handler.point_ids)
    N, FT = size(position, 1), default_float_type()
    sizes = SystemSizes{N,FT}(chunk_handler, length(neighbor))
    correction = get_correction(body.mat, sizes)
    kernels = find_kernels(body, chunk_handler, neighbor, bond_length, bond_ids)
    system = BondSystem{typeof(correction),N,FT}(position, volume, neighbor, bond_length,
                                                 fail_permit, n_neighbors, bond_ids,
                                                 kernels, correction, chunk_handler)
    return system
end

function get_bond_data(body::AbstractBody, pd::PointDecomposition, chunk_id)
    loc_points = pd.decomp[chunk_id]
    neighbor, bond_length, fail_permit, n_neighbors = find_bonds(body, loc_points)
    halo_points = find_halo_points(neighbor, loc_points)
    chunk_handler = ChunkHandler(pd, halo_points, chunk_id)
    localize!(neighbor, chunk_handler.localizer)
    bond_ids = find_bond_ids(n_neighbors)
    return neighbor, bond_length, fail_permit, n_neighbors, bond_ids, chunk_handler
end

function get_system(body::AbstractBody{Material}, pd::PointDecomposition,
                    chunk_id::Int) where {Material<:AbstractBondSystemMaterial}
    return BondSystem(body, pd, chunk_id)
end

@inline function system_type(mat::AbstractBondSystemMaterial,
                             ::Type{FT}=default_float_type(),
                             ::Val{N}=Val(3)) where {FT,N}
    return host_system_type(BondSystem, host_type(correction_type(mat), Val(N), FT), Val(N),
                            FT)
end

function check_system_compat(::Type{S}, ::M) where {S<:BondSystem,M<:AbstractMaterial}
    msg = "body with material `$(M)` incompatible to `BondSystem`!\n"
    msg *= "The material has to be a subtype of `AbstractBondSystemMaterial`!\n"
    return throw(ArgumentError(msg))
end

function check_system_compat(::Type{<:BondSystem}, ::AbstractBondSystemMaterial)
    return nothing
end

function find_bonds(body::AbstractBody, loc_points::AbstractVector{Int})
    δmax = maximum_horizon(body)
    nhs = GridNeighborhoodSearch{3}(search_radius=δmax, n_points=body.n_points)
    initialize_grid!(nhs, body.position)
    neighbor = Vector{Int}()
    bond_length = Vector{Float64}()
    fail_permit = Vector{Bool}()
    sizehint!(neighbor, body.n_points * 300)
    sizehint!(bond_length, body.n_points * 300)
    sizehint!(fail_permit, body.n_points * 300)
    n_neighbors = zeros(Int, length(loc_points))
    for (li, i) in enumerate(loc_points)
        n_neighbors[li] = find_bonds!(neighbor, bond_length, fail_permit, nhs, body.position,
                                      body.fail_permit, get_point_param(body, :δ, i), i)
    end
    filter_bonds!(neighbor, bond_length, fail_permit, n_neighbors, loc_points, body)
    return neighbor, bond_length, fail_permit, n_neighbors
end

function find_bonds!(neighbor::Vector{Int}, bond_length::Vector{Float64},
                     fail_permit::Vector{Bool}, nhs::PointNeighbors.GridNeighborhoodSearch,
                     position::Matrix{Float64}, body_fail_permit::Vector{Bool}, δ::Float64,
                     point_id::Int)
    n_bonds_pre = length(neighbor)
    foreach_neighbor(position, position, nhs, point_id; search_radius=δ) do i, j, _, L
        if i != j
            check_point_duplicates(L, i, j)
            push!(neighbor, j)
            push!(bond_length, L)
            push!(fail_permit, body_fail_permit[i] & body_fail_permit[j])
        end
    end
    n_neighbors = length(neighbor) - n_bonds_pre
    return n_neighbors
end

@inline function check_point_duplicates(L::Float64, i::Int, j::Int)
    if L < eps()
        msg = "point duplicate found!\n"
        msg *= "Point #$(i) has a duplicate #$(j) which will lead to `NaN`s!\n"
        error(msg)
    end
    return nothing
end

function filter_bonds!(neighbor::Vector{Int}, bond_length::Vector{Float64},
                       fail_permit::Vector{Bool}, n_neighbors::Vector{Int},
                       loc_points::AbstractVector{Int}, body::AbstractBody)
    for crack in body.point_sets_precracks
        filter_bonds_by_crack!(neighbor, bond_length, fail_permit, n_neighbors, loc_points,
                               crack, body)
    end
    return nothing
end

function filter_bonds_by_crack!(neighbor::Vector{Int}, bond_length::Vector{Float64},
                                fail_permit::Vector{Bool}, n_neighbors::Vector{Int},
                                loc_points::AbstractVector{Int}, crack::PointSetsPreCrack,
                                body::AbstractBody)
    filter_bonds(crack) || return nothing
    set_a, set_b = body.point_sets[crack.set_a], body.point_sets[crack.set_b]
    bond_ids = find_bond_ids(n_neighbors)
    bonds_to_delete = fill(false, length(neighbor))
    for (loc_point_id, point_id) in enumerate(loc_points)
        for bond_id in bond_ids[loc_point_id]
            neighbor_id = neighbor[bond_id]
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
    deleteat!(neighbor, bonds_to_delete)
    deleteat!(bond_length, bonds_to_delete)
    deleteat!(fail_permit, bonds_to_delete)
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

function find_kernels(body::AbstractBody, chunk_handler::ChunkHandler, neighbor::Vector{Int},
                      bond_length::Vector{Float64}, bond_ids::Vector{UnitRange{Int}})
    hasproperty(body.mat, :kernel) || return Vector{Float64}()
    kernels = zeros(length(neighbor))
    for i in each_point_idx(chunk_handler)
        params = get_point_param(body, i)
        for bond_id in bond_ids[i]
            kernels[bond_id] = get_kernel(body.mat, params, bond_length[bond_id])
        end
    end
    return kernels
end

function get_kernel(mat::AbstractMaterial, params::AbstractPointParameters, L)
    ω = mat.kernel(params.δ, L)
    return ω
end

"""
    kernel(system, bond_id)

$(extension_api_note())

Return the value of the influence function ``\\omega`` of bond `bond_id`. The kernel is
evaluated once during setup from the kernel function of the material, e.g. `linear_kernel`
or `cubic_b_spline_kernel`, and the initial bond length, so a material only has to look it
up.

# Example

```julia
for bond_id in Peridynamics.each_bond_idx(system, i)
    ωij = Peridynamics.kernel(system, bond_id)
end
```
"""
@inline function kernel(system::AbstractBondSystem, bond_id::Int)
    return system.kernels[bond_id]
end

function find_halo_points(neighbor::Vector{Int}, loc_points::AbstractVector{Int})
    halo_points = Vector{Int}()
    for j in neighbor
        if !in(j, loc_points) && !in(j, halo_points)
            push!(halo_points, j)
        end
    end
    return halo_points
end

function calc_timestep_point(system::AbstractBondSystem, params::AbstractPointParameters,
                             point_id::Int)
    dtsum = 0.0
    for bond_id in each_bond_idx(system, point_id)
        j = get_neighbor(system, bond_id)
        L = reference_bond_length(system, bond_id)
        dtsum += system.volume[j] * params.bc / L
    end
    return sqrt(2 * params.rho / dtsum)
end

function calc_force_density!(chunk::AbstractBodyChunk{<:AbstractBondSystem}, t, Δt)
    (; system, mat, paramsetup, storage) = chunk
    calc_force_density!(storage, system, mat, paramsetup, t, Δt)
    return nothing
end

"""
    update_bond_lengths!(storage, system, i)

$(extension_api_note())

Write the current length of every bond of point `i` into `storage.bond_length`, the cache a
material opts into by inheriting [`BondLengthCache`](@ref). Inside a simulation the package
fills the cache at the top of the point loop, before [`calc_failure!`](@ref) and before the
force density of the material, so a material and a damage model never have to call this.

Call it yourself in two situations:

- A unit test that calls a [`calc_failure!`](@ref) or a [`force_density_point!`](@ref) of
  your own on a chunk directly, because there the point loop of the package is not what
  runs.
- An entry point of your own that walks bonds outside of `calc_force_density!`, e.g. a
  `strain_energy_density_point!`, which [`export_field`](@ref) evaluates outside of a time
  step.

# Default

A no-op for a storage that does not declare `bond_length`, decided at compile time.

See also [`current_bond_length`](@ref), [`bond_stretch`](@ref), [`BondLengthCache`](@ref).
"""
@inline function update_bond_lengths!(storage::AbstractStorage, system::AbstractSystem, i)
    hasfield(typeof(storage), :bond_length) || return nothing
    (; position, bond_length) = storage
    for bond_id in each_bond_idx(system, i)
        j = get_neighbor(system, bond_id)
        Δxij = get_vector_diff(position, i, j, dims(system))
        @inbounds bond_length[bond_id] = norm(Δxij)
    end
    return nothing
end

"""
    current_bond_length(storage, system, i, bond_id)

$(extension_api_note())

The current length of bond `bond_id` of point `i`, i.e. the distance of its two points in the
deformed configuration. [`reference_bond_length`](@ref) is the length of the same bond in the
reference configuration, and the ratio of the two is [`bond_stretch`](@ref).

This is what a damage model and a force density use instead of gathering the two positions
and taking the norm themselves. A material that inherits [`BondLengthCache`](@ref) has the
length cached and this reads the cache, a material without the field gets the distance
computed, decided at compile time. It works on every bond system and on the
[`InteractionSystem`](@ref), where it is the current length of a one-neighbor interaction.

See also [`bond_stretch`](@ref), [`calc_failure!`](@ref), [`update_bond_lengths!`](@ref),
[`each_bond_idx`](@ref).
"""
@inline function current_bond_length(storage::AbstractStorage, system::AbstractSystem, i,
                                     bond_id)
    if hasfield(typeof(storage), :bond_length)
        return @inbounds storage.bond_length[bond_id]
    end
    j = get_neighbor(system, bond_id)
    return norm(get_vector_diff(storage.position, i, j, dims(system)))
end

"""
    bond_stretch(storage, system, i, bond_id)

$(extension_api_note())

The stretch of bond `bond_id` of point `i`, i.e. `(l - L) / L` with the current length `l` of
the bond and its length `L` in the reference configuration. It is what a damage criterion
compares against the critical stretch `εc` of the point parameters, and the strain that the
force density of a bond-based material multiplies with the bond constant.

The current length comes from [`current_bond_length`](@ref), so this reads the cache of a
material that keeps one and computes the distance for a material that does not. It works on
every bond system and on the [`InteractionSystem`](@ref), where it is the stretch of a
one-neighbor interaction.

Call this when the stretch is all you need, which is the case for a failure criterion and
for the strain energy density of a bond-based material. A force density that needs the
current length as well reads that once with [`current_bond_length`](@ref) and forms the
stretch from it and [`reference_bond_length`](@ref).

See also [`current_bond_length`](@ref), [`calc_failure!`](@ref), [`each_bond_idx`](@ref).
"""
@inline function bond_stretch(storage::AbstractStorage, system::AbstractSystem, i, bond_id)
    L = reference_bond_length(system, bond_id)
    return (current_bond_length(storage, system, i, bond_id) - L) / L
end

function calc_force_density!(storage::AbstractStorage, system::AbstractBondSystem,
                             mat::AbstractBondSystemMaterial,
                             paramsetup::AbstractParameterSetup, t, Δt)
    (; dmgmodel) = mat
    storage.b_int .= 0.0
    for i in each_point_idx(system)
        update_bond_lengths!(storage, system, i)
        calc_failure!(storage, system, mat, dmgmodel, paramsetup, t, Δt, i)
        calc_damage!(storage, system, mat, dmgmodel, paramsetup, i)
        force_density_point!(storage, system, mat, paramsetup, t, Δt, i)
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
        n_bonds += get_n_bonds(chunk.system)
    end
    return n_bonds
end

function calc_n_bonds(dh::AbstractMPIBodyDataHandler)
    n_bonds = MPI.Reduce(get_n_bonds(dh.chunk.system), MPI.SUM, mpi_comm())
    return n_bonds
end

# the neighbor count is system knowledge, so this initial value serves the field inside a
# damage state and as a flat storage field alike; `bond_active` and `damage` need no hook,
# their initial values follow from the declarations of `BondFracFields`
function init_field_system(system::AbstractBondSystem, ::Val{:n_active_bonds})
    return copy(system.n_neighbors)
end

function required_point_parameters(::Type{<:AbstractBondSystemMaterial})
    return (:δ, :rho, elasticity_parameters()...)
end

"""
    BondLengthCache

$(extension_api_note())

The bond length cache of a bond system, see [`@storage_fields`](@ref). `bond_length` holds
the current length of every bond, i.e. the distance of the two points of the bond in the
deformed configuration, and the package fills it once per point and per time step before
the damage model and the force density of the material run, see
[`update_bond_lengths!`](@ref).

Inheriting this block is a decision of the material, and of the shipped ones
[`BBMaterial`](@ref), [`DHBBMaterial`](@ref), [`GBBMaterial`](@ref) and
[`OSBMaterial`](@ref) do. Nothing reads the field directly. A damage model and a force
density go through [`current_bond_length`](@ref) and [`bond_stretch`](@ref), which is what
makes them run on a material with the cache and on one without it.

$(block_table(BondLengthCache))
"""
@storage_fields BondLengthCache begin
    bond_length::BondScalar
end

function log_material(mat::M; indentation::Int=2) where {M<:AbstractBondSystemMaterial}
    msg = msg_qty("material type", nameof(M); indentation)
    if !(correction_type(mat) <: Nothing)
        msg *= msg_qty("correction type", correction_type(mat); indentation)
    end
    for prop in fieldnames(M)
        msg *= log_material_property(Val(prop), mat; indentation)
    end
    return msg
end

function log_material_property(::Val{:dmgmodel}, mat::AbstractBondSystemMaterial;
                               indentation::Int=2)
    return log_dmgmodel(mat.dmgmodel; indentation)
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
