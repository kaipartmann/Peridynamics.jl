"""
TODO
"""
struct OSBMaterial <: AbstractMaterial end

struct OSBPointParameters <: AbstractPointParameters
    δ::Float64
    rho::Float64
    E::Float64
    nu::Float64
    G::Float64
    K::Float64
    λ::Float64
    μ::Float64
    Gc::Float64
    εc::Float64
    bc::Float64
end

@inline point_param_type(::OSBMaterial) = OSBPointParameters
@inline allowed_material_kwargs(::OSBMaterial) = DEFAULT_POINT_KWARGS

function get_point_params(::OSBMaterial, p::Dict{Symbol,Any})
    δ = get_horizon(p)
    rho = get_density(p)
    E, nu, G, K, λ, μ = get_elastic_params(p)
    Gc, εc = get_frac_params(p, δ, K)
    bc = 18 * K / (π * δ^4)
    return OSBPointParameters(δ, rho, E, nu, G, K, λ, μ, Gc, εc, bc)
end

@inline discretization_type(::OSBMaterial) = BondSystem

@inline function init_discretization(body::Body{OSBMaterial}, args...)
    return init_bond_discretization(body, args...)
end

struct OSBVerletStorage <: AbstractStorage
    position::Matrix{Float64}
    displacement::Matrix{Float64}
    velocity::Matrix{Float64}
    velocity_half::Matrix{Float64}
    acceleration::Matrix{Float64}
    b_int::Matrix{Float64}
    b_ext::Matrix{Float64}
    damage::Vector{Float64}
    bond_active::Vector{Bool}
    n_active_bonds::Vector{Int}
end

const OSBStorage = Union{OSBVerletStorage}

@inline storage_type(::OSBMaterial, ::VelocityVerlet) = OSBVerletStorage

function init_storage(::Body{OSBMaterial}, ::VelocityVerlet, bd::BondSystem,
                      ch::ChunkHandler)
    n_loc_points = length(ch.loc_points)
    position = copy(bd.position)
    displacement = zeros(3, n_loc_points)
    velocity = zeros(3, n_loc_points)
    velocity_half = zeros(3, n_loc_points)
    acceleration = zeros(3, n_loc_points)
    b_int = zeros(3, length(ch.point_ids))
    b_ext = zeros(3, n_loc_points)
    damage = zeros(n_loc_points)
    bond_active = ones(Bool, length(bd.bonds))
    n_active_bonds = copy(bd.n_neighbors)
    return OSBVerletStorage(position, displacement, velocity, velocity_half, acceleration,
                            b_int, b_ext, damage, bond_active, n_active_bonds)
end

@inline get_halo_read_fields(s::OSBStorage) = (s.position,)
@inline get_halo_write_fields(s::OSBStorage)  = (s.b_int,)

function force_density_point!(s::OSBStorage, bd::BondSystem, mat::OSBMaterial,
                              param::OSBPointParameters, i::Int)
    # weighted volume
    wvol = calc_weighted_volume(s, bd, param, i)
    iszero(wvol) && return nothing
    # dilatation
    dil = calc_dilatation(s, bd, param, wvol, i)
    # force density
    c1 = 15.0 * param.G / wvol
    c2 = dil * (3.0 * param.K / wvol - c1 / 3.0)
    for bond_id in each_bond_idx(bd, i)
        bond = bd.bonds[bond_id]
        j, L = bond.neighbor, bond.length
        Δxijx = s.position[1, j] - s.position[1, i]
        Δxijy = s.position[2, j] - s.position[2, i]
        Δxijz = s.position[3, j] - s.position[3, i]
        l = sqrt(Δxijx * Δxijx + Δxijy * Δxijy + Δxijz * Δxijz)
        ε = (l - L) / L

        # failure mechanism
        if ε > param.εc && bond.fail_permit
            s.bond_active[bond_id] = false
        end
        s.n_active_bonds[i] += s.bond_active[bond_id]

        # update of force density
        ωij = (1 + param.δ / L) * s.bond_active[bond_id]
        temp = ωij * (c2 * L + c1 * (l - L)) / l
        s.b_int[1, i] += temp * Δxijx * bd.volume[j]
        s.b_int[2, i] += temp * Δxijy * bd.volume[j]
        s.b_int[3, i] += temp * Δxijz * bd.volume[j]
        s.b_int[1, j] -= temp * Δxijx * bd.volume[i]
        s.b_int[2, j] -= temp * Δxijy * bd.volume[i]
        s.b_int[3, j] -= temp * Δxijz * bd.volume[i]
    end
    return nothing
end

function calc_weighted_volume(s::OSBStorage, bd::BondSystem,
                              param::OSBPointParameters, i::Int)
    wvol = 0.0
    for bond_id in each_bond_idx(bd, i)
        bond = bd.bonds[bond_id]
        j, L = bond.neighbor, bond.length
        ΔXijx = bd.position[1, j] - bd.position[1, i]
        ΔXijy = bd.position[2, j] - bd.position[2, i]
        ΔXijz = bd.position[3, j] - bd.position[3, i]
        ΔXij_sq = ΔXijx * ΔXijx + ΔXijy * ΔXijy + ΔXijz * ΔXijz
        ωij = (1 + param.δ / L) * s.bond_active[bond_id]
        wvol += ωij * ΔXij_sq * bd.volume[j]
    end
    return wvol
end

function calc_dilatation(s::OSBStorage, bd::BondSystem, param::OSBPointParameters,
                         wvol::Float64, i::Int)
    dil = 0.0
    c1 = 3.0 / wvol
    for bond_id in each_bond_idx(bd, i)
        bond = bd.bonds[bond_id]
        j, L = bond.neighbor, bond.length
        Δxijx = s.position[1, j] - s.position[1, i]
        Δxijy = s.position[2, j] - s.position[2, i]
        Δxijz = s.position[3, j] - s.position[3, i]
        l = sqrt(Δxijx * Δxijx + Δxijy * Δxijy + Δxijz * Δxijz)
        ωij = (1 + param.δ / L) * s.bond_active[bond_id]
        dil += ωij * c1 * L * (l - L) * bd.volume[j]
    end
    return dil
end
