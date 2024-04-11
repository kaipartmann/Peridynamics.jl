"""
    OSBMaterial <: AbstractMaterial

material type for ordinary state-based peridynamic simulations

# Allowed material parameters

- `:horizon::Float64`: radius of point interactions
- `:rho::Float64`: density
- `:E::Float64`: Young's modulus
- `:nu::Float64`: Poisson's ratio
- `:Gc::Float64`: critical energy release rate
- `:epsilon_c::Float64`: critical strain

# Allowed export fields

- `position::Matrix{Float64}`: position of each point
- `displacement::Matrix{Float64}`: displacement of each point
- `velocity::Matrix{Float64}`: velocity of each point
- `velocity_half::Matrix{Float64}`: velocity parameter for Verlet time solver
- `acceleration::Matrix{Float64}`: acceleration of each point
- `b_int::Matrix{Float64}`: internal force density of each point
- `b_ext::Matrix{Float64}`: external force density of each point
- `damage::Vector{Float64}`: damage of each point
- `n_active_bonds::Vector{Int}`: number of intact bonds for each point
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

function OSBPointParameters(::OSBMaterial, p::Dict{Symbol,Any})
    δ = get_horizon(p)
    rho = get_density(p)
    E, nu, G, K, λ, μ = get_elastic_params(p)
    Gc, εc = get_frac_params(p, δ, K)
    bc = 18 * K / (π * δ^4)
    return OSBPointParameters(δ, rho, E, nu, G, K, λ, μ, Gc, εc, bc)
end

@pointparams OSBMaterial OSBPointParameters

@system OSBMaterial BondSystem

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

function OSBVerletStorage(::AbstractBody{OSBMaterial}, ::VelocityVerlet, bd::BondSystem,
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
    s = OSBVerletStorage(position, displacement, velocity, velocity_half, acceleration,
                         b_int, b_ext, damage, bond_active, n_active_bonds)
    return s
end

@storage OSBMaterial VelocityVerlet OSBVerletStorage

@inline point_data_fields(::Type{OSBVerletStorage}) = (:position, :displacement, :velocity,
                                                       :velocity_half, :acceleration,
                                                       :b_int, :b_ext, :damage,
                                                       :n_active_bonds)

@inline halo_read_fields(::OSBVerletStorage) = (:position,)
@inline halo_write_fields(::OSBVerletStorage)  = (:b_int,)
@inline is_halo_field(::OSBVerletStorage, ::Val{:position}) = true
@inline is_halo_field(::OSBVerletStorage, ::Val{:b_int}) = true
@inline is_halo_field(::OSBVerletStorage, ::Val{F}) where {F} = false

function force_density_point!(s::OSBVerletStorage, bd::BondSystem, ::OSBMaterial,
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

function calc_weighted_volume(s::OSBVerletStorage, bd::BondSystem,
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

function calc_dilatation(s::OSBVerletStorage, bd::BondSystem, param::OSBPointParameters,
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
