
"""

"""
struct BBMaterial <: AbstractMaterial end

struct BBPointParameters <: AbstractPointParameters
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

@inline point_param_type(::BBMaterial) = BBPointParameters

@inline allowed_material_kwargs(::BBMaterial) = DEFAULT_POINT_KWARGS

function get_point_params(::BBMaterial, p::Dict{Symbol,Any})
    δ = get_horizon(p)
    rho = get_density(p)
    if haskey(p, :nu)
        msg = "keyword `nu` is not allowed for BBMaterial!\n"
        msg *= "Material BBMaterial has a limitation on the possion ration!\n"
        msg *= "When using a BBMaterial, `nu` is hardcoded as 1/4.\n"
        throw(ArgumentError(msg))
    else
        p[:nu] = 0.25
    end
    E, nu, G, K, λ, μ = get_elastic_params(p)
    Gc, εc = get_frac_params(p, δ, K)
    bc = 18 * K / (π * δ^4) # bond constant
    return BBPointParameters(δ, rho, E, nu, G, K, λ, μ, Gc, εc, bc)
end

@inline discretization_type(::BBMaterial) = BondDiscretization

@inline function init_discretization(body::Body{BBMaterial}, args...)
    return init_bond_discretization(body, args...)
end

struct BBVerletStorage <: AbstractStorage
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

const BBStorage = Union{BBVerletStorage}

@inline storage_type(::BBMaterial, ::VelocityVerlet) = BBVerletStorage

function init_storage(::Body{BBMaterial}, ::VelocityVerlet, bd::BondDiscretization,
                      ch::ChunkHandler)
    n_loc_points = length(ch.loc_points)
    position = copy(bd.position)
    displacement = zeros(3, n_loc_points)
    velocity = zeros(3, n_loc_points)
    velocity_half = zeros(3, n_loc_points)
    acceleration = zeros(3, n_loc_points)
    b_int = zeros(3, n_loc_points)
    b_ext = zeros(3, n_loc_points)
    damage = zeros(n_loc_points)
    bond_active = ones(Bool, length(bd.bonds))
    n_active_bonds = copy(bd.n_neighbors)
    return BBVerletStorage(position, displacement, velocity, velocity_half, acceleration,
                           b_int, b_ext, damage, bond_active, n_active_bonds)
end

@inline reads_from_halo(::BBMaterial) = (:position,)
@inline writes_to_halo(::BBMaterial) = ()

function force_density_point!(s::BBStorage, bd::BondDiscretization,
                              param::BBPointParameters, i::Int)
    for bond_id in each_bond_idx(bd, i)
        bond = bd.bonds[bond_id]
        j, L = bond.neighbor, bond.length

        # current bond length
        Δxijx = s.position[1, j] - s.position[1, i]
        Δxijy = s.position[2, j] - s.position[2, i]
        Δxijz = s.position[3, j] - s.position[3, i]
        l = sqrt(Δxijx * Δxijx + Δxijy * Δxijy + Δxijz * Δxijz)

        # bond strain
        ε = (l - L) / L

        # failure mechanism
        if ε > param.εc && bond.fail_permit
            s.bond_active[bond_id] = false
        end
        s.n_active_bonds[i] += s.bond_active[bond_id]

        # update of force density
        temp = s.bond_active[bond_id] * param.bc * ε / l * bd.volume[j]
        s.b_int[1, i] += temp * Δxijx
        s.b_int[2, i] += temp * Δxijy
        s.b_int[3, i] += temp * Δxijz
    end
    return nothing
end
