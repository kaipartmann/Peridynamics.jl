
"""
    BBMaterial <: AbstractMaterial

material type for bond-based peridynamic simulations

# Allowed material parameters

- `:horizon::Float64`: radius of point interactions
- `:rho::Float64`: density
- `:E::Float64`: Young's modulus
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
        msg *= "Bond-based peridynamics has a limitation on the possion ratio.\n"
        msg *= "Therefore, when using BBMaterial, `nu` is hardcoded as 1/4.\n"
        throw(ArgumentError(msg))
    else
        p[:nu] = 0.25
    end
    E, nu, G, K, λ, μ = get_elastic_params(p)
    Gc, εc = get_frac_params(p, δ, K)
    bc = 18 * K / (π * δ^4) # bond constant
    return BBPointParameters(δ, rho, E, nu, G, K, λ, μ, Gc, εc, bc)
end

@inline system_type(::BBMaterial) = BondSystem

@inline function get_system(body::AbstractBody{BBMaterial}, args...)
    return get_bond_system(body, args...)
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

function init_storage(::AbstractBody{BBMaterial}, ::VelocityVerlet, bd::BondSystem,
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

@inline get_halo_read_fields(s::BBStorage) = (s.position,)
@inline get_halo_write_fields(::BBStorage)  = ()

function force_density_point!(s::BBStorage, bd::BondSystem, ::BBMaterial,
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
