"""
    CKIMaterial <: AbstractMaterial

Material type for continuum-kinematics-inspired peridynamic simulations

# Allowed material parameters

- `horizon::Float64`: Radius of point interactions
- `rho::Float64`: Density
- `E::Float64`: Young's modulus
- `nu::Float64`: Poisson's ratio
- `Gc::Float64`: Critical energy release rate
- `epsilon_c::Float64`: Critical strain

# Allowed export fields

TODO struct
"""
struct CKIMaterial <: AbstractInteractionSystemMaterial end

struct CKIPointParameters <: AbstractPointParameters
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
    C1::Float64
    C2::Float64
    C3::Float64
end

function CKIPointParameters(::CKIMaterial, p::Dict{Symbol,Any})
    δ = get_horizon(p)
    rho = get_density(p)
    E, nu, G, K, λ, μ = get_elastic_params(p)
    Gc, εc = get_frac_params(p, δ, K)
    C1, C2, C3 = cki_parameters(p, δ, λ, μ)
    return CKIPointParameters(δ, rho, E, nu, G, K, λ, μ, Gc, εc, C1, C2, C3)
end

function cki_parameters(p::Dict{Symbol,Any}, δ, λ, μ)
    if haskey(p, :C1)
        C1::Float64 = float(p[:C1])
    else
        C1 = 0.0
    end

    if haskey(p, :C2)
        C2::Float64 = float(p[:C2])
    else
        C2 = 0.0
    end

    if haskey(p, :C3)
        C3::Float64 = float(p[:C3])
    else
        C3 = 0.0
    end

    if C1 ≈ 0 && C2 ≈ 0 && C3 ≈ 0
        C1 = 30 / π * μ / δ^4
        C2 = 0.0
        C3 = 32 / π^4 * (λ - μ) / δ^12
    else
        msg = "parameters for CKIMaterial specified manually!\n"
        msg *= "Be careful when adjusting these parameters to avoid unexpected outcomes!"
        @warn msg
    end

    return C1, C2, C3
end

function allowed_material_kwargs(::CKIMaterial)
    return (DEFAULT_POINT_KWARGS..., :C1, :C2, :C3)
end

@params CKIMaterial CKIPointParameters

@inline get_c2(params::CKIPointParameters) = params.C2
@inline get_c2(params) = 0.0
@inline get_c3(params::CKIPointParameters) = params.C3
@inline get_c3(params) = 0.0

struct CKIVerletStorage <: AbstractStorage
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

function CKIVerletStorage(::CKIMaterial, ::VelocityVerlet, system::InteractionSystem, ch)
    n_loc_points = length(ch.loc_points)
    position = copy(system.position)
    displacement = zeros(3, n_loc_points)
    velocity = zeros(3, n_loc_points)
    velocity_half = zeros(3, n_loc_points)
    acceleration = zeros(3, n_loc_points)
    b_int = zeros(3, n_loc_points)
    b_ext = zeros(3, n_loc_points)
    damage = zeros(n_loc_points)
    bond_active = ones(Bool, length(system.bonds))
    n_active_bonds = copy(system.n_neighbors)
    s = CKIVerletStorage(position, displacement, velocity, velocity_half, acceleration,
                         b_int, b_ext, damage, bond_active, n_active_bonds)
    return s
end

@storage CKIMaterial VelocityVerlet CKIVerletStorage

@loc_to_halo_fields CKIVerletStorage :position

struct CKIRelaxationStorage <: AbstractStorage
    position::Matrix{Float64}
    displacement::Matrix{Float64}
    velocity::Matrix{Float64}
    velocity_half::Matrix{Float64}
    velocity_half_old::Matrix{Float64}
    b_int::Matrix{Float64}
    b_int_old::Matrix{Float64}
    b_ext::Matrix{Float64}
    density_matrix::Matrix{Float64}
    damage::Vector{Float64}
    bond_active::Vector{Bool}
    n_active_bonds::Vector{Int}
end

function CKIRelaxationStorage(::CKIMaterial, ::DynamicRelaxation, system::InteractionSystem,
                              ch)
    n_loc_points = length(ch.loc_points)
    position = copy(system.position)
    displacement = zeros(3, n_loc_points)
    velocity = zeros(3, n_loc_points)
    velocity_half = zeros(3, n_loc_points)
    velocity_half_old = zeros(3, n_loc_points)
    b_int = zeros(3, n_loc_points)
    b_int_old = zeros(3, n_loc_points)
    b_ext = zeros(3, n_loc_points)
    density_matrix = zeros(3, n_loc_points)
    damage = zeros(n_loc_points)
    bond_active = ones(Bool, length(system.bonds))
    n_active_bonds = copy(system.n_neighbors)
    s = CKIRelaxationStorage(position, displacement, velocity, velocity_half,
                             velocity_half_old, b_int, b_int_old, b_ext, density_matrix,
                             damage, bond_active, n_active_bonds)
    return s
end

@storage CKIMaterial DynamicRelaxation CKIRelaxationStorage

@loc_to_halo_fields CKIRelaxationStorage :position

const CKIStorage = Union{CKIVerletStorage,CKIRelaxationStorage}
