struct TestMaterial <: Peridynamics.AbstractMaterial end
struct TestPointParameters <: Peridynamics.AbstractPointParameters
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
end
Peridynamics.point_param_type(::TestMaterial) = TestPointParameters
Peridynamics.allowed_material_kwargs(::TestMaterial) = Peridynamics.DEFAULT_POINT_KWARGS
function Peridynamics.get_point_params(::TestMaterial, p::Dict{Symbol,Any})
    δ = Peridynamics.get_horizon(p)
    rho = Peridynamics.get_density(p)
    E, nu, G, K, λ, μ = Peridynamics.get_elastic_params(p)
    Gc, εc = Peridynamics.get_frac_params(p, δ, K)
    return TestPointParameters(δ, rho, E, nu, G, K, λ, μ, Gc, εc)
end
Peridynamics.discretization_type(::TestMaterial) = Peridynamics.BondDiscretization
function Peridynamics.init_discretization(body::Body{TestMaterial}, args...)
    return Peridynamics.init_bond_discretization(body, args...)
end
struct TestVerletStorage <: Peridynamics.AbstractStorage
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
function Peridynamics.storage_type(::TestMaterial, ::Peridynamics.VelocityVerlet)
    return TestVerletStorage
end
function Peridynamics.init_storage(::Body{TestMaterial}, ::Peridynamics.VelocityVerlet,
                                    bd::Peridynamics.BondDiscretization,
                                    ch::Peridynamics.ChunkHandler)
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
    return TestVerletStorage(position, displacement, velocity, velocity_half,
                                acceleration, b_int, b_ext, damage, bond_active,
                                n_active_bonds)
end
Peridynamics.get_halo_read_fields(s::TestVerletStorage) = (s.position,)
Peridynamics.get_halo_write_fields(s::TestVerletStorage)  = (s.b_int,)
