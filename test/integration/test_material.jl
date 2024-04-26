struct TestMaterial <: Peridynamics.AbstractBondSystemMaterial{Peridynamics.NoCorrection} end
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
function TestPointParameters(::TestMaterial, p::Dict{Symbol,Any})
    δ = Peridynamics.get_horizon(p)
    rho = Peridynamics.get_density(p)
    E, nu, G, K, λ, μ = Peridynamics.get_elastic_params(p)
    Gc, εc = Peridynamics.get_frac_params(p, δ, K)
    return TestPointParameters(δ, rho, E, nu, G, K, λ, μ, Gc, εc)
end
Peridynamics.@params TestMaterial TestPointParameters
# Peridynamics.@system TestMaterial Peridynamics.BondSystem
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
function TestVerletStorage(::TestMaterial, ::Peridynamics.VelocityVerlet,
                           system::Peridynamics.BondSystem, ch::Peridynamics.ChunkHandler)
    n_loc_points = length(ch.loc_points)
    position = copy(system.position)
    displacement = zeros(3, n_loc_points)
    velocity = zeros(3, n_loc_points)
    velocity_half = zeros(3, n_loc_points)
    acceleration = zeros(3, n_loc_points)
    b_int = zeros(3, length(ch.point_ids))
    b_ext = zeros(3, n_loc_points)
    damage = zeros(n_loc_points)
    bond_active = ones(Bool, length(system.bonds))
    n_active_bonds = copy(system.n_neighbors)
    return TestVerletStorage(position, displacement, velocity, velocity_half,
                             acceleration, b_int, b_ext, damage, bond_active,
                             n_active_bonds)
end
Peridynamics.@storage TestMaterial VelocityVerlet TestVerletStorage
Peridynamics.@loc_to_halo_fields TestVerletStorage :position
Peridynamics.@halo_to_loc_fields TestVerletStorage :b_int
