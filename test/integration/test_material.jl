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
function TestPointParameters(mat::TestMaterial, p::Dict{Symbol,Any})
    (; δ, rho, E, nu, G, K, λ, μ) = Peridynamics.get_required_point_parameters(mat, p)
    (; Gc, εc) = Peridynamics.get_frac_params(p, δ, K)
    return TestPointParameters(δ, rho, E, nu, G, K, λ, μ, Gc, εc)
end
Peridynamics.@params TestMaterial TestPointParameters
struct TestVerletStorage <: Peridynamics.AbstractStorage
    position::Matrix{Float64}
    displacement::Matrix{Float64}
    velocity::Matrix{Float64}
    velocity_half::Matrix{Float64}
    velocity_half_old::Matrix{Float64}
    acceleration::Matrix{Float64}
    b_int::Matrix{Float64}
    b_int_old::Matrix{Float64}
    b_ext::Matrix{Float64}
    density_matrix::Matrix{Float64}
    damage::Vector{Float64}
    bond_active::Vector{Bool}
    n_active_bonds::Vector{Int}
end
Peridynamics.@storage TestMaterial TestVerletStorage
Peridynamics.@loc_to_halo_fields TestVerletStorage :position
Peridynamics.@halo_to_loc_fields TestVerletStorage :b_int
