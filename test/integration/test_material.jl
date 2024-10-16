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
    bc::Float64
end
function TestPointParameters(mat::TestMaterial, p::Dict{Symbol,Any})
    (; δ, rho, E, nu, G, K, λ, μ) = Peridynamics.get_required_point_parameters(mat, p)
    (; Gc, εc) = Peridynamics.get_frac_params(p, δ, K)
    bc = 18 * K / (π * δ^4) # bond constant
    return TestPointParameters(δ, rho, E, nu, G, K, λ, μ, Gc, εc, bc)
end
Peridynamics.@params TestMaterial TestPointParameters
struct TestStorage <: Peridynamics.AbstractStorage
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
Peridynamics.@storage TestMaterial TestStorage
Peridynamics.@loc_to_halo_fields TestStorage :position
Peridynamics.@halo_to_loc_fields TestStorage :b_int
function Peridynamics.init_field(::TestMaterial, ::Peridynamics.VelocityVerlet,
                                 system::Peridynamics.BondSystem, ::Val{:b_int})
    return zeros(3, Peridynamics.get_n_points(system))
end
function Peridynamics.force_density_point!(storage::TestStorage,
                                           system::Peridynamics.BondSystem, ::TestMaterial,
                                           params::TestPointParameters, i::Int)
    for bond_id in Peridynamics.each_bond_idx(system, i)
        bond = system.bonds[bond_id]
        j, L = bond.neighbor, bond.length
        Δxij = Peridynamics.get_coordinates_diff(storage, i, j)
        l = Peridynamics.LinearAlgebra.norm(Δxij)
        ε = (l - L) / L
        Peridynamics.stretch_based_failure!(storage, system, bond, params, ε, i, bond_id)
        b_int = Peridynamics.bond_failure(storage, bond_id) *
                Peridynamics.surface_correction_factor(system.correction, bond_id) *
                params.bc * ε / l * system.volume[j] .* Δxij
        Peridynamics.update_add_b_int!(storage, i, b_int)
    end
    return nothing
end
