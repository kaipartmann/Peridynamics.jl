"""
all calculated parameters that should not change ever inside the time loop!
"""
struct SimParameters
    pc::PointCloud
    n_bonds::Int
    bond_data::Vector{Tuple{Int, Int, Float64, Bool}}
    n_family_members::Vector{Int}
    mat::PDMaterial
end

"""

"""
struct GlobalStorage
    b_ext::Matrix{Float64}
    damage::Matrix{Float64}
    bond_failure::Vector{Int}
end

"""

"""
struct ThreadLocalStorage
    position::Matrix{Float64}
    displacement::Matrix{Float64}
    velocity::Matrix{Float64}
    velocity_half::Matrix{Float64}
    acceleration::Matrix{Float64}
    b_int::Matrix{Float64}
end

"""

"""
struct PDProblem
    sp::SimParameters
    gs::GlobalStorage
    tls::Vector{ThreadLocalStorage}
end
