"""
    CCMaterial(; maxdmg, maxjacobi, corr)

A material type used to assign the material of a [`Body`](@ref) with the local continuum
consistent (correspondence) formulation of non-ordinary state-based peridynamics.

# Keywords
- `maxdmg::Float64`: Maximum value of damage a point is allowed to obtain. If this value is
    exceeded, all bonds of that point are broken because the deformation gradient would then
    possibly contain `NaN` values.
    (default: `0.85`)
- `corr::Float64`: Correction factor used for zero-energy mode stabilization. The
    stabilization algorithm of Silling (2017) is used.
    (default: `100.0`)

!!! note "Stability of fracture simulations"
    This formulation is known to be not suitable for fracture simultations without
    stabilization of the zero-energy modes. Therefore be careful when doing fracture
    simulations and try out different parameters for `maxdmg` and `corr`.

# Examples

```julia-repl
julia> mat = CCMaterial()
CCMaterial(maxdmg=0.85, corr=100.0)
```

---

```julia
CCMaterial
```

Material type for the local continuum consistent (correspondence) formulation of
non-ordinary state-based peridynamics.

# Fields
- `maxdmg::Float64`: Maximum value of damage a point is allowed to obtain. See the
    constructor docs for more informations.
- `corr::Float64`: Correction factor used for zero-energy mode stabilization. See the
    constructor docs for more informations.

# Allowed material parameters
When using [`material!`](@ref) on a [`Body`](@ref) with `CCMaterial`, then the following
parameters are allowed:
- `horizon::Float64`: Radius of point interactions
- `rho::Float64`: Density
- `E::Float64`: Young's modulus
- `nu::Float64`: Poisson's ratio
- `Gc::Float64`: Critical energy release rate
- `epsilon_c::Float64`: Critical strain

# Allowed export fields
When specifying the `fields` keyword of [`Job`](@ref) for a [`Body`](@ref) with
`CCMaterial`, the following fields are allowed:
- `position::Matrix{Float64}`: Position of each point
- `displacement::Matrix{Float64}`: Displacement of each point
- `velocity::Matrix{Float64}`: Velocity of each point
- `velocity_half::Matrix{Float64}`: Velocity parameter for Verlet time solver
- `acceleration::Matrix{Float64}`: Acceleration of each point
- `b_int::Matrix{Float64}`: Internal force density of each point
- `b_ext::Matrix{Float64}`: External force density of each point
- `damage::Vector{Float64}`: Damage of each point
- `n_active_bonds::Vector{Int}`: Number of intact bonds of each point
"""
struct CCMaterial{CM,SI} <: AbstractCorrespondenceMaterial{CM,SI,NoCorrection}
    constitutive_model::CM
    stress_integration::SI
    maxdmg::Float64
    corr::Float64
    function CCMaterial(cm::CM, si::SI, maxdmg::Real, corr::Real) where {CM,SI}
        return new{CM,SI}(cm, si, maxdmg, corr)
    end
end

function CCMaterial(constitutive_model::AbstractConstitutiveModel=NeoHookeNonlinear();
                    stress_integration::AbstractStressIntegration=FlanaganTaylorRotation(),
                    maxdmg::Real=0.85, corr::Real=100.0)
    return CCMaterial(constitutive_model, stress_integration, maxdmg, corr)
end

function Base.show(io::IO, @nospecialize(mat::CCMaterial))
    print(io, typeof(mat))
    print(io, msg_fields_in_brackets(mat, (:maxdmg, :corr)))
    return nothing
end

struct CCPointParameters <: AbstractPointParameters
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
    C::MArray{NTuple{4,3},Float64,4,81}
end

function CCPointParameters(::CCMaterial, p::Dict{Symbol,Any})
    δ = get_horizon(p)
    rho = get_density(p)
    E, nu, G, K, λ, μ = get_elastic_params(p)
    Gc, εc = get_frac_params(p, δ, K)
    bc = 18 * K / (π * δ^4) # bond constant
    C = get_hooke_matrix(nu, λ, μ)
    return CCPointParameters(δ, rho, E, nu, G, K, λ, μ, Gc, εc, bc, C)
end

@params CCMaterial CCPointParameters

struct CCVerletStorage <: AbstractStorage
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
    rotation::Matrix{Float64}
    left_stretch::Matrix{Float64}
    stress::Matrix{Float64}
    C_1::MMatrix{3,3,Float64,9}
    C_rotated::MArray{NTuple{4,3},Float64,4,81}
end

function CCVerletStorage(mat::CCMaterial, ::VelocityVerlet, system::BondSystem, ch)
    n_loc_points = length(ch.loc_points)
    n_all_points = length(ch.point_ids)
    position = copy(system.position)
    displacement = zeros(3, n_loc_points)
    velocity = zeros(3, n_all_points)
    velocity_half = zeros(3, n_loc_points)
    acceleration = zeros(3, n_loc_points)
    b_int = zeros(3, n_all_points)
    b_ext = zeros(3, n_loc_points)
    damage = zeros(n_loc_points)
    bond_active = ones(Bool, length(system.bonds))
    n_active_bonds = copy(system.n_neighbors)
    rotation = init_rotation(mat, n_loc_points)
    left_stretch = init_left_stretch(mat, n_loc_points)
    stress = zeros(9, n_loc_points)
    C_1 = @MMatrix zeros(3, 3)
    C_rotated = @MArray zeros(3, 3, 3, 3)
    s = CCVerletStorage(position, displacement, velocity, velocity_half, acceleration,
                        b_int, b_ext, damage, bond_active, n_active_bonds, rotation,
                        left_stretch, stress, C_1, C_rotated)
    return s
end

@storage CCMaterial VelocityVerlet CCVerletStorage

@loc_to_halo_fields CCVerletStorage :position :velocity
@halo_to_loc_fields CCVerletStorage :b_int

struct CCRelaxationStorage <: AbstractStorage
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
    rotation::Matrix{Float64}
    left_stretch::Matrix{Float64}
    stress::Matrix{Float64}
    C_1::MMatrix{3,3,Float64,9}
    C_rotated::MArray{NTuple{4,3},Float64,4,81}
end

function CCRelaxationStorage(mat::CCMaterial, ::DynamicRelaxation, system::BondSystem, ch)
    n_loc_points = length(ch.loc_points)
    n_all_points = length(ch.point_ids)
    position = copy(system.position)
    displacement = zeros(3, n_loc_points)
    velocity = zeros(3, n_all_points)
    velocity_half = zeros(3, n_loc_points)
    velocity_half_old = zeros(3, n_loc_points)
    b_int = zeros(3, n_all_points)
    b_int_old = zeros(3, n_loc_points)
    b_ext = zeros(3, n_loc_points)
    density_matrix = zeros(3, n_loc_points)
    damage = zeros(n_loc_points)
    bond_active = ones(Bool, length(system.bonds))
    n_active_bonds = copy(system.n_neighbors)
    rotation = init_rotation(mat, n_loc_points)
    left_stretch = init_left_stretch(mat, n_loc_points)
    stress = zeros(9, n_loc_points)
    C_1 = @MMatrix zeros(3, 3)
    C_rotated = @MArray zeros(3, 3, 3, 3)
    s = CCRelaxationStorage(position, displacement, velocity, velocity_half,
                            velocity_half_old, b_int, b_int_old, b_ext, density_matrix,
                            damage, bond_active, n_active_bonds, rotation, left_stretch,
                            stress, C_1, C_rotated)
    return s
end

@storage CCMaterial DynamicRelaxation CCRelaxationStorage

@loc_to_halo_fields CCRelaxationStorage :position :velocity
@halo_to_loc_fields CCRelaxationStorage :b_int

const CCStorage = Union{CCVerletStorage,CCRelaxationStorage}

function init_rotation(::CCMaterial{CM,<:FlanaganTaylorRotation}, n::Int) where{CM}
    R = zeros(9, n)
    R[[1, 5, 9], :] .= 1.0
    return R
end

function init_rotation(::CCMaterial{CM,<:NoRotation}, n::Int) where {CM}
    return Array{Float64,2}(undef, 0, 0)
end

function init_left_stretch(::CCMaterial{CM,<:FlanaganTaylorRotation}, n::Int) where {CM}
    V = zeros(9, n)
    V[[1, 5, 9], :] .= 1.0
    return V
end

function init_left_stretch(::CCMaterial{CM,<:NoRotation}, n::Int) where {CM}
    return Array{Float64,2}(undef, 0, 0)
end

function point_data_fields(::Type{S}) where {S<:CCStorage}
    fields = (:position, :displacement, :velocity, :velocity_half, :acceleration, :b_int,
              :b_ext, :damage, :n_active_bonds, :rotation, :left_stretch, :stress)
    return fields
end

point_data_field(s::CCStorage, ::Val{:left_stretch}) = getfield(s, :left_stretch)
point_data_field(s::CCStorage, ::Val{:rotation}) = getfield(s, :rotation)
point_data_field(s::CCStorage, ::Val{:stress}) = getfield(s, :stress)

function force_density_point!(storage::CCStorage, system::BondSystem, mat::CCMaterial,
                              paramhandler::AbstractParameterHandler, t::Float64,
                              Δt::Float64, i::Int)
    params = get_params(paramhandler, i)
    force_density_point!(storage, system, mat, params, t, Δt, i)
    return nothing
end

function force_density_point!(storage::CCStorage, system::BondSystem, mat::CCMaterial,
                              params::CCPointParameters, t::Float64, Δt::Float64,
                              i::Int)
    (; bonds, volume) = system
    (; bond_active) = storage
    F, Ḟ, Kinv = calc_deformation_gradient(storage, system, mat, params, i)
    if storage.damage[i] > mat.maxdmg || containsnan(F)
        kill_point!(storage, system, i)
        return nothing
    end
    P = calc_first_piola_kirchhoff!(storage, mat, params, F, Ḟ, Δt, i)
    PKinv = P * Kinv

    # stabilization
    C₁ = get_zem_stiffness!(storage, params, Kinv, mat.stress_integration, i)

    for bond_id in each_bond_idx(system, i)
        bond = bonds[bond_id]
        j, L = bond.neighbor, bond.length
        ΔXij = get_coordinates_diff(system, i, j)
        Δxij = get_coordinates_diff(storage, i, j)
        l = norm(Δxij)
        ε = (l - L) / L
        stretch_based_failure!(storage, system, bond, params, ε, i, bond_id)

        ωij = influence_function(mat, params, L) * bond_active[bond_id]

        # previous stabilization from Silling:
        # Tij = mat.corr .* params.bc * ωij / ω0 .* (Δxij .- F * ΔXij)

        # improved stabilization from this article:
        # https://doi.org/10.1007/s10409-019-00873-y
        Tij = ωij * C₁ * (Δxij .- F * ΔXij)

        # update of force density
        tij = ωij * PKinv * ΔXij + Tij
        update_add_b_int!(storage, i, tij .* volume[j])
        update_add_b_int!(storage, j, -tij .* volume[i])
    end
    return nothing
end

@inline function influence_function(::CCMaterial, params::CCPointParameters, L::Float64)
    return params.δ / L
end

function calc_deformation_gradient(storage::AbstractStorage, system::BondSystem,
                                   mat::AbstractCorrespondenceMaterial,
                                   params::AbstractPointParameters, i::Int)
    (; bonds, volume) = system
    (; bond_active) = storage
    K = zeros(SMatrix{3,3})
    _F = zeros(SMatrix{3,3})
    _Ḟ = zeros(SMatrix{3,3})
    for bond_id in each_bond_idx(system, i)
        bond = bonds[bond_id]
        j, L = bond.neighbor, bond.length
        ΔXij = get_diff(system.position, i, j)
        Δxij = get_diff(storage.position, i, j)
        Δvij = get_diff(storage.velocity, i, j)
        ωij = influence_function(mat, params, L) * bond_active[bond_id]
        m = ωij * volume[j]
        K += m * ΔXij * ΔXij'
        _F += m * Δxij * ΔXij'
        _Ḟ += m * Δvij * ΔXij'
    end
    Kinv = inv(K)
    F = _F * Kinv
    Ḟ = _Ḟ * Kinv
    return F, Ḟ, Kinv
end

function calc_first_piola_kirchhoff!(storage::CCStorage, mat::CCMaterial,
                                     params::CCPointParameters, F::SMatrix{3,3},
                                     Ḟ::SMatrix{3,3}, Δt::Float64, i::Int)
    init_stress_integration!(storage, mat.stress_integration, F, Ḟ, Δt, i)
    σ = cauchy_stress(mat.constitutive_model, storage, params, F)
    T = stress_integration(storage, mat.stress_integration, σ, Δt, i)
    update_tensor!(storage.stress, i, T)
    P = det(F) * T * inv(F)'
    return P
end

@inline function get_tensor(T::AbstractMatrix, i::Int)
    tensor = SMatrix{3,3}(T[1, i], T[2, i], T[3, i], T[4, i], T[5, i], T[6, i], T[7, i],
                          T[8, i], T[9, i])
    return tensor
end

@inline function update_tensor!(Tₙ::AbstractMatrix, i::Int, Tₙ₊₁::SMatrix{3,3})
    Tₙ[1, i] = Tₙ₊₁[1, 1]
    Tₙ[2, i] = Tₙ₊₁[1, 2]
    Tₙ[3, i] = Tₙ₊₁[1, 3]
    Tₙ[4, i] = Tₙ₊₁[2, 1]
    Tₙ[5, i] = Tₙ₊₁[2, 2]
    Tₙ[6, i] = Tₙ₊₁[2, 3]
    Tₙ[7, i] = Tₙ₊₁[3, 1]
    Tₙ[8, i] = Tₙ₊₁[3, 2]
    Tₙ[9, i] = Tₙ₊₁[3, 3]
    return nothing
end

function containsnan(K::T) where {T<:AbstractArray}
    @simd for i in eachindex(K)
        isnan(K[i]) && return true
    end
    return false
end

function kill_point!(s::AbstractStorage, bd::BondSystem, i::Int)
    s.bond_active[each_bond_idx(bd, i)] .= false
    s.n_active_bonds[i] = 0
    return nothing
end
