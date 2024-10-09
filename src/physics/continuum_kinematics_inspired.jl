"""
    CKIMaterial()

A material type used to assign the material of a [`Body`](@ref) with the
continuum-kinematics-inspired peridynamics fomulation.

# Examples

```julia-repl
julia> mat = CKIMaterial()
CKIMaterial()
```

---

```julia
CKIMaterial
```

Material type for the continuum-kinematics-inspired peridynamics framework.

# Allowed material parameters
When using [`material!`](@ref) on a [`Body`](@ref) with `CKIMaterial`, then the following
parameters are allowed:
- `horizon::Float64`: Radius of point interactions
- `rho::Float64`: Density
- `E::Float64`: Young's modulus
- `nu::Float64`: Poisson's ratio
- `Gc::Float64`: Critical energy release rate
- `epsilon_c::Float64`: Critical strain
- `C1::Float64`: One-neighbor interaction parameter (default: `0.0`)
- `C2::Float64`: Two-neighbor interaction parameter (default: `0.0`)
- `C3::Float64`: Two-neighbor interaction parameter (default: `0.0`)

!!! warning "Specification of interaction parameters"
    If any of the interaction parameters is used with [`material!`](@ref), the the Young's
    modulus and Poisson's ratio are ignored and only the specified interaction parameters
    will influence the force density calculated from that interaction.

    If no interaction parameter is specified, then the Young's modulus and Poisson's ratio
    are used to calculate these parameters accordingly to Ekiz, Steinmann, and Javili
    (2022).

# Allowed export fields
When specifying the `fields` keyword of [`Job`](@ref) for a [`Body`](@ref) with
`CKIMaterial`, the following fields are allowed:
- `position::Matrix{Float64}`: Position of each point
- `displacement::Matrix{Float64}`: Displacement of each point
- `velocity::Matrix{Float64}`: Velocity of each point
- `velocity_half::Matrix{Float64}`: Velocity parameter for Verlet time solver
- `acceleration::Matrix{Float64}`: Acceleration of each point
- `b_int::Matrix{Float64}`: Internal force density of each point
- `b_ext::Matrix{Float64}`: External force density of each point
- `damage::Vector{Float64}`: Damage of each point
- `n_active_one_nis::Vector{Int}`: Number of intact one-neighbor interactions of each point
"""
struct CKIMaterial <: AbstractInteractionSystemMaterial end

@params CKIMaterial struct CKIPointParameters <: AbstractPointParameters
    @required_params
    C1::Float64
    C2::Float64
    C3::Float64
end

function calc_parameter(::CKIMaterial, ::Val{:C1}, params::NamedTuple)
    return 30 / π * params.μ / params.δ^4
end

function calc_parameter(::CKIMaterial, ::Val{:C2}, params::NamedTuple)
    return 0.0
end

function calc_parameter(::CKIMaterial, ::Val{:C3}, params::NamedTuple)
    return 32 / π^4 * (params.λ - params.μ) / params.δ^12
end

# struct CKIPointParameters <: AbstractPointParameters
#     δ::Float64
#     rho::Float64
#     E::Float64
#     nu::Float64
#     G::Float64
#     K::Float64
#     λ::Float64
#     μ::Float64
#     Gc::Float64
#     εc::Float64
#     C1::Float64
#     C2::Float64
#     C3::Float64
# end

# function CKIPointParameters(::CKIMaterial, p::Dict{Symbol,Any})
#     δ = get_horizon(p)
#     rho = get_density(p)
#     E, nu, G, K, λ, μ = get_elastic_params(p)
#     Gc, εc = get_frac_params(p, δ, K)
#     C1, C2, C3 = cki_parameters(p, δ, λ, μ)
#     return CKIPointParameters(δ, rho, E, nu, G, K, λ, μ, Gc, εc, C1, C2, C3)
# end

# function cki_parameters(p::Dict{Symbol,Any}, δ, λ, μ)
#     if haskey(p, :C1)
#         C1::Float64 = float(p[:C1])
#     else
#         C1 = 0.0
#     end

#     if haskey(p, :C2)
#         C2::Float64 = float(p[:C2])
#     else
#         C2 = 0.0
#     end

#     if haskey(p, :C3)
#         C3::Float64 = float(p[:C3])
#     else
#         C3 = 0.0
#     end

#     if C1 ≈ 0 && C2 ≈ 0 && C3 ≈ 0
#         C1 = 30 / π * μ / δ^4
#         C2 = 0.0
#         C3 = 32 / π^4 * (λ - μ) / δ^12
#     else
#         msg = "parameters for CKIMaterial specified manually!\n"
#         msg *= "Be careful when adjusting these parameters to avoid unexpected outcomes!"
#         @mpiroot @warn msg
#     end

#     return C1, C2, C3
# end

# function allowed_material_kwargs(::CKIMaterial)
#     return (DEFAULT_POINT_KWARGS..., :C1, :C2, :C3)
# end

# @params CKIMaterial CKIPointParameters

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
    one_ni_active::Vector{Bool}
    n_active_one_nis::Vector{Int}
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
    one_ni_active = ones(Bool, length(system.one_nis))
    n_active_one_nis = copy(system.n_one_nis)
    s = CKIVerletStorage(position, displacement, velocity, velocity_half, acceleration,
                         b_int, b_ext, damage, one_ni_active, n_active_one_nis)
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
    one_ni_active::Vector{Bool}
    n_active_one_nis::Vector{Int}
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
    one_ni_active = ones(Bool, length(system.one_nis))
    n_active_one_nis = copy(system.n_one_nis)
    s = CKIRelaxationStorage(position, displacement, velocity, velocity_half,
                             velocity_half_old, b_int, b_int_old, b_ext, density_matrix,
                             damage, one_ni_active, n_active_one_nis)
    return s
end

@storage CKIMaterial DynamicRelaxation CKIRelaxationStorage

@loc_to_halo_fields CKIRelaxationStorage :position

const CKIStorage = Union{CKIVerletStorage,CKIRelaxationStorage}

function force_density_point!(storage::CKIStorage, system::InteractionSystem,
                              mat::CKIMaterial, params::AbstractParameterSetup, i::Int)
    force_density_point_one_ni!(storage, system, mat, params, i)
    has_two_nis(params) && force_density_point_two_ni!(storage, system, mat, params, i)
    has_three_nis(params) && force_density_point_three_ni!(storage, system, mat, params, i)
    return nothing
end

function force_density_point_one_ni!(storage::CKIStorage, system::InteractionSystem,
                                     ::CKIMaterial, params::CKIPointParameters, i::Int)
    for one_ni_id in each_one_ni_idx(system, i)
        one_ni = system.one_nis[one_ni_id]
        j, L = one_ni.neighbor, one_ni.length
        Δxij = get_coordinates_diff(storage, i, j)
        l = norm(Δxij)
        ε = (l - L) / L
        stretch_based_failure!(storage, system, one_ni, params, ε, i, one_ni_id)
        b_int = one_ni_failure(storage, one_ni_id) * params.C1 * (1 / L - 1 / l) *
                system.volume_one_nis[i] .* Δxij
        update_add_b_int!(storage, i, b_int)
    end
    return nothing
end

function force_density_point_one_ni!(storage::CKIStorage, system::InteractionSystem,
                                     ::CKIMaterial, paramhandler::ParameterHandler, i::Int)
    params_i = get_params(paramhandler, i)
    for one_ni_id in each_one_ni_idx(system, i)
        one_ni = system.one_nis[one_ni_id]
        j, L = one_ni.neighbor, one_ni.length
        Δxij = get_coordinates_diff(storage, i, j)
        l = norm(Δxij)
        ε = (l - L) / L
        stretch_based_failure!(storage, system, one_ni, params_i, ε, i, one_ni_id)
        params_j = get_params(paramhandler, j)
        b_int = one_ni_failure(storage, one_ni_id) * (params_i.C1 + params_j.C1) / 2 *
                (1 / L - 1 / l) * system.volume_one_nis[i] .* Δxij
        update_add_b_int!(storage, i, b_int)
    end
    return nothing
end

function force_density_point_two_ni!(storage::CKIStorage, system::InteractionSystem,
                                     ::CKIMaterial, params::CKIPointParameters, i::Int)
    for two_ni_id in each_two_ni_idx(system, i)
        two_ni = system.two_nis[two_ni_id]
        oni_j_id, oni_k_id, surface_ref = two_ni.oni_j, two_ni.oni_k, two_ni.surface
        j = system.one_nis[oni_j_id].neighbor
        k = system.one_nis[oni_k_id].neighbor
        Δxijx = storage.position[1, j] - storage.position[1, i]
        Δxijy = storage.position[2, j] - storage.position[2, i]
        Δxijz = storage.position[3, j] - storage.position[3, i]
        Δxikx = storage.position[1, k] - storage.position[1, i]
        Δxiky = storage.position[2, k] - storage.position[2, i]
        Δxikz = storage.position[3, k] - storage.position[3, i]
        aijkx = Δxijy * Δxikz - Δxijz * Δxiky
        aijky = Δxijz * Δxikx - Δxijx * Δxikz
        aijkz = Δxijx * Δxiky - Δxijy * Δxikx
        surface = sqrt(aijkx * aijkx + aijky * aijky + aijkz * aijkz)
        if surface == 0 # to avoid divide by zero error for failed interactions
            surface = 1e-40
        end
        failure = storage.one_ni_active[oni_j_id] * storage.one_ni_active[oni_k_id]
        _temp = failure * 2 * params.C2 * (1 / surface_ref - 1 / surface)
        temp = _temp * system.volume_two_nis[i]
        storage.b_int[1, i] += temp * (Δxiky * aijkz - Δxikz * aijky)
        storage.b_int[2, i] += temp * (Δxikz * aijkx - Δxikx * aijkz)
        storage.b_int[3, i] += temp * (Δxikx * aijky - Δxiky * aijkx)
        # variable switch: i,j = j,i
        storage.b_int[1, i] += temp * (Δxijy * (Δxikx * Δxijy - Δxiky * Δxijx) -
                                Δxijz * (Δxikz * Δxijx - Δxikx * Δxijz))
        storage.b_int[2, i] += temp * (Δxijz * (Δxiky * Δxijz - Δxikz * Δxijy) -
                                Δxijx * (Δxikx * Δxijy - Δxiky * Δxijx))
        storage.b_int[3, i] += temp * (Δxijx * (Δxikz * Δxijx - Δxikx * Δxijz) -
                                Δxijy * (Δxiky * Δxijz - Δxikz * Δxijy))
    end
    return nothing
end

function force_density_point_two_ni!(storage::CKIStorage, system::InteractionSystem,
                                     ::CKIMaterial, paramhandler::ParameterHandler, i::Int)
    params_i = get_params(paramhandler, i)
    for two_ni_id in each_two_ni_idx(system, i)
        two_ni = system.two_nis[two_ni_id]
        oni_j_id, oni_k_id, surface_ref = two_ni.oni_j, two_ni.oni_k, two_ni.surface
        j = system.one_nis[oni_j_id].neighbor
        k = system.one_nis[oni_k_id].neighbor
        Δxijx = storage.position[1, j] - storage.position[1, i]
        Δxijy = storage.position[2, j] - storage.position[2, i]
        Δxijz = storage.position[3, j] - storage.position[3, i]
        Δxikx = storage.position[1, k] - storage.position[1, i]
        Δxiky = storage.position[2, k] - storage.position[2, i]
        Δxikz = storage.position[3, k] - storage.position[3, i]
        aijkx = Δxijy * Δxikz - Δxijz * Δxiky
        aijky = Δxijz * Δxikx - Δxijx * Δxikz
        aijkz = Δxijx * Δxiky - Δxijy * Δxikx
        surface = sqrt(aijkx * aijkx + aijky * aijky + aijkz * aijkz)
        if surface == 0 # to avoid divide by zero error for failed interactions
            surface = 1e-40
        end
        failure = storage.one_ni_active[oni_j_id] * storage.one_ni_active[oni_k_id]
        params_j = get_params(paramhandler, j)
        params_k = get_params(paramhandler, k)
        C2_effective = (params_i.C2 + params_j.C2 + params_k.C2) / 3
        _temp = failure * 2 * C2_effective * (1 / surface_ref - 1 / surface)
        temp = _temp * system.volume_two_nis[i]
        storage.b_int[1, i] += temp * (Δxiky * aijkz - Δxikz * aijky)
        storage.b_int[2, i] += temp * (Δxikz * aijkx - Δxikx * aijkz)
        storage.b_int[3, i] += temp * (Δxikx * aijky - Δxiky * aijkx)
        # variable switch: i,j = j,i
        storage.b_int[1, i] += temp * (Δxijy * (Δxikx * Δxijy - Δxiky * Δxijx) -
                                Δxijz * (Δxikz * Δxijx - Δxikx * Δxijz))
        storage.b_int[2, i] += temp * (Δxijz * (Δxiky * Δxijz - Δxikz * Δxijy) -
                                Δxijx * (Δxikx * Δxijy - Δxiky * Δxijx))
        storage.b_int[3, i] += temp * (Δxijx * (Δxikz * Δxijx - Δxikx * Δxijz) -
                                Δxijy * (Δxiky * Δxijz - Δxikz * Δxijy))
    end
    return nothing
end

function force_density_point_three_ni!(storage::CKIStorage, system::InteractionSystem,
                                       ::CKIMaterial, params::CKIPointParameters, i::Int)
    for three_ni_id in each_three_ni_idx(system, i)
        three_ni = system.three_nis[three_ni_id]
        oni_j_id = three_ni.oni_j
        oni_k_id = three_ni.oni_k
        oni_l_id = three_ni.oni_l
        volume_ref = three_ni.volume
        j = system.one_nis[oni_j_id].neighbor
        k = system.one_nis[oni_k_id].neighbor
        l = system.one_nis[oni_l_id].neighbor
        Δxijx = storage.position[1, j] - storage.position[1, i]
        Δxijy = storage.position[2, j] - storage.position[2, i]
        Δxijz = storage.position[3, j] - storage.position[3, i]
        Δxikx = storage.position[1, k] - storage.position[1, i]
        Δxiky = storage.position[2, k] - storage.position[2, i]
        Δxikz = storage.position[3, k] - storage.position[3, i]
        Δxilx = storage.position[1, l] - storage.position[1, i]
        Δxily = storage.position[2, l] - storage.position[2, i]
        Δxilz = storage.position[3, l] - storage.position[3, i]
        # ijk
        aijkx = Δxijy * Δxikz - Δxijz * Δxiky
        aijky = Δxijz * Δxikx - Δxijx * Δxikz
        aijkz = Δxijx * Δxiky - Δxijy * Δxikx
        volume = aijkx * Δxilx + aijky * Δxily + aijkz * Δxilz
        if volume == 0 # avoid to divide by zero error for failed interactions
            volume = 1e-40
        end
        abs_volume = abs(volume)
        failure = storage.one_ni_active[oni_j_id] * storage.one_ni_active[oni_k_id] *
                  storage.one_ni_active[oni_l_id]
        _temp = failure * 3 * params.C3 * (1 / volume_ref - 1 / abs_volume) * volume
        temp = _temp * system.volume_three_nis[i]
        storage.b_int[1, i] += temp * (Δxiky * Δxilz - Δxikz * Δxily)
        storage.b_int[2, i] += temp * (Δxikz * Δxilx - Δxikx * Δxilz)
        storage.b_int[3, i] += temp * (Δxikx * Δxily - Δxiky * Δxilx)
        # kij  |  i->k, j->i, k->j  |  vkij == v
        storage.b_int[1, i] += temp * (Δxijy * Δxikz - Δxijz * Δxiky)
        storage.b_int[2, i] += temp * (Δxijz * Δxikx - Δxijx * Δxikz)
        storage.b_int[3, i] += temp * (Δxijx * Δxiky - Δxijy * Δxikx)
        # jki  |  i->j, j->k, k->i  |  vjki == v
        storage.b_int[1, i] += temp * (Δxily * Δxijz - Δxilz * Δxijy)
        storage.b_int[2, i] += temp * (Δxilz * Δxijx - Δxilx * Δxijz)
        storage.b_int[3, i] += temp * (Δxilx * Δxijy - Δxily * Δxijx)
        # ikj  |  i->i, j->k, k->j  |  vikj == -v
        storage.b_int[1, i] += temp * (Δxily * Δxikz - Δxilz * Δxiky)
        storage.b_int[2, i] += temp * (Δxilz * Δxikx - Δxilx * Δxikz)
        storage.b_int[3, i] += temp * (Δxilx * Δxiky - Δxily * Δxikx)
        # kji  |  i->k, j->j, k->i  |  vkji == -v
        storage.b_int[1, i] += temp * (Δxiky * Δxijz - Δxikz * Δxijy)
        storage.b_int[2, i] += temp * (Δxikz * Δxijx - Δxikx * Δxijz)
        storage.b_int[3, i] += temp * (Δxikx * Δxijy - Δxiky * Δxijx)
        # jik  |  i->j, j->i, k->k  |  vjik == -v
        storage.b_int[1, i] += temp * (Δxijy * Δxilz - Δxijz * Δxily)
        storage.b_int[2, i] += temp * (Δxijz * Δxilx - Δxijx * Δxilz)
        storage.b_int[3, i] += temp * (Δxijx * Δxily - Δxijy * Δxilx)
    end
    return nothing
end

function force_density_point_three_ni!(storage::CKIStorage, system::InteractionSystem,
                                       ::CKIMaterial, paramhandler::ParameterHandler,
                                       i::Int)
    params_i = get_params(paramhandler, i)
    for three_ni_id in each_three_ni_idx(system, i)
        three_ni = system.three_nis[three_ni_id]
        oni_j_id = three_ni.oni_j
        oni_k_id = three_ni.oni_k
        oni_l_id = three_ni.oni_l
        volume_ref = three_ni.volume
        j = system.one_nis[oni_j_id].neighbor
        k = system.one_nis[oni_k_id].neighbor
        l = system.one_nis[oni_l_id].neighbor
        Δxijx = storage.position[1, j] - storage.position[1, i]
        Δxijy = storage.position[2, j] - storage.position[2, i]
        Δxijz = storage.position[3, j] - storage.position[3, i]
        Δxikx = storage.position[1, k] - storage.position[1, i]
        Δxiky = storage.position[2, k] - storage.position[2, i]
        Δxikz = storage.position[3, k] - storage.position[3, i]
        Δxilx = storage.position[1, l] - storage.position[1, i]
        Δxily = storage.position[2, l] - storage.position[2, i]
        Δxilz = storage.position[3, l] - storage.position[3, i]
        # ijk
        aijkx = Δxijy * Δxikz - Δxijz * Δxiky
        aijky = Δxijz * Δxikx - Δxijx * Δxikz
        aijkz = Δxijx * Δxiky - Δxijy * Δxikx
        volume = aijkx * Δxilx + aijky * Δxily + aijkz * Δxilz
        if volume == 0 # avoid to divide by zero error for failed interactions
            volume = 1e-40
        end
        abs_volume = abs(volume)
        failure = storage.one_ni_active[oni_j_id] * storage.one_ni_active[oni_k_id] *
                  storage.one_ni_active[oni_l_id]
        params_j = get_params(paramhandler, j)
        params_k = get_params(paramhandler, k)
        params_l = get_params(paramhandler, l)
        C3_effective = (params_i.C3 + params_j.C3 + params_k.C3 + params_l.C3) / 4
        _temp = failure * 3 * C3_effective * (1 / volume_ref - 1 / abs_volume) * volume
        temp = _temp * system.volume_three_nis[i]
        storage.b_int[1, i] += temp * (Δxiky * Δxilz - Δxikz * Δxily)
        storage.b_int[2, i] += temp * (Δxikz * Δxilx - Δxikx * Δxilz)
        storage.b_int[3, i] += temp * (Δxikx * Δxily - Δxiky * Δxilx)
        # kij  |  i->k, j->i, k->j  |  vkij == v
        storage.b_int[1, i] += temp * (Δxijy * Δxikz - Δxijz * Δxiky)
        storage.b_int[2, i] += temp * (Δxijz * Δxikx - Δxijx * Δxikz)
        storage.b_int[3, i] += temp * (Δxijx * Δxiky - Δxijy * Δxikx)
        # jki  |  i->j, j->k, k->i  |  vjki == v
        storage.b_int[1, i] += temp * (Δxily * Δxijz - Δxilz * Δxijy)
        storage.b_int[2, i] += temp * (Δxilz * Δxijx - Δxilx * Δxijz)
        storage.b_int[3, i] += temp * (Δxilx * Δxijy - Δxily * Δxijx)
        # ikj  |  i->i, j->k, k->j  |  vikj == -v
        storage.b_int[1, i] += temp * (Δxily * Δxikz - Δxilz * Δxiky)
        storage.b_int[2, i] += temp * (Δxilz * Δxikx - Δxilx * Δxikz)
        storage.b_int[3, i] += temp * (Δxilx * Δxiky - Δxily * Δxikx)
        # kji  |  i->k, j->j, k->i  |  vkji == -v
        storage.b_int[1, i] += temp * (Δxiky * Δxijz - Δxikz * Δxijy)
        storage.b_int[2, i] += temp * (Δxikz * Δxijx - Δxikx * Δxijz)
        storage.b_int[3, i] += temp * (Δxikx * Δxijy - Δxiky * Δxijx)
        # jik  |  i->j, j->i, k->k  |  vjik == -v
        storage.b_int[1, i] += temp * (Δxijy * Δxilz - Δxijz * Δxily)
        storage.b_int[2, i] += temp * (Δxijz * Δxilx - Δxijx * Δxilz)
        storage.b_int[3, i] += temp * (Δxijx * Δxily - Δxijy * Δxilx)
    end
    return nothing
end
