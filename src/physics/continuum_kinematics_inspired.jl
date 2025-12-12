"""
    CKIMaterial(; dmgmodel)

A material type used to assign the material of a [`Body`](@ref) with the
continuum-kinematics-inspired peridynamics formulation.

# Keywords
- `dmgmodel::AbstractDamageModel`: Damage model defining the damage behavior. \\
    (default: `CriticalStretch()`)

# Examples

```julia-repl
julia> mat = CKIMaterial()
CKIMaterial(dmgmodel=CriticalStretch())
```

---

```julia
CKIMaterial{DM}
```

Material type for the continuum-kinematics-inspired peridynamics framework.

# Type Parameters
- `DM`: A damage model type. See the constructor docs for more informations.

# Fields
- `dmgmodel::AbstractDamageModel`: Damage model defining the damage behavior. See the
    constructor docs for more informations.

# Allowed material parameters
When using [`material!`](@ref) on a [`Body`](@ref) with `CKIMaterial`, then the following
parameters are allowed:
Material parameters:
- `horizon::Float64`: Radius of point interactions.
- `rho::Float64`: Density.
Elastic parameters:
- `E::Float64`: Young's modulus.
- `nu::Float64`: Poisson's ratio.
- `G::Float64`: Shear modulus.
- `K::Float64`: Bulk modulus.
- `lambda::Float64`: 1st Lamé parameter.
- `mu::Float64`: 2nd Lamé parameter.
Fracture parameters:
- `Gc::Float64`: Critical energy release rate.
- `epsilon_c::Float64`: Critical strain.
Interaction parameters:
- `C1::Float64`: One-neighbor interaction parameter. (default: `0.0`)
- `C2::Float64`: Two-neighbor interaction parameter. (default: `0.0`)
- `C3::Float64`: Two-neighbor interaction parameter. (default: `0.0`)

!!! warning "Specification of interaction parameters"
    If any of the interaction parameters is used with [`material!`](@ref), the Young's
    modulus and Poisson's ratio are ignored and only the specified interaction parameters
    will influence the force density calculated from that interaction.

    If no interaction parameter is specified, then the Young's modulus and Poisson's ratio
    are used to calculate these parameters accordingly to Ekiz, Steinmann, and Javili
    (2022).

!!! note "Elastic parameters"
    Note that exactly two elastic parameters are required to specify a material.
    Please choose two out of the six allowed elastic parameters.

# Allowed export fields
When specifying the `fields` keyword of [`Job`](@ref) for a [`Body`](@ref) with
`CKIMaterial`, the following fields are allowed:
- `position::Matrix{Float64}`: Position of each point.
- `displacement::Matrix{Float64}`: Displacement of each point.
- `velocity::Matrix{Float64}`: Velocity of each point.
- `velocity_half::Matrix{Float64}`: Velocity parameter for Verlet time solver.
- `acceleration::Matrix{Float64}`: Acceleration of each point.
- `b_int::Matrix{Float64}`: Internal force density of each point.
- `b_ext::Matrix{Float64}`: External force density of each point.
- `damage::Vector{Float64}`: Damage of each point.
- `n_active_one_nis::Vector{Int}`: Number of intact one-neighbor interactions of each point.
- `strain_energy_density::Vector{Float64}`: Strain energy density of each point.
"""
struct CKIMaterial{DM} <: AbstractInteractionSystemMaterial
    dmgmodel::DM
    function CKIMaterial(dmgmodel::DM) where DM
        return new{DM}(dmgmodel)
    end
end

function CKIMaterial(; dmgmodel::AbstractDamageModel=CriticalStretch())
    return CKIMaterial(dmgmodel)
end

"""
    CKIPointParameters

$(internal_api_warning())

Type containing the material parameters for a continuum-kinematics-inspired peridynamics
model.

# Fields

- `δ::Float64`: Horizon.
- `rho::Float64`: Density.
- `E::Float64`: Young's modulus.
- `nu::Float64`: Poisson's ratio.
- `G::Float64`: Shear modulus.
- `K::Float64`: Bulk modulus.
- `λ::Float64`: 1st Lamé parameter.
- `μ::Float64`: 2nd Lamé parameter.
- `Gc::Float64`: Critical energy release rate.
- `εc::Float64`: Critical strain.
- `C1::Float64`: Material constant for one-neighbor interactions.
- `C2::Float64`: Material constant for two-neighbor interactions.
- `C3::Float64`: Material constant for three-neighbor interactions.
"""
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

function CKIPointParameters(mat::CKIMaterial, p::Dict{Symbol,Any})
    (; δ, rho, E, nu, G, K, λ, μ, C1, C2, C3) = get_required_point_parameters(mat, p)
    (; Gc, εc) = get_frac_params(mat.dmgmodel, p, δ, K)
    return CKIPointParameters(δ, rho, E, nu, G, K, λ, μ, Gc, εc, C1, C2, C3)
end

@params CKIMaterial CKIPointParameters

@storage CKIMaterial struct CKIStorage <: AbstractStorage
    @lthfield position::Matrix{Float64}
    @pointfield displacement::Matrix{Float64}
    @pointfield velocity::Matrix{Float64}
    @pointfield velocity_half::Matrix{Float64}
    @pointfield velocity_half_old::Matrix{Float64}
    @pointfield acceleration::Matrix{Float64}
    @pointfield b_int::Matrix{Float64}
    @pointfield b_int_old::Matrix{Float64}
    @pointfield b_ext::Matrix{Float64}
    @pointfield density_matrix::Matrix{Float64}
    @pointfield damage::Vector{Float64}
    @pointfield n_active_one_nis::Vector{Int}
    @pointfield strain_energy_density::Vector{Float64}
    one_ni_active::Vector{Bool}
    residual::Vector{Float64}
    displacement_copy::Matrix{Float64}
    b_int_copy::Matrix{Float64}
    temp_force::Vector{Float64}
    Δu::Vector{Float64}
    v_temp::Vector{Float64}
    Jv_temp::Vector{Float64}
end

function init_field(::CKIMaterial, ::AbstractTimeSolver, system::InteractionSystem,
                    ::Val{:strain_energy_density})
    return zeros(get_n_loc_points(system))
end

function force_density_point!(storage::CKIStorage, system::InteractionSystem,
                              mat::CKIMaterial, params::AbstractParameterSetup, t, Δt, i)
    force_density_point_one_ni!(storage, system, mat, params, t, Δt, i)
    has_two_nis(params) && force_density_point_two_ni!(storage, system, mat, params, t, Δt, i)
    has_three_nis(params) && force_density_point_three_ni!(storage, system, mat, params, t, Δt, i)
    return nothing
end

function force_density_point_one_ni!(storage::CKIStorage, system::InteractionSystem,
                                     ::CKIMaterial, params::CKIPointParameters, t, Δt, i)
    for one_ni_id in each_one_ni_idx(system, i)
        one_ni = system.one_nis[one_ni_id]
        j, L = one_ni.neighbor, one_ni.length
        Δxij = get_vector_diff(storage.position, i, j)
        l = norm(Δxij)
        b_int = one_ni_failure(storage, one_ni_id) * params.C1 * (1 / L - 1 / l) *
                system.volume_one_nis[i] .* Δxij
        update_add_vector!(storage.b_int, i, b_int)
    end
    return nothing
end

function force_density_point_one_ni!(storage::CKIStorage, system::InteractionSystem,
                                     ::CKIMaterial, paramhandler::ParameterHandler, t, Δt, i)
    params_i = get_params(paramhandler, i)
    for one_ni_id in each_one_ni_idx(system, i)
        one_ni = system.one_nis[one_ni_id]
        j, L = one_ni.neighbor, one_ni.length
        Δxij = get_vector_diff(storage.position, i, j)
        l = norm(Δxij)
        params_j = get_params(paramhandler, j)
        b_int = one_ni_failure(storage, one_ni_id) * (params_i.C1 + params_j.C1) / 2 *
                (1 / L - 1 / l) * system.volume_one_nis[i] .* Δxij
        update_add_vector!(storage.b_int, i, b_int)
    end
    return nothing
end

function force_density_point_two_ni!(storage::CKIStorage, system::InteractionSystem,
                                     ::CKIMaterial, params::CKIPointParameters, t, Δt, i)
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
                                     ::CKIMaterial, paramhandler::ParameterHandler, t, Δt, i)
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
        # avoid to divide by zero error for failed interactions
        isapprox(surface, 0; atol=eps()) && (surface = 1e-40)
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
                                       ::CKIMaterial, params::CKIPointParameters, t, Δt, i)
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
        # avoid to divide by zero error for failed interactions
        isapprox(volume, 0; atol=eps()) && (volume = 1e-40)
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
                                       t, Δt, i)
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

function strain_energy_density_point!(storage::CKIStorage, system::InteractionSystem,
                                      mat::CKIMaterial, paramsetup::AbstractParameterSetup,
                                      i)
    storage.strain_energy_density[i] = 0.0
    params = get_params(paramsetup, i)
    strain_energy_density_point_one_ni!(storage, system, mat, params, i)
    has_two_nis(params) && strain_energy_density_point_two_ni!(storage, system, mat, params, i)
    has_three_nis(params) && strain_energy_density_point_three_ni!(storage, system, mat, params, i)
    return nothing
end

function strain_energy_density_point_one_ni!(storage::CKIStorage, system::InteractionSystem,
                                             ::CKIMaterial, params::CKIPointParameters, i)
    Ψ = 0.0
    for one_ni_id in each_one_ni_idx(system, i)
        one_ni = system.one_nis[one_ni_id]
        j, L = one_ni.neighbor, one_ni.length
        Δxij = get_vector_diff(storage.position, i, j)
        l = norm(Δxij)
        εl = (l - L) / L
        ψij = 0.5 * params.C1 * εl * εl * L
        failure = one_ni_failure(storage, one_ni_id)
        Ψ += 0.5 * failure * ψij * system.volume_one_nis[i]
    end
    storage.strain_energy_density[i] += Ψ
    return nothing
end

function strain_energy_density_point_two_ni!(storage::CKIStorage, system::InteractionSystem,
                                             ::CKIMaterial, params::CKIPointParameters, i)
    Ψ = 0.0
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
        εa = (surface - surface_ref) / surface_ref
        ψijk = 0.5 * params.C2 * εa * εa * surface_ref
        failure = one_ni_failure(storage, oni_j_id) * one_ni_failure(storage, oni_k_id)
        Ψ += 1/3 * failure * ψijk * system.volume_two_nis[i]
    end
    storage.strain_energy_density[i] += Ψ
    return nothing
end

function strain_energy_density_point_three_ni!(storage::CKIStorage,
                                               system::InteractionSystem, ::CKIMaterial,
                                               params::CKIPointParameters, i)
    Ψ = 0.0
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
        abs_volume = abs(volume)
        abs_volume_ref = abs(volume_ref)
        εv = (abs_volume - abs_volume_ref) / abs_volume_ref
        ψijkl = 0.5 * params.C3 * εv * εv * abs_volume_ref
        failure = one_ni_failure(storage, oni_j_id) * one_ni_failure(storage, oni_k_id) *
                  one_ni_failure(storage, oni_l_id)
        Ψ += 1/4 * failure * ψijkl * system.volume_three_nis[i]
    end
    storage.strain_energy_density[i] += Ψ
    return nothing
end

function export_field(::Val{:strain_energy_density}, mat::CKIMaterial,
                      system::InteractionSystem, storage::CKIStorage,
                      paramsetup::CKIPointParameters, t)
    for i in each_point_idx(system)
        strain_energy_density_point!(storage, system, mat, paramsetup, i)
    end
    return storage.strain_energy_density
end
