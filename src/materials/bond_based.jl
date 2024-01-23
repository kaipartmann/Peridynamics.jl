
"""
    BBMaterial <: AbstractMaterial

Bond based peridynamic material model.

# Fields
- `δ::Float64`: horizon
- `rho::Float64`: density
- `E::Float64`: young's modulus
- `nu::Float64`: poisson ratio
- `G::Float64`: shear modulus
- `K::Float64`: bulk modulus
- `bc::Float64`: bond constant
- `Gc::Float64`: critical energy release rate
- `εc::Float64`: critical bond strain

---
```julia
BBMaterial(; horizon::Real, rho::Real, E::Real, [Gc::Real, epsilon_c::Real])
```

Specify a material with `horizon`, density `rho`, Young's modulus `E` and either the
critical energy release rate `Gc` or the critical bond strain `epsilon_c`.

# Keywords
- `horizon::Real`: horizon
- `rho::Real`:density
- `E::Real`: young's modulus
- `Gc::Real`: critical energy release rate
- `epsilon_c::Real`: critical bond strain
"""
struct BBMaterial <: Material
    δ::Float64
    rho::Float64
    E::Float64
    nu::Float64
    G::Float64
    K::Float64
    bc::Float64
    Gc::Float64
    εc::Float64
end

function BBMaterial(; horizon::Real, rho::Real, E::Real, Gc::Real=-1, epsilon_c::Real=-1)
    nu = 0.25 # limitation of the bond-based formulation
    G = E / (2 * (1 + nu))
    K = E / (3 * (1 - 2 * nu))
    bc = 18 * K / (π * horizon^4)
    # TODO: this logic should be a seperate function
    if (Gc !== -1) && (epsilon_c == -1)
        epsilon_c = sqrt(5.0 * Gc / (9.0 * K * horizon))
    elseif (Gc == -1) && (epsilon_c !== -1)
        Gc = 9.0 / 5.0 * K * horizon * epsilon_c^2
    elseif (Gc !== -1) && (epsilon_c !== -1)
        msg = "Duplicate definition: define either Gc or epsilon_c, not both!"
        throw(ArgumentError(msg))
    elseif (Gc == -1) && (epsilon_c == -1)
        throw(ArgumentError("Either Gc or epsilon_c have to be defined!"))
    end
    return BBMaterial(horizon, rho, E, nu, G, K, bc, Gc, epsilon_c)
end

struct BBStorage <: AbstractStorage
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

function init_storage(::BBMaterial, pbd::PointBondDiscretization,
                      loc_points::UnitRange{Int}, halo_points::Vector{Int})
    n_loc_points = length(loc_points)
    position = copy(pbd.position)
    displacement = zeros(3, n_loc_points)
    velocity = zeros(3, n_loc_points)
    velocity_half = zeros(3, n_loc_points)
    acceleration = zeros(3, n_loc_points)
    b_int = zeros(3, n_loc_points)
    b_ext = zeros(3, n_loc_points)
    damage = zeros(n_loc_points)
    bond_active = ones(Bool, length(pbd.bonds))
    n_active_bonds = copy(pbd.n_neighbors)
    BBStorage(position, displacement, velocity, velocity_half, acceleration, b_int, b_ext,
              damage, bond_active, n_active_bonds)
end

function force_density!(s::BBStorage, d::PointBondDiscretization, mat::BBMaterial,
                        i::Int)
    for bond_id in d.bond_range[i]
        bond = d.bonds[bond_id]
        j, L = bond.j, bond.L

        # current bond length
        Δxijx = s.position[1, j] - s.position[1, i]
        Δxijy = s.position[2, j] - s.position[2, i]
        Δxijz = s.position[3, j] - s.position[3, i]
        l = sqrt(Δxijx * Δxijx + Δxijy * Δxijy + Δxijz * Δxijz)

        # bond strain
        ε = (l - L) / L

        # failure mechanism
        if ε > mat.εc && bond.failure_allowed
            s.bond_active[bond_id] = false
        end
        s.n_active_bonds[i] += s.bond_failure[bond_id]

        # update of force density
        temp = s.bond_active[bond_id] * mat.bc * ε / l * d.volume[j]
        s.b_int[1, i] += temp * Δxijx
        s.b_int[2, i] += temp * Δxijy
        s.b_int[3, i] += temp * Δxijz
    end
    return nothing
end
