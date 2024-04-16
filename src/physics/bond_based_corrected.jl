
"""
    BBCMaterial <: AbstractMaterial

material type for bond-based peridynamic simulations

# Allowed material parameters

- `:horizon::Float64`: radius of point interactions
- `:rho::Float64`: density
- `:E::Float64`: Young's modulus
- `:Gc::Float64`: critical energy release rate
- `:epsilon_c::Float64`: critical strain

# Allowed export fields

- `position::Matrix{Float64}`: position of each point
- `displacement::Matrix{Float64}`: displacement of each point
- `velocity::Matrix{Float64}`: velocity of each point
- `velocity_half::Matrix{Float64}`: velocity parameter for Verlet time solver
- `acceleration::Matrix{Float64}`: acceleration of each point
- `b_int::Matrix{Float64}`: internal force density of each point
- `b_ext::Matrix{Float64}`: external force density of each point
- `damage::Vector{Float64}`: damage of each point
- `n_active_bonds::Vector{Int}`: number of intact bonds for each point
"""
struct BBCMaterial <: AbstractMaterial end

struct BBCPointParameters <: AbstractPointParameters
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

function BBCPointParameters(::BBCMaterial, p::Dict{Symbol,Any})
    δ = get_horizon(p)
    rho = get_density(p)
    if haskey(p, :nu)
        msg = "keyword `nu` is not allowed for BBCMaterial!\n"
        msg *= "Bond-based peridynamics has a limitation on the possion ratio.\n"
        msg *= "Therefore, when using BBCMaterial, `nu` is hardcoded as 1/4.\n"
        throw(ArgumentError(msg))
    else
        p[:nu] = 0.25
    end
    E, nu, G, K, λ, μ = get_elastic_params(p)
    Gc, εc = get_frac_params(p, δ, K)
    bc = 18 * K / (π * δ^4) # bond constant
    return BBCPointParameters(δ, rho, E, nu, G, K, λ, μ, Gc, εc, bc)
end

@params BBCMaterial BBCPointParameters

@system BBCMaterial BondSystem

struct BBCVerletStorage <: AbstractStorage
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
    h::Matrix{Float64}
    sucofa::Vector{Float64}
end

function BBCVerletStorage(::BBCMaterial, ::VelocityVerlet, system::BondSystem, ch)
    n_loc_points = length(ch.loc_points)
    n_points = length(ch.point_ids)
    n_bonds = length(system.bonds)
    position = copy(system.position)
    displacement = zeros(3, n_loc_points)
    velocity = zeros(3, n_loc_points)
    velocity_half = zeros(3, n_loc_points)
    acceleration = zeros(3, n_loc_points)
    b_int = zeros(3, n_loc_points)
    b_ext = zeros(3, n_loc_points)
    damage = zeros(n_loc_points)
    bond_active = ones(Bool, n_bonds)
    n_active_bonds = copy(system.n_neighbors)
    h = zeros(3, n_points)
    sucofa = zeros(n_bonds)
    s = BBCVerletStorage(position, displacement, velocity, velocity_half, acceleration,
                        b_int, b_ext, damage, bond_active, n_active_bonds, h, sucofa)
    return s
end

@storage BBCMaterial VelocityVerlet BBCVerletStorage

@halo_read_fields BBCVerletStorage :position :h

point_data_field(s::BBCVerletStorage, ::Val{:h}) = getfield(s, :h)

function init_chunk!(chunk::AbstractBodyChunk{BBCMaterial})
    chunk.storage.h .= init_surface_correction_factors(chunk)
    return nothing
end

function init_surface_correction_factors(chunk)
    n_points = length(chunk.ch.point_ids)
    # n_bonds = length(chunk.system.bonds)

    h = zeros(3, n_points)
    stendens = zeros(3, n_points)
    defposition = zeros(3, n_points, 3)
    # sucofa = zeros(n_bonds)

    @views for d in 1:3
        defposition[:,:,d] .= copy(chunk.system.position)
        defposition[d,:,d] .*= 1.001
        calc_stendens!(stendens, defposition[:,:,d], chunk, d)
        for i in each_point_idx(chunk)
            param = get_param(chunk, i)
            h[d,i] = 0.5 * param.G * 1e-6 / stendens[d,i]
        end
    end

    # for i in each_point_idx(chunk)
    #     calc_sucofa!(sucofa, chunk.system, h, i)
    # end

    # return sucofa

    return h
end

function calc_stendens!(stendens, defposition, chunk, d)
    system = chunk.system
    for i in each_point_idx(chunk)
        param = get_param(chunk, i)
        temp = 15 * param.G /(2π * param.δ * param.δ * param.δ * param.δ)
        for bond_id in each_bond_idx(system, i)
            bond = system.bonds[bond_id]
            j, L = bond.neighbor, bond.length
            Δxijx = defposition[1, j] - defposition[1, i]
            Δxijy = defposition[2, j] - defposition[2, i]
            Δxijz = defposition[3, j] - defposition[3, i]
            l = sqrt(Δxijx * Δxijx + Δxijy * Δxijy + Δxijz * Δxijz)
            ε = (l - L) / L
            stendens[d, i] += temp * ε * ε * L * system.volume[j]
        end
    end
    return nothing
end

# function calc_sucofa!(sucofa, system, h, i)
#     for bond_id in each_bond_idx(system, i)
#         bond = system.bonds[bond_id]
#         j, L = bond.neighbor, bond.length
#         Δxijx = system.position[1, j] - system.position[1, i]
#         Δxijy = system.position[2, j] - system.position[2, i]
#         Δxijz = system.position[3, j] - system.position[3, i]
#         if abs(Δxijz) <= 1e-10
#             if abs(Δxijy) <= 1e-10
#                 θ = 0.0
#             elseif abs(Δxijx) <= 1e-10
#                 θ = 90 * π / 180
#             else
#                 θ = atan(abs(Δxijy)/abs(Δxijx))
#             end
#             ϕ = 90 * π / 180
#             scx = (h[1,i] + h[1,j]) / 2
#             scy = (h[2,i] + h[2,j]) / 2
#             scz = (h[3,i] + h[3,j]) / 2
#             scr = sqrt(
#                 1/(
#                     (cos(θ) * sin(ϕ))^2 / scx^2 +
#                     (sin(θ) * sin(ϕ))^2 / scy^2 +
#                     cos(ϕ)^2 / scz^2
#                 )
#             )
#         elseif abs(Δxijx) <= 1e-10 && abs(Δxijy) <= 1e-10
#             scz = (h[3,i] + h[3,j]) / 2
#             scr = scz
#         else
#             θ = atan(abs(Δxijy)/abs(Δxijx))
#             ϕ = acos(abs(Δxijz)/L)
#             scx = (h[1,i] + h[1,j]) / 2
#             scy = (h[2,i] + h[2,j]) / 2
#             scz = (h[3,i] + h[3,j]) / 2
#             scr = sqrt(
#                 1/(
#                     (cos(θ) * sin(ϕ))^2 / scx^2 +
#                     (sin(θ) * sin(ϕ))^2 / scy^2 +
#                     cos(ϕ)^2 / scz^2
#                 )
#             )
#         end
#         sucofa[bond_id] = scr
#     end
# end

function calc_scr(Δxijx, Δxijy, Δxijz, L, h, i, j)
    if abs(Δxijz) <= 1e-10
        if abs(Δxijy) <= 1e-10
            θ = 0.0
        elseif abs(Δxijx) <= 1e-10
            θ = 90 * π / 180
        else
            θ = atan(abs(Δxijy)/abs(Δxijx))
        end
        ϕ = 90 * π / 180
        scx = (h[1,i] + h[1,j]) / 2
        scy = (h[2,i] + h[2,j]) / 2
        scz = (h[3,i] + h[3,j]) / 2
        scr = sqrt(
            1/(
                (cos(θ) * sin(ϕ))^2 / scx^2 +
                (sin(θ) * sin(ϕ))^2 / scy^2 +
                cos(ϕ)^2 / scz^2
            )
        )
    elseif abs(Δxijx) <= 1e-10 && abs(Δxijy) <= 1e-10
        scz = (h[3,i] + h[3,j]) / 2
        scr = scz
    else
        θ = atan(abs(Δxijy)/abs(Δxijx))
        ϕ = acos(abs(Δxijz)/L)
        scx = (h[1,i] + h[1,j]) / 2
        scy = (h[2,i] + h[2,j]) / 2
        scz = (h[3,i] + h[3,j]) / 2
        scr = sqrt(
            1/(
                (cos(θ) * sin(ϕ))^2 / scx^2 +
                (sin(θ) * sin(ϕ))^2 / scy^2 +
                cos(ϕ)^2 / scz^2
            )
        )
    end
    return scr
end

function force_density_point!(s::BBCVerletStorage, system::BondSystem, ::BBCMaterial,
                              param::BBCPointParameters, i::Int)
    for bond_id in each_bond_idx(system, i)
        bond = system.bonds[bond_id]
        j, L = bond.neighbor, bond.length

        ΔXijx = system.position[1, j] - system.position[1, i]
        ΔXijy = system.position[2, j] - system.position[2, i]
        ΔXijz = system.position[3, j] - system.position[3, i]
        scr = calc_scr(ΔXijx, ΔXijy, ΔXijz, L, s.h, i, j)

        # current bond length
        Δxijx = s.position[1, j] - s.position[1, i]
        Δxijy = s.position[2, j] - s.position[2, i]
        Δxijz = s.position[3, j] - s.position[3, i]
        l = sqrt(Δxijx * Δxijx + Δxijy * Δxijy + Δxijz * Δxijz)

        # bond strain
        ε = (l - L) / L

        # failure mechanism
        if ε > param.εc && bond.fail_permit
            s.bond_active[bond_id] = false
        end
        s.n_active_bonds[i] += s.bond_active[bond_id]

        # update of force density
        factors = s.bond_active[bond_id] * scr
        temp = factors * param.bc * ε / l * system.volume[j]
        s.b_int[1, i] += temp * Δxijx
        s.b_int[2, i] += temp * Δxijy
        s.b_int[3, i] += temp * Δxijz
    end
    return nothing
end
