
"""
    NoCorrection

A correction handler for materials that use the bond system. If `NoCorrection` is used,
then no correction will be applied.

See also [`BBMaterial`](@ref), [`OSBMaterial`](@ref) for further information on how to use
the correction type.
"""
struct NoCorrection <: AbstractCorrection end

function correction_type(::AbstractBondSystemMaterial{Correction}) where {Correction}
    return Correction
end

function get_correction(::AbstractBondSystemMaterial{NoCorrection}, ::Int, ::Int, ::Int)
    return NoCorrection()
end

@inline function surface_correction_factor(::NoCorrection, ::Int)
    return 1
end

#===========================================================================================
    ENERGY BASED SURFACE CORRECTION
===========================================================================================#

"""
    EnergySurfaceCorrection

A correction handler for materials that use the bond system. If `EnergySurfaceCorrection`
is used, then the energy based surface correction method of Le and Bobaru (2018) is used.

See also [`BBMaterial`](@ref), [`OSBMaterial`](@ref) for further information on how to use
the correction type.
"""
struct EnergySurfaceCorrection <: AbstractCorrection
    mfactor::Matrix{Float64} # multiplication factor mfactor[ndims, npoints]
    scfactor::Vector{Float64} # surface correction factor scfactor[nbonds]
end

function get_correction(::AbstractBondSystemMaterial{EnergySurfaceCorrection},
                        n_loc_points::Int, n_points::Int, n_bonds::Int)
    mfactor = zeros(3, n_points)
    scfactor = ones(n_bonds)
    return EnergySurfaceCorrection(mfactor, scfactor)
end

@inline function surface_correction_factor(correction::EnergySurfaceCorrection,
                                           bond_id::Int)
    return correction.scfactor[bond_id]
end

function initialize!(dh::AbstractThreadsBodyDataHandler{BondSystem{EnergySurfaceCorrection}},
                     solver::AbstractTimeSolver)
    @threads :static for chunk in dh.chunks
        calc_mfactor!(chunk)
    end
    @threads :static for chunk_id in eachindex(dh.chunks)
        exchange_loc_to_halo!(get_mfactor, dh, chunk_id)
        calc_scfactor!(dh.chunks[chunk_id])
    end
    calc_force_density!(dh, 0.0, solver.Δt)
    return nothing
end

function initialize!(dh::AbstractMPIBodyDataHandler{BondSystem{EnergySurfaceCorrection}},
                     solver::AbstractTimeSolver)
    calc_mfactor!(dh.chunk)
    exchange_loc_to_halo!(get_mfactor, dh)
    calc_scfactor!(dh.chunk)
    calc_force_density!(dh, 0.0, solver.Δt)
    return nothing
end

function calc_mfactor!(chunk::AbstractBodyChunk{BondSystem{EnergySurfaceCorrection}})
    (; system, storage, paramsetup) = chunk
    (; mfactor) = system.correction
    ka, a = 1.001, 0.001
    for d in 1:3
        defposition = copy(system.position)
        defposition[d, :] .*= ka
        for i in each_point_idx(chunk)
            stendens, E_mean, nu_mean = stendens_point(system, paramsetup, storage,
                                                       defposition, i)
            mfactor[d, i] = analytical_stendens(E_mean, nu_mean, a) / stendens
        end
    end
    return nothing
end

function stendens_point(system::BondSystem{C}, params::AbstractPointParameters,
                        storage::AbstractStorage, defposition::AbstractMatrix{Float64},
                        i::Int) where {C<:EnergySurfaceCorrection}
    stendens = 0.0
    for bond_id in each_bond_idx(system, i)
        bond = system.bonds[bond_id]
        j, L = bond.neighbor, bond.length
        Δxij = get_vector_diff(defposition, i, j)
        l = norm(Δxij)
        ε = (l - L) / L
        failure = storage.bond_active[bond_id]
        stendens += failure * 0.25 * params.bc * ε * ε * L * system.volume[j]
    end
    return stendens, params.E, params.nu
end

function stendens_point(system::BondSystem{C}, paramhandler::AbstractParameterHandler,
                        storage::AbstractStorage, defposition::AbstractMatrix{Float64},
                        i::Int) where {C<:EnergySurfaceCorrection}
    stendens = 0.0
    params_i = get_params(paramhandler, i)
    E_mean, nu_mean = 0.0, 0.0
    for bond_id in each_bond_idx(system, i)
        bond = system.bonds[bond_id]
        j, L = bond.neighbor, bond.length
        Δxij = get_vector_diff(defposition, i, j)
        l = norm(Δxij)
        ε = (l - L) / L
        params_j = get_params(paramhandler, j)
        bc_mean = (params_i.bc + params_j.bc) / 2
        failure = storage.bond_active[bond_id]
        stendens += failure * 0.25 * bc_mean * ε * ε * L * system.volume[j]
        E_mean += failure * (params_i.E + params_j.E) / 2
        nu_mean += failure * (params_i.nu + params_j.nu) / 2
    end
    E_mean /= storage.n_active_bonds[i]
    nu_mean /= storage.n_active_bonds[i]
    return stendens, E_mean, nu_mean
end

@inline function analytical_stendens(E, nu, a)
    return E / (2 * (1 + nu)) * (nu / (1 - 2 * nu) + 1) * a^2
end

function calc_stendens!(stendens, defposition, chunk)
    system = chunk.system
    for i in each_point_idx(chunk)
        params = get_params(chunk, i)
        for bond_id in each_bond_idx(system, i)
            bond = system.bonds[bond_id]
            j, L = bond.neighbor, bond.length
            Δxijx = defposition[1, j] - defposition[1, i]
            Δxijy = defposition[2, j] - defposition[2, i]
            Δxijz = defposition[3, j] - defposition[3, i]
            l = sqrt(Δxijx * Δxijx + Δxijy * Δxijy + Δxijz * Δxijz)
            ε = (l - L) / L
            stendens[i] += 0.25 * params.bc * ε * ε * L * system.volume[j]
        end
    end
    return nothing
end

@inline get_mfactor(chunk::AbstractBodyChunk) = chunk.system.correction.mfactor

function calc_scfactor!(chunk::AbstractBodyChunk)
    system = chunk.system
    mfactor = system.correction.mfactor
    scfactor = system.correction.scfactor
    for i in each_point_idx(chunk)
        for bond_id in each_bond_idx(system, i)
            bond = system.bonds[bond_id]
            j, L = bond.neighbor, bond.length
            Δxijx = system.position[1, j] - system.position[1, i]
            Δxijy = system.position[2, j] - system.position[2, i]
            Δxijz = system.position[3, j] - system.position[3, i]
            if abs(Δxijz) <= 1e-10
                if abs(Δxijy) <= 1e-10
                    θ = 0.0
                elseif abs(Δxijx) <= 1e-10
                    θ = 90 * π / 180
                else
                    θ = atan(abs(Δxijy) / abs(Δxijx))
                end
                ϕ = 90 * π / 180
                scx = (mfactor[1, i] + mfactor[1, j]) / 2
                scy = (mfactor[2, i] + mfactor[2, j]) / 2
                scz = (mfactor[3, i] + mfactor[3, j]) / 2
                scr = sqrt(1 / ((cos(θ) * sin(ϕ))^2 / scx^2 +
                            (sin(θ) * sin(ϕ))^2 / scy^2 +
                            cos(ϕ)^2 / scz^2))
            elseif abs(Δxijx) <= 1e-10 && abs(Δxijy) <= 1e-10
                scz = (mfactor[3, i] + mfactor[3, j]) / 2
                scr = scz
            else
                θ = atan(abs(Δxijy) / abs(Δxijx))
                ϕ = acos(abs(Δxijz) / L)
                scx = (mfactor[1, i] + mfactor[1, j]) / 2
                scy = (mfactor[2, i] + mfactor[2, j]) / 2
                scz = (mfactor[3, i] + mfactor[3, j]) / 2
                scr = sqrt(1 / ((cos(θ) * sin(ϕ))^2 / scx^2 +
                            (sin(θ) * sin(ϕ))^2 / scy^2 +
                            cos(ϕ)^2 / scz^2))
            end
            scfactor[bond_id] = scr
        end
    end
    return nothing
end
