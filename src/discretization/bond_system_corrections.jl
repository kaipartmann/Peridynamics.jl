
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
    (; mat, system, storage, paramsetup) = chunk
    (; mfactor) = system.correction
    λ = 1.001 # deformation stretch factor
    for d in 1:3
        # reset to undeformed positions
        storage.position .= system.position
        # apply uniform deformation only in dimension d
        @views storage.position[d, :] .= λ .* system.position[d, :]
        # calculate the mfactor for dimension d
        for i in each_point_idx(chunk)
            lame = get_averaged_lame_parameters(system, storage, paramsetup, i)
            Ψ_theory = stendens_uniext_small_strain(lame, λ)
            strain_energy_density_point!(storage, system, mat, paramsetup, i)
            Ψ_pd = storage.strain_energy_density[i]
            mfactor[d, i] = Ψ_theory / Ψ_pd
        end
        # this is needed since the position is modified during the calculation
        storage.position .= system.position
        # also the artificial strain energy density and b_int should be reset to zero
        storage.strain_energy_density .= 0.0
        storage.b_int .= 0.0
    end
    return nothing
end

function get_averaged_lame_parameters(::BondSystem, ::AbstractStorage,
                                      params::AbstractPointParameters, i)
    return SVector{2,Float64}(params.λ, params.μ)
end

function get_averaged_lame_parameters(system::BondSystem, storage::AbstractStorage,
                                      paramsetup::AbstractParameterHandler, i)
    (; bonds) = system
    (; bond_active) = storage
    lame = zero(SVector{2,Float64}) # lame[1] = λ, lame[2] = μ
    params_i = get_params(paramsetup, i)
    n_active_bonds = 0
    for bond_id in each_bond_idx(system, i)
        if bond_active[bond_id]
            bond = bonds[bond_id]
            j = bond.neighbor
            params_j = get_params(paramsetup, j)
            λ = (params_i.λ + params_j.λ) / 2
            μ = (params_i.μ + params_j.μ) / 2
            lame += SVector{2,Float64}(λ, μ)
            n_active_bonds += 1
        end
    end
    lame /= n_active_bonds
    return lame
end

@inline function stendens_uniext_small_strain(lame, λ)
    return (0.5 * lame[1] + lame[2]) * (λ - 1)^2
end

@inline function stendens_uniext_finite_strain(lame, λ)
    return (lame[1] + 2 * lame[2]) / 8 * (λ^2 - 1)^2
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
