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

function initialize!(chunk::BodyChunk{BondSystem{EnergySurfaceCorrection}})
    inititalize_mfactor!(chunk)
    return nothing
end

function inititalize_mfactor!(chunk::BodyChunk{BondSystem{EnergySurfaceCorrection}})
    system = chunk.system
    mfactor = system.correction.mfactor
    stendens = zeros(3, n_points)
    for d in 1:3
        defposition .= copy(system.position)
        defposition[d,:] .*= 1.001
        @views calc_stendens!(stendens[d,:], defposition, chunk)
        for i in each_point_idx(chunk)
            param = get_params(chunk, i)
            mfactor[d,i] = 0.5 * param.G * 1e-6 / stendens[d,i]
        end
    end
    return nothing
end

function calc_stendens!(stendens, defposition, chunk)
    system = chunk.system
    for i in each_point_idx(chunk)
        param = get_params(chunk, i)
        temp = 15 * param.G /(2π * param.δ * param.δ * param.δ * param.δ)
        for bond_id in each_bond_idx(system, i)
            bond = system.bonds[bond_id]
            j, L = bond.neighbor, bond.length
            Δxijx = defposition[1, j] - defposition[1, i]
            Δxijy = defposition[2, j] - defposition[2, i]
            Δxijz = defposition[3, j] - defposition[3, i]
            l = sqrt(Δxijx * Δxijx + Δxijy * Δxijy + Δxijz * Δxijz)
            ε = (l - L) / L
            stendens[i] += temp * ε * ε * L * system.volume[j]
        end
    end
    return nothing
end
