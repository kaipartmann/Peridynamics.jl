struct NoZEMStabilization <: AbstractZEMStabilization end

function calc_zem_force_density(::NoZEMStabilization, args...)
    return zeros(SVector{3})
end

function get_correction(::AbstractBondSystemMaterial{NoZEMStabilization},
                        n_loc_points::Int, n_points::Int, n_bonds::Int)
    return NoZEMStabilization()
end

struct SillingZEMStabilization <: AbstractZEMStabilization end
