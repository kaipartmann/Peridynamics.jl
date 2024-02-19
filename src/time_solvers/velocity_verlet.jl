mutable struct VelocityVerlet <: AbstractTimeSolver
    end_time::Float64
    n_steps::Int
    Δt::Float64
    safety_factor::Float64

    function VelocityVerlet(; time::Real=-1, steps::Int=-1, stepsize::Real=-1,
                            safety_factor::Real=0.7)
        if time > 0 && steps > 0
            msg = "specify either time or number of steps, not both!"
            throw(ArgumentError(msg))
        elseif time < 0 && steps < 0
            msg = "specify either time or number of steps!"
            throw(ArgumentError(msg))
        end
        if !(0 < safety_factor < 1)
            msg = "wrong safety factor specified! condition: 0 < safety_factor < 1"
            throw(ArgumentError(msg))
        end
        if stepsize > 0
            @warn "stepsize specified! Please be sure that the CFD-condition holds!"
        end
        new(time, steps, stepsize, safety_factor)
    end
end

function init_time_solver!(vv::VelocityVerlet, dh::AbstractDataHandler)
    if vv.Δt < 0
        vv.Δt = calc_stable_timestep(dh, vv.safety_factor)
    end
    if vv.end_time < 0
        vv.end_time = vv.n_steps * vv.Δt
    elseif vv.n_steps < 0
        vv.n_steps = vv.end_time ÷ vv.Δt + 1
    end
    velocity_verlet_check(vv)
    return nothing
end

function velocity_verlet_check(vv::VelocityVerlet)
    if vv.end_time < 0
        error("`end_time` of VelocityVerlet smaller than zero!\n")
    end
    if vv.n_steps < 0
        error("`n_steps` of VelocityVerlet smaller than zero!\n")
    end
    if vv.Δt < 0
        error("`Δt` of VelocityVerlet smaller than zero!\n")
    end
    return nothing
end

function calc_timestep(b::AbstractBodyChunk)
    Δt = fill(Inf, length(each_point_idx(b.ch)))
    for point_id in each_point_idx(b.ch)
        pp = get_point_param(b, point_id)
        Δt[point_id] = _calc_timestep(b.discret, pp, point_id)
    end
    return minimum(Δt)
end

function _calc_timestep(bd::BondDiscretization, pp::AbstractPointParameters, point_id::Int)
    dtsum = 0.0
    for bond_id in each_bond_idx(bd, point_id)
        bond = bd.bonds[bond_id]
        dtsum += bd.volume[bond.neighbor] * pp.bc / bond.length
    end
    return sqrt(2 * pp.rho / dtsum)
end
