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
