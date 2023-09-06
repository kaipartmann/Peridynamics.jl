@timeit TO function apply_bcs!(pdp::PDProblem, t::Float64)
    @sync for bc in pdp.sp.bcs
        Threads.@spawn apply_boundarycondition!(pdp.gs, bc, t)
    end
    return nothing
end

function apply_boundarycondition!(gs::GlobalStorage, bc::VelocityBC, t::Float64)
    value = bc.fun(t)
    dim = bc.dim
    @inbounds @simd for i in bc.point_id_set
        gs.velocity_half[dim, i] = value
    end
    return nothing
end

function apply_boundarycondition!(gs::GlobalStorage, bc::ForceDensityBC, t::Float64)
    value = bc.fun(t)
    dim = bc.dim
    @inbounds @simd for i in bc.point_id_set
        gs.b_ext[dim, i] = value
    end
    return nothing
end

function apply_boundarycondition!(gs::GlobalStorage, bc::PosDepVelBC, t::Float64)
    dim = bc.dim
    @inbounds @simd for i in bc.point_id_set
        gs.velocity_half[dim, i] = bc.fun(gs.position[1, i], gs.position[2, i],
                                          gs.position[3, i], t)
    end
    return nothing
end

@timeit TO function apply_ics!(pdp::PDProblem)
    @sync for ic in pdp.sp.ics
        Threads.@spawn apply_initialcondition!(pdp.gs, ic)
    end
    return nothing
end

function apply_initialcondition!(gs::GlobalStorage, ic::VelocityIC)
    dim = ic.dim
    @inbounds @simd for i in ic.point_id_set
        gs.velocity[dim, i] = ic.val
    end
    return nothing
end
