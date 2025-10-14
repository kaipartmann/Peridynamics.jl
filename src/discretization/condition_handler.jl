struct ConditionHandler
    loc_point_sets::Dict{Symbol,Vector{Int}}
    single_dim_bcs::Vector{SingleDimBC}
    posdep_single_dim_bcs::Vector{PosDepSingleDimBC}
    pos_single_dim_bcs::Vector{PosSingleDimBC}
    data_bcs::Vector{DataBC}
    free_dofs::Vector{Int}
    constrained_dofs::Vector{Int}
end

function ConditionHandler(body::AbstractBody, system::AbstractSystem)
    (; single_dim_bcs, posdep_single_dim_bcs, pos_single_dim_bcs, data_bcs) = body
    (; point_sets) = body
    loc_point_sets = localized_point_sets(point_sets, system.chunk_handler)
    free_dofs, constrained_dofs = filter_dofs(body, system, loc_point_sets)
    condhandler = ConditionHandler(loc_point_sets, single_dim_bcs, posdep_single_dim_bcs,
                                   pos_single_dim_bcs, data_bcs, free_dofs,
                                   constrained_dofs)
    return condhandler
end

function filter_dofs(body::AbstractBody, system::AbstractSystem,
                     loc_point_sets::Dict{Symbol,Vector{Int}})
    constrained_dofs = Vector{Int}()
    # Collect constrained DOFs from all boundary conditions
    constrained_dofs!(constrained_dofs, body, system, loc_point_sets)
    # Collect free DOFs (those not constrained by any BC)
    free_dofs = Vector{Int}()
    for dof in each_loc_dof(system)
        if !(dof in constrained_dofs)
            push!(free_dofs, dof)
        end
    end
    @assert sort(free_dofs) == sort(setdiff(1:get_n_loc_dof(system), constrained_dofs))
    return free_dofs, constrained_dofs
end

@inline function free_dofs(condhandler::ConditionHandler)
    return condhandler.free_dofs
end
@inline function constrained_dofs(condhandler::ConditionHandler)
    return condhandler.constrained_dofs
end
