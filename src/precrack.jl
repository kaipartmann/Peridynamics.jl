@timeit TO function define_precrack!(pdp::PDProblem, precrack::PreCrack)
    @inbounds @threads :static for tid in 1:pdp.sp.n_threads
        pdp.tls[tid].n_active_family_members .= 0
        for current_one_ni in pdp.sp.owned_bonds[tid]
            i, j, _, _ = pdp.sp.bond_data[current_one_ni]
            i_is_in_set_a = in(i, precrack.point_id_set_a)
            i_is_in_set_b = in(i, precrack.point_id_set_b)
            j_is_in_set_a = in(j, precrack.point_id_set_a)
            j_is_in_set_b = in(j, precrack.point_id_set_b)
            if i_is_in_set_a && j_is_in_set_b || i_is_in_set_b && j_is_in_set_a
                pdp.gs.bond_failure[current_one_ni] = 0
            end
            li = pdp.tls[tid].pointmap[i]
            pdp.tls[tid].n_active_family_members[li] += pdp.gs.bond_failure[current_one_ni]
            if pdp.sp.unique_bonds
                lj = pdp.tls[tid].pointmap[j]
                pdp.tls[tid].n_active_family_members[lj] += pdp.gs.bond_failure[current_one_ni]
            end
        end
    end
    return nothing
end

@timeit TO function calc_damage!(pdp::PDProblem)
    @inbounds @threads :static for i in 1:pdp.sp.pc.n_points
        pdp.gs.damage[i] = 1 - pdp.gs.n_active_family_members[i] / pdp.sp.n_family_members[i]
    end
    return nothing
end
