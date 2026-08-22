# The storage of every material has the fields every time solver needs, allocated with the right
# sizes: point fields for the local points, halo fields (`@lth`/`@htl`) for all points
# of the chunk including its halo, bond fields for the bonds, and the fields of the other solvers
# empty. The expectations are a table: the fields all bond-system storages share, plus what each
# material adds. Adding a material means adding one line to `MATERIAL_FIELDS`.

@testmodule StorageFields begin
    using Peridynamics, Test

    # A shape spec is a tuple of (d, count): d rows and `count` columns, where `count` is
    # `:all` (points incl. halo), `:loc` (local points), `:bonds`, `:dof` (local dof),
    # `:neighbors` (max. neighbors of a point); `(d,)` with d === 0 is empty, `(:state,)` is
    # the nested state of a constitutive model without state. Vectors use d = 1.
    const COMMON_FIELDS = Dict(
        :position => (3, :all),
        :displacement => (3, :loc),
        :damage => (1, :loc),
    )

    # fields the bond systems share, beyond `COMMON_FIELDS`
    const BOND_SYSTEM_FIELDS = Dict(
        :n_active_bonds => (1, :loc),
        :bond_active => (1, :bonds),
    )

    # per solver: the fields it needs and how, and the fields of the others it leaves empty
    const SOLVER_FIELDS = Dict(
        VelocityVerlet => Dict(
            :velocity => (3, :loc), :velocity_half => (3, :loc), :acceleration => (3, :loc),
            :b_ext => (3, :loc), :velocity_half_old => (0,), :b_int_old => (0,),
            :density_matrix => (0,), :residual => (0,), :displacement_copy => (0,),
            :b_int_copy => (0,), :temp_force => (0,), :Δu => (0,), :v_temp => (0,),
            :Jv_temp => (0,),
        ),
        DynamicRelaxation => Dict(
            :velocity => (3, :loc), :velocity_half => (3, :loc), :acceleration => (3, :loc),
            :b_ext => (3, :loc), :velocity_half_old => (3, :loc), :b_int_old => (3, :loc),
            :density_matrix => (3, :loc), :residual => (0,), :displacement_copy => (0,),
            :b_int_copy => (0,), :temp_force => (0,), :Δu => (0,), :v_temp => (0,),
            :Jv_temp => (0,),
        ),
        NewtonKrylov => Dict(
            :velocity => (0,), :velocity_half => (0,), :acceleration => (0,),
            :b_ext => (3, :all), :velocity_half_old => (0,), :b_int_old => (0,),
            :density_matrix => (0,), :residual => (1, :dof), :displacement_copy => (3, :loc),
            :b_int_copy => (3, :all), :temp_force => (1, :dof), :Δu => (1, :dof),
            :v_temp => (1, :dof), :Jv_temp => (1, :dof),
        ),
    )

    # The internal force density is a point field unless the material declares it `@htl`,
    # and the Newton-Krylov solver always allocates it for all points.
    b_int_spec(mat, solver) = (3, :all)
    b_int_spec(::Union{BBMaterial,GBBMaterial,CKIMaterial}, ::Union{VelocityVerlet,DynamicRelaxation}) = (3, :loc)

    # `velocity_half` is `@lth` for the rotated formulations, which need the velocities of
    # the halo for the velocity gradient.
    velocity_half_spec(mat, solver) = SOLVER_FIELDS[typeof(solver)][:velocity_half]
    velocity_half_spec(::Union{CRMaterial,RKCRMaterial}, ::Union{VelocityVerlet,DynamicRelaxation}) = (3, :all)

    # what each material adds to the common fields
    const MATERIAL_FIELDS = [
        BBMaterial() => Dict(:strain_energy_density => (1, :loc), :bond_length => (1, :bonds)),
        DHBBMaterial() => Dict(:strain_energy_density => (1, :loc), :bond_length => (1, :bonds)),
        GBBMaterial() => Dict(:strain_energy_density => (1, :loc), :bond_length => (1, :bonds),
                              :weighted_volume => (1, :loc)),
        OSBMaterial() => Dict(:strain_energy_density => (1, :loc), :bond_length => (1, :bonds)),
        CMaterial() => Dict(:strain_energy_density => (1, :loc), :defgrad => (9, :loc),
                            :cauchy_stress => (9, :loc), :von_mises_stress => (1, :loc),
                            :cm_state => (:state,)),
        CRMaterial() => Dict(:strain_energy_density => (1, :loc), :defgrad => (9, :loc),
                             :cauchy_stress => (9, :loc), :von_mises_stress => (1, :loc),
                             :unrotated_stress => (9, :loc), :left_stretch => (9, :loc),
                             :rotation => (9, :loc), :zem_stiffness_rotated => (3, 3, 3, 3)),
        RKCMaterial() => Dict(:strain_energy_density => (1, :loc), :defgrad => (9, :all),
                              :weighted_volume => (1, :all), :update_gradients => (1, :loc),
                              :cauchy_stress => (9, :loc), :von_mises_stress => (1, :loc),
                              :gradient_weight => (3, :bonds),
                              :bond_first_piola_kirchhoff => (9, :bonds),
                              :cm_state => (:state,)),
        RKCRMaterial() => Dict(:strain_energy_density => (1, :loc), :defgrad => (9, :all),
                               :defgrad_dot => (9, :all), :weighted_volume => (1, :all),
                               :update_gradients => (1, :loc), :cauchy_stress => (9, :loc),
                               :von_mises_stress => (1, :loc), :gradient_weight => (3, :bonds),
                               :bond_first_piola_kirchhoff => (9, :bonds),
                               :left_stretch => (9, :bonds), :rotation => (9, :bonds),
                               :bond_unrot_cauchy_stress => (9, :bonds)),
        BACMaterial() => Dict(:stress => (9, :loc), :von_mises_stress => (1, :loc),
                              :bond_stress => (9, :neighbors), :cm_state => (:state,)),
        CKIMaterial() => Dict(:strain_energy_density => (1, :loc), :n_active_one_nis => (1, :loc),
                              :one_ni_active => (1, :bonds)),
    ]

    # the expected size of a field with shape `spec`, given the `counts` of the system
    function expected_size(spec, counts)
        spec[2] isa Symbol || return spec # a fixed size, e.g. a stiffness tensor
        d = spec[1]
        n = counts[spec[2]]
        return d == 1 ? (n,) : (d, n)
    end

    # number of points, local points, bonds, local dof and max. neighbors of `system`
    function counts(system::Peridynamics.AbstractBondSystem)
        return Dict(:all => Peridynamics.get_n_points(system),
                    :loc => Peridynamics.get_n_loc_points(system),
                    :bonds => Peridynamics.get_n_bonds(system),
                    :dof => Peridynamics.get_n_loc_dof(system),
                    :neighbors => maximum(system.n_neighbors))
    end
    function counts(system::Peridynamics.InteractionSystem)
        return Dict(:all => Peridynamics.get_n_points(system),
                    :loc => Peridynamics.get_n_loc_points(system),
                    :bonds => Peridynamics.get_n_one_nis(system),
                    :dof => Peridynamics.get_n_loc_dof(system))
    end

    # the first of two chunks of a ten-point line body of `mat`, with its storage and system
    function first_chunk(mat, solver)
        position = zeros(3, 10)
        position[1, :] .= 0.0:9.0
        body = Body(mat, position, ones(10))
        material!(body, horizon=1.5, rho=1.0, E=1.0, nu=0.25, Gc=1.0)
        pd = Peridynamics.PointDecomposition(body, 2)
        ps = Peridynamics.get_param_spec(body)
        return Peridynamics.BodyChunk(body, solver, pd, 1, ps)
    end

    # check every field of the storage of `mat` with `solver` against the table
    function check(mat, solver)
        chunk = first_chunk(mat, solver)
        (; storage, system) = chunk
        n = counts(system)
        expected = merge(COMMON_FIELDS, SOLVER_FIELDS[typeof(solver)],
                         Dict(MATERIAL_FIELDS)[mat])
        system isa Peridynamics.AbstractBondSystem && merge!(expected, BOND_SYSTEM_FIELDS)
        expected[:b_int] = b_int_spec(mat, solver)
        expected[:velocity_half] = velocity_half_spec(mat, solver)
        # the rotated formulations are not supported by the Newton-Krylov solver and have no
        # fields for it
        if mat isa Union{CRMaterial,RKCRMaterial}
            for field in (:residual, :displacement_copy, :b_int_copy, :temp_force, :Δu,
                          :v_temp, :Jv_temp)
                delete!(expected, field)
            end
        end
        # nothing is forgotten: every field of the storage has an expectation and vice versa
        @test Set(fieldnames(typeof(storage))) == Set(keys(expected))
        for (field, spec) in expected
            hasfield(typeof(storage), field) || continue
            value = getfield(storage, field)
            if spec[1] === :state
                # the nested state of a stateless constitutive model is `nothing`
                @test isnothing(value) || (field, typeof(value)) === nothing
            elseif spec[1] == 0
                @test isempty(value) || (field, size(value)) === nothing
            else
                @test size(value) == expected_size(spec, n) || (field, size(value)) === nothing
            end
        end
        return nothing
    end
end

@testitem "storage fields: sizes for every material and solver" setup=[Fixtures, StorageFields] begin
    for (mat, _) in StorageFields.MATERIAL_FIELDS
        @testset "$(nameof(typeof(mat)))" begin
            for solver in (VelocityVerlet(steps=1), DynamicRelaxation(steps=1), NewtonKrylov(steps=1))
                Fixtures.supports(solver, mat) || continue
                @testset "$(nameof(typeof(solver)))" begin
                    StorageFields.check(mat, solver)
                end
            end
        end
    end
end
