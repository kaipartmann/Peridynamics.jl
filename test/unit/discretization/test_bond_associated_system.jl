@testitem "bond-associated family" begin
    # setup
    position = [0.0 1.0 0.0 0.0
                0.0 0.0 1.0 0.0
                0.0 0.0 0.0 1.0]
    volume = [1.1, 1.2, 1.3, 1.4]
    mat = BACMaterial()
    body = Body(mat, position, volume)
    material!(body, horizon=2, rho=1, E=1, nu=0.25, Gc=1)
    pd = Peridynamics.PointDecomposition(body, 1)

    # 1
    system = Peridynamics.BondAssociatedSystem(body, pd, 1)

    @test system.position == position
    @test system.volume == volume
    @test system.neighbor == [2, 3, 4, 1, 3, 4, 1, 2, 4, 1, 2, 3]
    @test system.bond_length == [1.0, 1.0, 1.0, 1.0, √2, √2, 1.0, √2, √2, 1.0, √2, √2]
    @test system.fail_permit == [true, true, true, true, true, true, true, true, true, true, true, true]
    @test system.n_neighbors == [3, 3, 3, 3]
    @test system.bond_ids == [1:3, 4:6, 7:9, 10:12]

end

@testitem "intersection bond ids" begin
    pos, vol = uniform_box(1, 0.25, 0.25, 0.25)
    body = Body(BACMaterial(), pos, vol)
    material!(body, horizon=0.5, bond_horizon=0.5, rho=1, E=1, nu=0.25, Gc=1)

    pd = Peridynamics.PointDecomposition(body, 1)

    bas = Peridynamics.get_system(body, pd, 1)
    (; neighbor, bond_length, fail_permit, bond_ids, intersection_bond_ids,
       intersection_ids) = bas

    @test neighbor == [2, 3, 1, 3, 4, 1, 2, 4, 2, 3]
    @test bond_length == [0.25, 0.5, 0.25, 0.25, 0.5, 0.5, 0.25, 0.25, 0.5, 0.25]
    @test fail_permit == fill(true, 10)
    @test bond_ids == [1:2, 3:5, 6:8, 9:10]

    # the families are stored flat: `intersection_bond_ids` holds the bond indices of the
    # chunk of all families one after the other, `intersection_ids[bond_id]` is the range of
    # that vector which belongs to `bond_id`
    @test intersection_bond_ids == [1, 2, 1, 2, 3, 4, 5, 4, 5, 6, 7, 6, 7, 8, 9, 10, 9, 10]
    @test intersection_ids == [1:2, 3:4, 5:5, 6:7, 8:9, 10:11, 12:13, 14:14, 15:16, 17:18]

    # the family of a bond is a family of the bonds of its own point
    for i in Peridynamics.each_point_idx(bas)
        for bond_idx in Peridynamics.each_bond_idx(bas, i)
            family = Peridynamics.each_intersecting_bond_idx(bas, i, bond_idx)
            @test issubset(family, Peridynamics.each_bond_idx(bas, i))
            @test in(bond_idx, family)
        end
    end
end

@testitem "bond-associated linear momentum consistency" begin
    # neither is a test dependency, so both come through the package
    using Peridynamics: SMatrix
    using Peridynamics.LinearAlgebra: I, inv, opnorm

    # A homogeneous deformation makes every `P_b` the same `P`, and the force state then gives
    # `Σ_j T_ij ⊗ ΔX_ij V_j = P Σ_b w_b`. So the weights have to sum to one and the operator
    # applied to `P` has to be the identity, also at the surface with incomplete families.
    # the identities are exact per point and hold at any resolution, so the bodies are small;
    # the third case keeps five points per edge so that the bond horizon gives non-trivial
    # intersections
    for (npyz, m, bond_horizon) in ((4, 3.015, nothing), (4, 2.015, nothing),
                                    (5, 3.015, 3.015 / 5 * 1.5))
        Δx = 1.0 / npyz
        δ = m * Δx
        pos, vol = uniform_box(1, 1, 1, Δx)
        body = Body(BACMaterial(), pos, vol)
        if isnothing(bond_horizon)
            material!(body; horizon=δ, rho=1, E=1, nu=0.25, Gc=1)
        else
            material!(body; horizon=δ, bond_horizon, rho=1, E=1, nu=0.25, Gc=1)
        end
        velocity_bc!(t -> 0.0, body, :all_points, :x)
        solver = VelocityVerlet(steps=1)
        dh = Peridynamics.threads_data_handler(body, solver, 1)
        Peridynamics.init_time_solver!(solver, dh)
        Peridynamics.initialize!(dh, solver)
        chunk = dh.chunks[1]
        system = chunk.system
        (; volume) = system

        worst_weights, worst_operator = 0.0, 0.0
        for i in 1:Peridynamics.get_n_loc_points(chunk)
            weights = 0.0
            A = zero(SMatrix{3,3,Float64,9})
            for bond_idx in Peridynamics.each_bond_idx(system, i)
                K = zero(SMatrix{3,3,Float64,9})
                for bond_id in Peridynamics.each_intersecting_bond_idx(system, i, bond_idx)
                    jj = Peridynamics.get_neighbor(system, bond_id)
                    ΔX = Peridynamics.get_vector_diff(system.position, i, jj, Peridynamics.dims(system))
                    K += Peridynamics.kernel(system, bond_id) * volume[jj] * (ΔX * ΔX')
                end
                Kinv = inv(K)
                w = Peridynamics.volume_fraction_factor(system, i, bond_idx)
                weights += w
                for bond_id in Peridynamics.each_intersecting_bond_idx(system, i, bond_idx)
                    j = Peridynamics.get_neighbor(system, bond_id)
                    ΔX = Peridynamics.get_vector_diff(system.position, i, j, Peridynamics.dims(system))
                    ω = Peridynamics.kernel(system, bond_id)
                    A += w * ω * volume[j] * (Kinv * (ΔX * ΔX'))
                end
            end
            worst_weights = max(worst_weights, abs(weights - 1))
            worst_operator = max(worst_operator, opnorm(A - I))
        end
        @test worst_weights < 1e-12
        @test worst_operator < 1e-10
    end
end

@testitem "bond-associated force state is the gradient of the energy" begin
    # neither is a test dependency, so both come through the package
    using Peridynamics: SMatrix
    using Peridynamics.LinearAlgebra: inv

    # The force state has to be `-∂/∂x_k Σ_i V_i W_i` with `W_i = Σ_b w_b W(F_b)`, otherwise no
    # energy is conserved. Compares the force density against a central difference of the energy,
    # so this covers the force loop and not only the weights.
    npyz = 4
    Δx = 1.0 / npyz
    pos, vol = uniform_box(1, 1, 1, Δx)
    mat = BACMaterial()
    body = Body(mat, pos, vol)
    material!(body; horizon=3.015Δx, rho=7850.0, E=210e9, nu=0.25, Gc=2.7)
    no_failure!(body)
    velocity_bc!(t -> 0.0, body, :all_points, :x)
    solver = VelocityVerlet(steps=1)
    dh = Peridynamics.threads_data_handler(body, solver, 1)
    Peridynamics.init_time_solver!(solver, dh)
    Peridynamics.initialize!(dh, solver)
    chunk = dh.chunks[1]
    system, storage = chunk.system, chunk.storage
    params = Peridynamics.get_params(chunk, 1)

    # total strain energy `Σ_i V_i Σ_b w_b W(F_b)` of the current configuration
    function total_energy()
        total = 0.0
        for i in Peridynamics.each_point_idx(system)
            for bond_idx in Peridynamics.each_bond_idx(system, i)
                K = zero(SMatrix{3,3,Float64,9})
                _F = zero(SMatrix{3,3,Float64,9})
                for bond_id in Peridynamics.each_intersecting_bond_idx(system, i, bond_idx)
                    j = Peridynamics.get_neighbor(system, bond_id)
                    ΔX = Peridynamics.get_vector_diff(system.position, i, j, Peridynamics.dims(system))
                    Δx_ = Peridynamics.get_vector_diff(storage.position, i, j, Peridynamics.dims(system))
                    ωV = Peridynamics.kernel(system, bond_id) * system.volume[j]
                    K += ωV * (ΔX * ΔX')
                    _F += ωV * (Δx_ * ΔX')
                end
                F = _F * inv(K)
                w = Peridynamics.volume_fraction_factor(system, i, bond_idx)
                Ψ = Peridynamics.strain_energy_density(mat.constitutive_model, storage,
                                                       params, F)
                total += system.volume[i] * w * Ψ
            end
        end
        return total
    end

    # not the identity and with a perturbation, so that no symmetry can hide a wrong force
    for i in axes(system.position, 2)
        X1, X2, X3 = system.position[1, i], system.position[2, i], system.position[3, i]
        storage.position[1, i] = 1.02X1 + 0.004X2 + 1e-3 * sinpi(2X2)
        storage.position[2, i] = 0.99X2 + 1e-3 * sinpi(2X3)
        storage.position[3, i] = 1.01X3 + 1e-3 * sinpi(2X1)
    end
    Peridynamics.calc_force_density!(chunk, 0.0, solver.Δt)
    b_int = copy(storage.b_int)

    # compromise between truncation error and cancellation in the energy (which is ~1e5 here)
    h = 1e-9
    scale = maximum(abs, b_int)
    errors = [begin
                  x0 = storage.position[d, i]
                  storage.position[d, i] = x0 + h
                  Ep = total_energy()
                  storage.position[d, i] = x0 - h
                  Em = total_energy()
                  storage.position[d, i] = x0
                  numerical = -(Ep - Em) / (2h) / system.volume[i]
                  abs(numerical - b_int[d, i]) / scale
              end
              for i in (1, 22, 43, 64), d in 1:3] # two corners, two interior points
    @test maximum(errors) < 1e-6
end

@testitem "bond horizon" begin
    # setup
    pos, vol = uniform_box(1, 0.25, 0.25, 0.25)
    body = Body(BACMaterial(), pos, vol)
    @test_throws ArgumentError material!(body, horizon=0.5, bond_horizon=-1, rho=1, E=1,
        nu=0.25, Gc=1)
    @test_logs (:warn,) material!(body, horizon=0.5, bond_horizon=0.1, rho=1, E=1, nu=0.25,
        Gc=1)
end

@testitem "bond-associated compatibility" begin
    # setup
    pos, vol = uniform_box(1, 0.25, 0.25, 0.25)
    body = Body(BBMaterial(), pos, vol)
    @test_throws ArgumentError Peridynamics.check_system_compat(Peridynamics.BondAssociatedSystem,
                                                                body.mat)
end

@testitem "bond-associated required parameters" begin
    params = Peridynamics.required_point_parameters(BACMaterial)
    @test :δ in params && :δb in params && :rho in params && :E in params
end
