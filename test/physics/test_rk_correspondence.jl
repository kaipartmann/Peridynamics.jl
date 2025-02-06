# Write a unit test for the qmatrix calculation of the RKMaterial
@testitem "Monomial vector Q" begin
    using Peridynamics.StaticArrays

    # CPP code as reference:
    # // calculate dimension of Q vector
    # Qdim = 0;
    # for(thisOrder=1; thisOrder<=accuracyOrder; thisOrder++){
    #   for(p1=thisOrder; p1>=0; p1--){ // x-power
    #     for(p2=thisOrder-p1; p2>=0; p2--){ // y-power
    #       p3=thisOrder-p1-p2; //z-power
    #       Qdim++;
    #     }
    #   }
    # }
    # // Calculate Q for this bond
    # counter = 0;
    # for(thisOrder=1; thisOrder<=accuracyOrder; thisOrder++){
    #   for(p1=thisOrder; p1>=0; p1--){ // x-power
    #     for(p2=thisOrder-p1; p2>=0; p2--){ // y-power
    #       p3=thisOrder-p1-p2; //z-power

    #       Q[counter] = 1.0;
    #       for(i=0; i<p1; i++)
    #         Q[counter] *= undeformedBondX / *delta;
    #       for(i=0; i<p2; i++)
    #         Q[counter] *= undeformedBondY / *delta;
    #       for(i=0; i<p3; i++)
    #         Q[counter] *= undeformedBondZ / *delta;

    #       counter++;
    #     }
    #   }
    # }
    # Julia function as test:
    function q_vector(ΔX, δ, m=2)
        # Calculate dimension of Q vector
        q_dim = 0
        for thisOrder in 1:m
            for p1 in thisOrder:-1:0
                for p2 in thisOrder-p1:-1:0
                    p3 = thisOrder - p1 - p2
                    q_dim += 1
                end
            end
        end
        # Calculate Q for this bond
        Q = zeros(q_dim)
        counter = 1
        for thisOrder in 1:m
            for p1 in thisOrder:-1:0
                for p2 in thisOrder-p1:-1:0
                    p3 = thisOrder - p1 - p2
                    Q[counter] = 1.0
                    for i in 1:p1
                        Q[counter] *= ΔX[1] / δ
                    end
                    for i in 1:p2
                        Q[counter] *= ΔX[2] / δ
                    end
                    for i in 1:p3
                        Q[counter] *= ΔX[3] / δ
                    end
                    counter += 1
                end
            end
        end
        return Q
    end

    # perilab function as a test:
    function calculate_Q(
        accuracy_order::Int64,
        dof::Int64,
        bond_geometry::Vector{Float64},
        horizon::Union{Int64,Float64},
        Q::Vector{Float64},
    )
        counter = 0
        p = @MVector zeros(Int64, dof)
        for this_order = 1:accuracy_order
            for p[1] = this_order:-1:0
                if dof == 3
                    for p[2] = this_order-p[1]:-1:0
                        p[3] = this_order - p[1] - p[2]
                        # Calculate the product for Q[counter]
                        counter += 1
                        Q[counter] = prod_Q(bond_geometry, horizon, p, Q[counter])
                    end
                else
                    p[2] = this_order - p[1]
                    counter += 1
                    Q[counter] = prod_Q(bond_geometry, horizon, p, Q[counter])
                end
            end
        end
        return Q
    end

    function prod_Q(bond_geometry, horizon, p, Q)
        Q = 1
        @inbounds @fastmath for m ∈ axes(p, 1)
            @views @inbounds @fastmath for pc = 1:p[m]
                Q *= bond_geometry[m] / horizon
            end
        end
        return Q
    end



    # test cases for q_vector
    undeformed_bond_and_horizon = [
        ([0.0, 0.0, 0.0], 1.0),
        ([1.0, 1.0, 1.0], 1.0),
        ([1.0, 0.0, 0.0], 2.0),
        (rand(3), rand()),
    ]
    for m in (1,2,3)
        for (ΔX, δ) in undeformed_bond_and_horizon
            Q_peridigm = q_vector(ΔX, δ, m)
            Q_perilab = calculate_Q(m, 3, ΔX, δ, zeros(size(Q_peridigm)))
            Q_peridynamics = Peridynamics.get_monomial_vector(Val(m), ΔX, δ)
            @test Q_peridigm ≈ Q_perilab
            @test Q_peridigm ≈ Q_peridynamics
            @test Q_perilab ≈ Q_peridynamics
        end
    end
end

@testitem "RKCMaterial kernel" begin
    mat = RKCMaterial()
    δ = 1.0
    kwargs = Dict{Symbol,Any}(:horizon => δ, :rho => 1, :E => 1, :nu => 0.25, :Gc => 1.0)
    params = Peridynamics.RKCPointParameters(mat, kwargs)

    L = 0.0
    ω_a = 0.0
    @test Peridynamics.get_kernel(mat, params, L) ≈ ω_a

    L = 1/4
    ω_a = 2/3 - 4 * (L/δ)^2 + 4 * (L/δ)^3
    @test Peridynamics.get_kernel(mat, params, L) ≈ ω_a

    L = 1/2
    ω_a = 2/3 - 4 * (L/δ)^2 + 4 * (L/δ)^3
    @test Peridynamics.get_kernel(mat, params, L) ≈ ω_a

    L = 3/4
    ω_a = 4/3 - 4 * (L/δ) + 4 * (L/δ)^2 - 4/3 * (L/δ)^3
    @test Peridynamics.get_kernel(mat, params, L) ≈ ω_a

    L = 1
    ω_a = 4/3 - 4 * (L/δ) + 4 * (L/δ)^2 - 4/3 * (L/δ)^3
    @test Peridynamics.get_kernel(mat, params, L) ≈ ω_a

    @test Peridynamics.get_kernel(mat, params, 1.5) ≈ 0.0
    @test Peridynamics.get_kernel(mat, params, -1.0) ≈ 0.0
end

@testitem "Gradient weights Φ" begin
    using Peridynamics, Test
    using Peridynamics: get_diff, get_monomial_vector, each_bond_idx, influence_function,
                        threads_data_handler, Bond
    using Peridynamics.StaticArrays
    # setup
    ref_position = [0.0 1.0 0.0
                    0.0 0.0 1.0
                    0.0 0.0 0.0]
    volume = fill(1.0, 3)
    δ = 1.5
    mat = RKCMaterial()
    body = Body(mat, ref_position, volume)
    material!(body, horizon=δ, rho=1, E=1, nu=0.25, Gc=1.0)
    failure_permit!(body, false)
    dh = threads_data_handler(body, VelocityVerlet(steps=1), 1)
    chunk = dh.chunks[1]
    (; mat, storage, system, paramsetup) = chunk

    @test system.bonds == [
        Bond(2, 1.0, false), # point 1
        Bond(3, 1.0, false),
        Bond(1, 1.0, false), # point 2
        Bond(3, √2, false),
        Bond(1, 1.0, false), # point 3
        Bond(2, √2, false),
    ]

    @test storage.gradient_weight == zeros(3, 6)

    Peridynamics.calc_gradient_weights!(storage, system, mat, paramsetup)

    # my old version with temp = ωij * volume[j]
    # Φ = storage.gradient_weight
    # @test Φ[:,1] ≈ [0.5000000000000001, 0.0, 0.0]
    # @test Φ[:,2] ≈ [0.0, 0.5000000000000001, 0.0]
    # @test Φ[:,3] ≈ [-0.5000000000000004, -0.3333333333333337, 0.0]
    # @test Φ[:,4] ≈ [-3.469446951953614e-18, 0.33333333333333365, 0.0]
    # @test Φ[:,5] ≈ [-0.33333333333333054, -0.5000000000000009, 0.0]
    # @test Φ[:,6] ≈ [0.3333333333333338, 7.728193085476676e-16, 0.0]

    # PeriLab version temp = ωij / δ
    # Φ = storage.gradient_weight
    # @test Φ[:,1] ≈ [0.0006958021262525432, 0.0, 0.0]
    # @test Φ[:,2] ≈ [0.0, 0.0006958021262525432, 0.0]
    # @test Φ[:,3] ≈ [-0.0006993163680475011, 3.5142417949579743e-6, 0.0]
    # @test Φ[:,4] ≈ [-3.5552013688439347e-6, 4.095957388592146e-8, 0.0]
    # @test Φ[:,5] ≈ [3.514241794958013e-6, -0.0006993163680475009, 0.0]
    # @test Φ[:,6] ≈ [4.095957388592165e-8, -3.555201368843933e-6, 0.0]

    #=
    using BenchmarkTools
    @btime Peridynamics.calc_gradient_weights!($storage, $system, $mat, $paramsetup)
    =#

    # storage.position[:, 1] += [-0.001, -0.001, -0.001]
    # storage.position[:, 2] += [0.001, 0.0, 0.0]
    # storage.position[:, 3] += [0.0, 0.001, 0.0]

    # F = storage.defgrad
    # @test F == zeros(9, 3)

    # Peridynamics.calc_deformation_gradients!(storage, system, mat, paramsetup, 0.0, 0.0)

    # my old version
    # @test F[:, 1] ≈ [1.0000013916042525, 6.958021262525431e-7, 0.0,
    #                  6.958021262525431e-7, 1.0000013916042525, 0.0,
    #                  6.958021262525431e-7, 6.958021262525431e-7, 1.0]
    # @test F[:, 2] ≈ [1.0000014021879375, -7.069443163801091e-9, 0.0,
    #                  6.957611666786575e-7, 0.9999999965267178, 0.0,
    #                  6.993163680475011e-7, -3.514241794957974e-9, 1.0]
    # @test F[:, 3] ≈ [0.9999999965267178, 6.957611666786573e-7, 0.0,
    #                  -7.069443163801169e-9, 1.0000014021879375, 0.0,
    #                  -3.514241794958013e-9, 6.993163680475009e-7, 1.0]

    # PeriLab version
    # @test F[:, 1] ≈ [1.0000013916042525, 6.958021262525431e-7, 6.958021262525431e-7,
    #                  6.958021262525431e-7, 1.0000013916042525, 6.958021262525431e-7,
    #                  0.0, 0.0, 1.0]
    # @test F[:, 2] ≈ [1.0000014021879375, 6.957611666786575e-7, 6.993163680475011e-7,
    #                  -7.069443163801091e-9, 0.9999999965267178, -3.514241794957974e-9,
    #                  0.0, 0.0, 1.0]
    # @test F[:, 3] ≈ [0.9999999965267178, -7.069443163801169e-9, -3.514241794958013e-9,
    #                  6.957611666786573e-7, 1.0000014021879375, 6.993163680475009e-7,
    #                  0.0, 0.0, 1.0]


    #=
    using BenchmarkTools
    @btime Peridynamics.calc_deformation_gradients!($storage, $system, $mat, $paramsetup, $0.0, $0.0)
    =#


    #--- My own testing

    # ξ12 = 1/1.5
    # ω12 = ω13 = 4/3 - 4 * ξ12 + 4 * ξ12^2 - 4/3 * ξ12^3
    # ΔX12 = SVector{3}([1.0, 0.0, 0.0])
    # ΔX13 = SVector{3}([0.0, 1.0, 0.0])
    # Q12 = SVector{9}([1.0, 0.0, 0.0, 1, 0.0, 0.0, 0.0, 0.0, 0.0])
    # # Q12 = copy(ΔX12)
    # Q13 = SVector{9}([0.0, 1.0, 0.0, 0.0, 1, 0.0, 0.0, 0.0, 0.0])
    # # Q13 = copy(ΔX13)
    # V1 = V2 = V3 = 1.0
    # M1 = ω12 * (Q12 * Q12') * V2 + ω13 * (Q13 * Q13') * V3

    # threshold = 1e-6 * δ^3
    # M1inv = Peridynamics.invert_moment_matrix(M1, threshold)

    # Q∇ᵀ = SMatrix{3,9,Int,27}(1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    #                           0, 0, 0, 0, 0, 0)
    # # Q∇ᵀ = SMatrix{3,3,Int,9}(1, 0, 0, 0, 1, 0, 0, 0, 1)


    # Φ12 = ω12 * Q∇ᵀ * M1inv * Q12
    # Φ13 = ω13 * Q∇ᵀ * M1inv * Q13
    # @test Φ[:,1] ≈ Φ12
    # @test Φ[:,2] ≈ Φ13

end

@testitem "Damage changed" begin
    # setup
    ref_position = [0.0 1.0 0.0
                    0.0 0.0 1.0
                    0.0 0.0 0.0]
    volume = fill(1.0, 3)
    δ = 1.5
    mat = RKCMaterial()
    body = Body(mat, ref_position, volume)
    material!(body, horizon=δ, rho=1, E=1, nu=0.25, Gc=1.0)
    failure_permit!(body, false)
    dh = Peridynamics.threads_data_handler(body, VelocityVerlet(steps=1), 1)
    chunk = dh.chunks[1]
    (; mat, storage, system, paramsetup) = chunk

    # initial value should be `true`
    @test storage.damage_changed == [true, true, true]

    # after first damage calculation this should be `false`
    Peridynamics.calc_damage!(chunk)
    @test storage.damage_changed == [false, false, false]

    # damage a bond:
    @test storage.bond_active[1] == true
    storage.bond_active[1] = false
    @test storage.bond_active[1] == false
    @test storage.n_active_bonds[1] == 2
    storage.n_active_bonds[1] -= 1
    @test storage.n_active_bonds[1] == 1

    Peridynamics.calc_damage!(chunk)
    @test storage.damage_changed == [true, false, false]
end

@testitem "accuracy_order = 1" begin
    # setup
    ref_position = [0.0 1.0 0.0
                    0.0 0.0 1.0
                    0.0 0.0 0.0]
    volume = fill(1.0, 3)
    δ = 1.5
    mat_rk = RKCMaterial(accuracy_order=1)
    body_kr = Body(mat_rk, ref_position, volume)
    material!(body_kr, horizon=δ, rho=1, E=1, nu=0.25, Gc=1.0)
    failure_permit!(body_kr, false)
    dh_rk = Peridynamics.threads_data_handler(body_kr, VelocityVerlet(steps=1), 1)
    chunk_rk = dh_rk.chunks[1]
    body_c = Body(CMaterial(), ref_position, volume)
    material!(body_c, horizon=δ, rho=1, E=1, nu=0.25, Gc=1.0)
    failure_permit!(body_c, false)
    dh_c = Peridynamics.threads_data_handler(body_c, VelocityVerlet(steps=1), 1)
    chunk_c = dh_c.chunks[1]

    (; mat, storage, system, paramsetup) = chunk_rk
    Peridynamics.calc_gradient_weights!(storage, system, mat, paramsetup)

end


##----

Φ = rand(3)
ΔU = rand(3)
F1 = ΔU * Φ'
F2 = zeros(3,3)
for i in 1:3
    F2[i,1] = ΔU[i] * Φ[1]
    F2[i,2] = ΔU[i] * Φ[2]
    F2[i,3] = ΔU[i] * Φ[3]
end
F1 ≈ F2
