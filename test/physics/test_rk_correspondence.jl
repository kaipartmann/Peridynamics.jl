# Write a unit test for the qmatrix calculation of the RKMaterial
@testitem "Monomial vector Q" begin
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
            Q_peridynamics = Peridynamics.get_monomial_vector(Val(m), ΔX, δ)
            @test Q_peridigm ≈ Q_peridynamics
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

    # @test system.bonds == [
    #     Bond(2, 1.0, false), # point 1
    #     Bond(3, 1.0, false),
    #     Bond(1, 1.0, false), # point 2
    #     Bond(3, √2, false),
    #     Bond(1, 1.0, false), # point 3
    #     Bond(2, √2, false),
    # ]

    # @test storage.gradient_weight == zeros(3, 6)

    Peridynamics.calc_gradient_weights!(storage, system, mat, paramsetup)

    # my old version with temp = ωij * volume[j]
    # Φ = storage.gradient_weight
    # @test Φ[:,1] ≈ [0.0010437031893788147, 0.0, 0.0]
    # @test Φ[:,2] ≈ [0.0, 0.0010437031893788147, 0.0]
    # @test Φ[:,3] ≈ [-0.0010489745520712516, 5.271362692436961e-6, 0.0]
    # @test Φ[:,4] ≈ [-5.332802053265902e-6, 6.143936082888219e-8, 0.0]
    # @test Φ[:,5] ≈ [5.27136269243702e-6, -0.0010489745520712514, 0.0]
    # @test Φ[:,6] ≈ [6.143936082888247e-8, -5.332802053265901e-6, 0.0]

    # PeriLab version temp = ωij / δ
    Φ = storage.gradient_weight
    @test Φ[:,1] ≈ [0.0006958021262525432, 0.0, 0.0]
    @test Φ[:,2] ≈ [0.0, 0.0006958021262525432, 0.0]
    @test Φ[:,3] ≈ [-0.0006993163680475011, 3.5142417949579743e-6, 0.0]
    @test Φ[:,4] ≈ [-3.5552013688439347e-6, 4.095957388592146e-8, 0.0]
    @test Φ[:,5] ≈ [3.514241794958013e-6, -0.0006993163680475009, 0.0]
    @test Φ[:,6] ≈ [4.095957388592165e-8, -3.555201368843933e-6, 0.0]

    #=
    using BenchmarkTools
    @btime Peridynamics.calc_gradient_weights!($storage, $system, $mat, $paramsetup)
    =#

    storage.position[:, 1] += [-0.001, -0.001, -0.001]
    storage.position[:, 2] += [0.001, 0.0, 0.0]
    storage.position[:, 3] += [0.0, 0.001, 0.0]

    F = storage.defgrad
    @test F == zeros(9, 3)

    Peridynamics.calc_deformation_gradients!(storage, system, mat, paramsetup, 0.0, 0.0)

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
    @test F[:, 1] ≈ [1.0000013916042525, 6.958021262525431e-7, 6.958021262525431e-7,
                     6.958021262525431e-7, 1.0000013916042525, 6.958021262525431e-7,
                     0.0, 0.0, 1.0]
    @test F[:, 2] ≈ [1.0000014021879375, 6.957611666786575e-7, 6.993163680475011e-7,
                     -7.069443163801091e-9, 0.9999999965267178, -3.514241794957974e-9,
                     0.0, 0.0, 1.0]
    @test F[:, 3] ≈ [0.9999999965267178, -7.069443163801169e-9, -3.514241794958013e-9,
                     6.957611666786573e-7, 1.0000014021879375, 6.993163680475009e-7,
                     0.0, 0.0, 1.0]


    #=
    using BenchmarkTools
    @btime Peridynamics.calc_deformation_gradients!($storage, $system, $mat, $paramsetup, $0.0, $0.0)
    =#
end
