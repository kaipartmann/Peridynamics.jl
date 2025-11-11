@testitem "get and update tensors" begin
    using Peridynamics.StaticArrays

    A = zeros(9, 2)
    for i in eachindex(A)
        A[i] = i
    end
    B = Peridynamics.get_tensor(A, 1)
    C = Peridynamics.get_tensor(A, 2)
    for i in eachindex(B,C)
        @test B[i] ≈ A[i, 1]
        @test C[i] ≈ A[i, 2]
    end

    Peridynamics.update_tensor!(A, 2, reverse(B))
    Peridynamics.update_tensor!(A, 1, reverse(C))
    @test A[:, 1] ≈ [18:-1:10;]
    @test A[:, 2] ≈ [9:-1:1;]
end

@testitem "get and update vectors" begin
    using Peridynamics.StaticArrays

    a = zeros(3, 2)
    for i in eachindex(a)
        a[i] = i
    end
    @test Peridynamics.get_vector(a, 1) ≈ [1, 2, 3]
    @test Peridynamics.get_vector(a, 2) ≈ [4, 5, 6]

    Peridynamics.update_vector!(a, 1, SVector{3}(3.0, 2.0, 1.0))
    @test a[:, 1] ≈ [3, 2, 1]

    Peridynamics.update_vector!(a, 2, SVector{3}(6.0, 5.0, 4.0))
    @test a[:, 2] ≈ [6, 5, 4]

    Peridynamics.update_add_vector!(a, 1, SVector{3}(1.0, 1.0, 1.0))
    Peridynamics.update_add_vector!(a, 2, SVector{3}(1.0, 1.0, 1.0))
    @test a ≈ [4; 3; 2;; 7; 6; 5]

    @test Peridynamics.get_vector_diff(a, 1, 2) ≈ [3, 3, 3]
end

@testitem "invreg" begin
    using Peridynamics.StaticArrays
    using Peridynamics.LinearAlgebra

    # Test with well-conditioned matrices
    A = @SMatrix rand(3, 3)
    @test Peridynamics.invreg(A, 0, 0) * A ≈ I
    @test Peridynamics.invreg(A, 0, sqrt(eps())) * A ≈ I
    @test Peridynamics.invreg(A, 1e-10, sqrt(eps())) * A ≈ I

    A = @SMatrix rand(4, 4)
    @test Peridynamics.invreg(A, 0, 0) * A ≈ I
    @test Peridynamics.invreg(A, 0, sqrt(eps())) * A ≈ I
    @test Peridynamics.invreg(A, 1e-10, sqrt(eps())) * A ≈ I

    A = @SMatrix rand(5, 5)
    @test Peridynamics.invreg(A, 0, 0) * A ≈ I
    @test Peridynamics.invreg(A, 0, sqrt(eps())) * A ≈ I
    @test Peridynamics.invreg(A, 1e-10, sqrt(eps())) * A ≈ I

    A = @SMatrix rand(6, 6)
    @test Peridynamics.invreg(A, 0, 0) * A ≈ I
    @test Peridynamics.invreg(A, 0, sqrt(eps())) * A ≈ I
    @test Peridynamics.invreg(A, 1e-10, sqrt(eps())) * A ≈ I

    B = @SMatrix [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    @test Peridynamics.invreg(B, 0, 0) ≈ I
    @test Peridynamics.invreg(B, 0, sqrt(eps())) ≈ I
    @test Peridynamics.invreg(B, 1e-10, sqrt(eps())) ≈ I

    # Test with moderately ill-conditioned matrices
    A = @SMatrix [1.0 0.99 0.0; 0.99 1.0 0.0; 0.0 0.0 1.0]
    @test Peridynamics.invreg(A, 0, 1e-8) * A ≈ I
    @test Peridynamics.invreg(A, 1e-8, 1e-8) * A ≈ I

    A = @SMatrix [2.0 1.9 0.0 0.0; 1.9 2.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]
    @test Peridynamics.invreg(A, 0, 1e-8) * A ≈ I
    @test Peridynamics.invreg(A, 1e-8, 1e-8) * A ≈ I

    A = @SMatrix [1.0 0.5 0.3; 0.5 1.0 0.4; 0.3 0.4 1.0]
    @test Peridynamics.invreg(A, 0, 1e-8) * A ≈ I
    @test Peridynamics.invreg(A, 1e-8, 1e-8) * A ≈ I

    # Test with very ill-conditioned matrices
    A = @SMatrix [1.0 1.0 1.0; 1.0 1.0+1e-10 1.0; 1.0 1.0 1.0+1e-10]
    @test Peridynamics.containsnan(inv(A)) # the normal inverse contains NaNs
    # the regularized inverse should not contain NaNs
    @test !Peridynamics.containsnan(Peridynamics.invreg(A, 0, 0))
    @test !Peridynamics.containsnan(Peridynamics.invreg(A, 1e-8, 1e-8))

    # Test with larger lambda values
    # Do not test with random matrices as the errors can get very large
    A = @SMatrix [0.5872 0.8188 0.792;
                  0.5830 0.9880 0.202;
                  0.1727 0.1524 0.354]
    @test Peridynamics.invreg(A, 1e-4, sqrt(eps())) * A ≈ I atol=√(2e-4)
    A = @SMatrix [0.3729 0.3725 0.3804 0.7302;
                  0.9697 0.4572 0.2289 0.2289;
                  0.5170 0.8261 0.3050 0.0196;
                  0.1578 0.5413 0.3667 0.3920]
    @test Peridynamics.invreg(A, 1e-5, sqrt(eps())) * A ≈ I atol=√(2e-5)
end
