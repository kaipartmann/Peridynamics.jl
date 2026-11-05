@testitem "get_sym_tensor / update_sym_tensor!: Voigt round trip of a symmetric tensor" begin
    using Peridynamics.StaticArrays
    using Peridynamics.LinearAlgebra

    A = SMatrix{3,3,Float64,9}(1.0, 4.0, 5.0,
                               4.0, 2.0, 6.0,
                               5.0, 6.0, 3.0)

    M = zeros(6, 3)
    Peridynamics.update_sym_tensor!(M, 2, A)

    # Voigt order: (11, 22, 33, 23, 13, 12)
    @test M[:, 2] == [1.0, 2.0, 3.0, 6.0, 5.0, 4.0]
    # the neighbouring columns must be untouched
    @test all(iszero, M[:, 1])
    @test all(iszero, M[:, 3])

    # the round trip is exact and gives back a symmetric tensor
    B = Peridynamics.get_sym_tensor(M, 2)
    @test B == A
    @test B == B'

    # only the upper triangle is read and mirrored back, so an asymmetry of round-off size
    # survives the round trip as round-off and never leaks into another component
    Aperturbed = A + 1e-16 * SMatrix{3,3,Float64,9}(0.0, -1.0, 0.0, 1.0, 0.0, 0.0,
                                                    0.0, 0.0, 0.0)
    Peridynamics.update_sym_tensor!(M, 2, Aperturbed)
    @test Peridynamics.get_sym_tensor(M, 2) ≈ A atol=1e-15

    # a tensor with a real skew part is not symmetrized: the upper triangle wins
    Askew = A + SMatrix{3,3,Float64,9}(0.0, -0.5, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0)
    Peridynamics.update_sym_tensor!(M, 2, Askew)
    @test Peridynamics.get_sym_tensor(M, 2)[1, 2] == Askew[1, 2]
    @test Peridynamics.get_sym_tensor(M, 2)[2, 1] == Askew[1, 2]

    Peridynamics.update_sym_tensor!(M, 2, A)

    # a symmetric field has six rows, so `get_tensor` must not be used on it
    @test_throws BoundsError Peridynamics.get_tensor(M, 2)

    # the element type follows the field
    M32 = zeros(Float32, 6, 1)
    Peridynamics.update_sym_tensor!(M32, 1, SMatrix{3,3,Float32,9}(A))
    @test eltype(Peridynamics.get_sym_tensor(M32, 1)) === Float32
end

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

@testitem "invreg" setup=[Fixtures] begin
    using Peridynamics.StaticArrays
    using Peridynamics.LinearAlgebra
    rng = Fixtures.rng()

    # Test with well-conditioned matrices
    A = SMatrix{3,3}(rand(rng, 3, 3))
    @test Peridynamics.invreg(A, 0, 0) * A ≈ I
    @test Peridynamics.invreg(A, 0, sqrt(eps())) * A ≈ I
    @test Peridynamics.invreg(A, 1e-10, sqrt(eps())) * A ≈ I

    A = SMatrix{4,4}(rand(rng, 4, 4))
    @test Peridynamics.invreg(A, 0, 0) * A ≈ I
    @test Peridynamics.invreg(A, 0, sqrt(eps())) * A ≈ I
    @test Peridynamics.invreg(A, 1e-10, sqrt(eps())) * A ≈ I

    A = SMatrix{5,5}(rand(rng, 5, 5))
    @test Peridynamics.invreg(A, 0, 0) * A ≈ I
    @test Peridynamics.invreg(A, 0, sqrt(eps())) * A ≈ I
    @test Peridynamics.invreg(A, 1e-10, sqrt(eps())) * A ≈ I

    A = SMatrix{6,6}(rand(rng, 6, 6))
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

@testitem "invreg: a singular value floor damps only what is below it" setup=[Fixtures] begin
    using Peridynamics.StaticArrays
    using Peridynamics.LinearAlgebra
    rng = Fixtures.rng()

    # A well-conditioned matrix is inverted exactly: every singular value is above the floor,
    # so the damping is zero and there is no bias at all.
    for n in (3, 4, 5, 6)
        A = SMatrix{n,n}(rand(rng, n, n)) + n * I
        @test Peridynamics.invreg(A, 1e-3) * A ≈ I
        @test Peridynamics.invreg(A, 1e-3) ≈ inv(A)
        @test Peridynamics.invreg(A, 0) ≈ inv(A)
    end

    B = @SMatrix [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    @test Peridynamics.invreg(B, 1e-3) ≈ I
    @test Peridynamics.invreg(B, 0.5) ≈ I

    # Above the floor the result is bit-identical to the two-parameter method with `λ = 0`,
    # which is what keeps a floor from changing well-conditioned simulations.
    A = @SMatrix [1.0 0.5 0.3; 0.5 1.0 0.4; 0.3 0.4 1.0]
    @test Peridynamics.invreg(A, 1e-3) == Peridynamics.invreg(A, 0, 1e-3)

    # A diagonal matrix makes the per-singular-value behavior directly readable.
    # With a floor of σ_f = 1 (= 1e-1 * σ_max): σ = 10 and σ = 1 are above or at the floor and
    # are inverted exactly, σ = 0.5 is damped to (0.5)/(0.5² + 0.5²) = 1, and σ = 0 goes to 0.
    D = @SMatrix [10.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 0.5 0.0; 0.0 0.0 0.0 0.0]
    Dinv = Peridynamics.invreg(D, 1e-1)
    @test Dinv ≈ Diagonal(SVector(0.1, 1.0, 1.0, 0.0))

    # The inverse stays bounded no matter how singular the matrix gets. The damped branch
    # s / (s² + (σ_f - s)²) peaks at s = σ_f/√2, where it is 1 / ((2√2 - 2) σ_f).
    σmax, ε = 3.0, 1e-2
    bound = 1 / ((2 * sqrt(2) - 2) * ε * σmax)
    for σmin in (1e-1, 1e-3, ε * σmax / sqrt(2), 1e-6, 1e-12, 0.0)
        M = @SMatrix [σmax 0.0; 0.0 σmin]
        Minv = Peridynamics.invreg(M, ε)
        @test !Peridynamics.containsnan(Minv)
        @test opnorm(Array(Minv)) ≤ bound * (1 + 1e-12)
    end
    # and the bound is attained, so it is tight
    M = @SMatrix [σmax 0.0; 0.0 ε*σmax/sqrt(2)]
    @test opnorm(Array(Peridynamics.invreg(M, ε))) ≈ bound

    # The transition at the floor is C¹, so no force impulse when a singular value crosses it.
    σ_f = 1e-1
    f(s) = Peridynamics.invreg((@SMatrix [1.0 0.0; 0.0 s]), σ_f)[2, 2]
    @test f(σ_f) ≈ 1 / σ_f                                 # exact inversion at the floor
    @test f(σ_f / 2) ≈ 1 / σ_f                             # damped, still finite
    @test f(0.0) ≈ 0.0                                     # bounded in the singular limit
    h = 1e-6
    @test f(σ_f - h) ≈ f(σ_f + h) rtol=1e-4                # C⁰
    lhs = (f(σ_f) - f(σ_f - h)) / h
    rhs = (f(σ_f + h) - f(σ_f)) / h
    @test lhs ≈ rhs rtol=1e-4                              # C¹
    @test rhs ≈ -1 / σ_f^2 rtol=1e-4                       # and equal to d(1/s)/ds

    # The rank-deficient case that makes the plain inverse produce NaNs
    A = @SMatrix [1.0 1.0 1.0; 1.0 1.0+1e-10 1.0; 1.0 1.0 1.0+1e-10]
    @test Peridynamics.containsnan(inv(A))
    @test !Peridynamics.containsnan(Peridynamics.invreg(A, 1e-3))
    @test !Peridynamics.containsnan(Peridynamics.invreg(A, 0))

    # A zero matrix has no largest singular value to scale the floor with, so it must not
    # divide by zero either
    Z = zero(SMatrix{3,3,Float64,9})
    @test !Peridynamics.containsnan(Peridynamics.invreg(Z, 1e-3))
    @test Peridynamics.invreg(Z, 1e-3) ≈ Z
end
