@testsnippet DimMocks begin
    # the accessors of `static_arrays.jl` take the dimension as a `Val{N}`, and `dims` is
    # what produces it. These stubs answer `get_n_dim` the way a real system, a real storage
    # and a real model state do, through a type parameter that is never a literal `3`.
    struct Dim3System <: Peridynamics.AbstractSystem end
    struct Dim2System <: Peridynamics.AbstractSystem end
    Peridynamics.get_n_dim(::Dim3System) = 3
    Peridynamics.get_n_dim(::Dim2System) = 2

    struct Dim3Storage <: Peridynamics.AbstractStorage end
    struct Dim2Storage <: Peridynamics.AbstractStorage end
    Peridynamics.get_n_dim(::Dim3Storage) = 3
    Peridynamics.get_n_dim(::Dim2Storage) = 2
end

@testitem "dims: the Val of a system, a storage and a real chunk" setup=[DimMocks, Fixtures] begin
    using Peridynamics: dims, get_n_dim

    @test dims(Dim3System()) === Val(3)
    @test dims(Dim2System()) === Val(2)
    @test dims(Dim3Storage()) === Val(3)
    @test dims(Dim2Storage()) === Val(2)

    # a real chunk answers the same from its system and from its storage, and the storage of
    # a chunk is built with the dimension of its system, so the two can never disagree
    chunk = Fixtures.chunk(Fixtures.cube())
    @test dims(chunk.system) === Val(3)
    @test dims(chunk.storage) === Val(3)
    @test get_n_dim(chunk.storage) == get_n_dim(chunk.system)

    # a body knows its dimension too, but off the size of its position matrix, so it is a
    # runtime value and never a kernel argument
    @test dims(Fixtures.cube()) === Val(3)
end

@testitem "get_sym_tensor / update_sym_tensor!: Voigt round trip of a symmetric tensor" setup=[DimMocks] begin
    using Peridynamics.StaticArrays
    using Peridynamics.LinearAlgebra
    using Peridynamics: dims

    A = SMatrix{3,3,Float64,9}(1.0, 4.0, 5.0,
                               4.0, 2.0, 6.0,
                               5.0, 6.0, 3.0)

    M = zeros(6, 3)
    Peridynamics.update_sym_tensor!(M, 2, A, Val(3))

    # Voigt order: (11, 22, 33, 23, 13, 12)
    @test M[:, 2] == [1.0, 2.0, 3.0, 6.0, 5.0, 4.0]
    # the neighbouring columns must be untouched
    @test all(iszero, M[:, 1])
    @test all(iszero, M[:, 3])

    # the round trip is exact and gives back a symmetric tensor
    B = Peridynamics.get_sym_tensor(M, 2, Val(3))
    @test B == A
    @test B == B'

    # the same through `dims` of a system, of a storage and of a real storage
    @test Peridynamics.get_sym_tensor(M, 2, dims(Dim3System())) == A
    @test Peridynamics.get_sym_tensor(M, 2, dims(Dim3Storage())) == A

    # only the upper triangle is read and mirrored back, so an asymmetry of round-off size
    # survives the round trip as round-off and never leaks into another component
    Aperturbed = A + 1e-16 * SMatrix{3,3,Float64,9}(0.0, -1.0, 0.0, 1.0, 0.0, 0.0,
                                                    0.0, 0.0, 0.0)
    Peridynamics.update_sym_tensor!(M, 2, Aperturbed, Val(3))
    @test Peridynamics.get_sym_tensor(M, 2, Val(3)) ≈ A atol=1e-15

    # a tensor with a real skew part is not symmetrized: the upper triangle wins
    Askew = A + SMatrix{3,3,Float64,9}(0.0, -0.5, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0)
    Peridynamics.update_sym_tensor!(M, 2, Askew, Val(3))
    @test Peridynamics.get_sym_tensor(M, 2, Val(3))[1, 2] == Askew[1, 2]
    @test Peridynamics.get_sym_tensor(M, 2, Val(3))[2, 1] == Askew[1, 2]

    Peridynamics.update_sym_tensor!(M, 2, A, Val(3))

    # a symmetric field has six rows, so `get_tensor` must not be used on it
    @test_throws BoundsError Peridynamics.get_tensor(M, 2, Val(3))

    # the element type follows the field
    M32 = zeros(Float32, 6, 1)
    Peridynamics.update_sym_tensor!(M32, 1, SMatrix{3,3,Float32,9}(A), Val(3))
    @test eltype(Peridynamics.get_sym_tensor(M32, 1, Val(3))) === Float32

    # a write of the wrong static size is a `MethodError`, not a silent truncation
    @test_throws MethodError Peridynamics.update_sym_tensor!(M, 2,
                                                             SMatrix{2,2,Float64,4}(A[1:2,
                                                                                      1:2]),
                                                             Val(3))
end

@testitem "get_sym_tensor / update_sym_tensor!: N = 2 uses the (11, 22, 12) Voigt order" setup=[DimMocks] begin
    using Peridynamics.StaticArrays
    using Peridynamics: dims

    A = SMatrix{2,2,Float64,4}(1.0, 4.0,
                               4.0, 2.0)
    M = zeros(3, 1)
    Peridynamics.update_sym_tensor!(M, 1, A, Val(2))

    # Voigt order: (11, 22, 12)
    @test M[:, 1] == [1.0, 2.0, 4.0]

    B = Peridynamics.get_sym_tensor(M, 1, Val(2))
    @test B == A
    @test B == B'

    @test Peridynamics.get_sym_tensor(M, 1, dims(Dim2System())) == A
    @test Peridynamics.get_sym_tensor(M, 1, dims(Dim2Storage())) == A

    # the 3D write does not fit a 2D field either
    @test_throws MethodError Peridynamics.update_sym_tensor!(M, 1,
                                                             zero(SMatrix{3,3,Float64,9}),
                                                             Val(2))
end

@testitem "get and update tensors" setup=[DimMocks] begin
    using Peridynamics.StaticArrays
    using Peridynamics: dims

    A = zeros(9, 2)
    for i in eachindex(A)
        A[i] = i
    end
    B = Peridynamics.get_tensor(A, 1, Val(3))
    C = Peridynamics.get_tensor(A, 2, dims(Dim3System()))
    @test C == Peridynamics.get_tensor(A, 2, dims(Dim3Storage()))
    for i in eachindex(B, C)
        @test B[i] ≈ A[i, 1]
        @test C[i] ≈ A[i, 2]
    end

    Peridynamics.update_tensor!(A, 2, reverse(B), Val(3))
    Peridynamics.update_tensor!(A, 1, reverse(C), Val(3))
    @test A[:, 1] ≈ [18:-1:10;]
    @test A[:, 2] ≈ [9:-1:1;]

    before = A[:, 1]
    Peridynamics.update_add_tensor!(A, 1, Peridynamics.get_tensor(A, 2, Val(3)), Val(3))
    @test A[:, 1] ≈ before + A[:, 2]

    Peridynamics.zero_tensor!(A, 1, Val(3))
    @test all(iszero, A[:, 1])
    @test !all(iszero, A[:, 2])

    # a write of the wrong static size is a `MethodError`, not a silent truncation
    @test_throws MethodError Peridynamics.update_tensor!(A, 1,
                                                         zero(SMatrix{2,2,Float64,4}),
                                                         Val(3))
end

@testitem "get and update tensors: N = 2 uses a 2×2 static matrix" setup=[DimMocks] begin
    using Peridynamics.StaticArrays
    using Peridynamics: dims

    A = zeros(4, 1)
    for i in eachindex(A)
        A[i] = i
    end
    B = Peridynamics.get_tensor(A, 1, Val(2))
    @test B isa SMatrix{2,2,Float64,4}
    @test B[:] ≈ A[:, 1]
    @test Peridynamics.get_tensor(A, 1, dims(Dim2System())) == B
    @test Peridynamics.get_tensor(A, 1, dims(Dim2Storage())) == B

    Peridynamics.update_tensor!(A, 1, reverse(B), Val(2))
    @test A[:, 1] ≈ [4, 3, 2, 1]

    Peridynamics.zero_tensor!(A, 1, Val(2))
    @test all(iszero, A[:, 1])
end

@testitem "get and update vectors" setup=[DimMocks] begin
    using Peridynamics.StaticArrays
    using Peridynamics: dims

    a = zeros(3, 2)
    for i in eachindex(a)
        a[i] = i
    end
    @test Peridynamics.get_vector(a, 1, Val(3)) ≈ [1, 2, 3]
    @test Peridynamics.get_vector(a, 2, dims(Dim3System())) ≈ [4, 5, 6]
    @test Peridynamics.get_vector(a, 2, dims(Dim3Storage())) ≈ [4, 5, 6]

    Peridynamics.update_vector!(a, 1, SVector{3}(3.0, 2.0, 1.0), Val(3))
    @test a[:, 1] ≈ [3, 2, 1]

    Peridynamics.update_vector!(a, 2, SVector{3}(6.0, 5.0, 4.0), Val(3))
    @test a[:, 2] ≈ [6, 5, 4]

    Peridynamics.update_add_vector!(a, 1, SVector{3}(1.0, 1.0, 1.0), Val(3))
    Peridynamics.update_add_vector!(a, 2, SVector{3}(1.0, 1.0, 1.0), Val(3))
    @test a ≈ [4; 3; 2;; 7; 6; 5]

    @test Peridynamics.get_vector_diff(a, 1, 2, Val(3)) ≈ [3, 3, 3]
    @test Peridynamics.get_vector_diff(a, 1, 2, dims(Dim3System())) ≈ [3, 3, 3]

    # a write of the wrong static size is a `MethodError`, not a silent truncation
    @test_throws MethodError Peridynamics.update_vector!(a, 1, SVector{2}(1.0, 2.0), Val(3))
    @test_throws MethodError Peridynamics.update_add_vector!(a, 1, SVector{2}(1.0, 2.0),
                                                             Val(3))
end

@testitem "get and update vectors: N = 2 uses a 2-vector" setup=[DimMocks] begin
    using Peridynamics.StaticArrays
    using Peridynamics: dims

    a = zeros(2, 2)
    for i in eachindex(a)
        a[i] = i
    end
    @test Peridynamics.get_vector(a, 1, Val(2)) ≈ [1, 2]
    @test Peridynamics.get_vector(a, 2, dims(Dim2System())) ≈ [3, 4]
    @test Peridynamics.get_vector(a, 2, dims(Dim2Storage())) ≈ [3, 4]

    Peridynamics.update_vector!(a, 1, SVector{2}(2.0, 1.0), Val(2))
    @test a[:, 1] ≈ [2, 1]

    @test Peridynamics.get_vector_diff(a, 1, 2, Val(2)) ≈ [1, 3]

    @test_throws MethodError Peridynamics.update_vector!(a, 1, SVector{3}(1.0, 2.0, 3.0),
                                                         Val(2))
end

@testitem "the accessors exist only in the Val form" begin
    using Peridynamics.StaticArrays
    a = zeros(3, 2)
    A = zeros(9, 2)
    S = zeros(6, 2)
    Sys = Peridynamics.BondSystem{Peridynamics.NoCorrection,3}

    # no plain-matrix method, i.e. no silent 3D fallback
    @test !hasmethod(Peridynamics.get_vector, Tuple{typeof(a),Int})
    @test !hasmethod(Peridynamics.get_vector_diff, Tuple{typeof(a),Int,Int})
    @test !hasmethod(Peridynamics.update_vector!, Tuple{typeof(a),Int,SVector{3,Float64}})
    @test !hasmethod(Peridynamics.update_add_vector!,
                     Tuple{typeof(a),Int,SVector{3,Float64}})
    @test !hasmethod(Peridynamics.get_tensor, Tuple{typeof(A),Int})
    @test !hasmethod(Peridynamics.update_tensor!,
                     Tuple{typeof(A),Int,SMatrix{3,3,Float64,9}})
    @test !hasmethod(Peridynamics.update_add_tensor!,
                     Tuple{typeof(A),Int,SMatrix{3,3,Float64,9}})
    @test !hasmethod(Peridynamics.zero_tensor!, Tuple{typeof(A),Int})
    @test !hasmethod(Peridynamics.get_sym_tensor, Tuple{typeof(S),Int})
    @test !hasmethod(Peridynamics.update_sym_tensor!,
                     Tuple{typeof(S),Int,SMatrix{3,3,Float64,9}})

    # and no method that takes the system as the second argument either
    @test !hasmethod(Peridynamics.get_vector, Tuple{typeof(a),Sys,Int})
    @test !hasmethod(Peridynamics.get_vector_diff, Tuple{typeof(a),Sys,Int,Int})
    @test !hasmethod(Peridynamics.update_vector!,
                     Tuple{typeof(a),Sys,Int,SVector{3,Float64}})
    @test !hasmethod(Peridynamics.update_add_vector!,
                     Tuple{typeof(a),Sys,Int,SVector{3,Float64}})
    @test !hasmethod(Peridynamics.get_tensor, Tuple{typeof(A),Sys,Int})
    @test !hasmethod(Peridynamics.update_tensor!,
                     Tuple{typeof(A),Sys,Int,SMatrix{3,3,Float64,9}})
    @test !hasmethod(Peridynamics.update_add_tensor!,
                     Tuple{typeof(A),Sys,Int,SMatrix{3,3,Float64,9}})
    @test !hasmethod(Peridynamics.zero_tensor!, Tuple{typeof(A),Sys,Int})
    @test !hasmethod(Peridynamics.get_sym_tensor, Tuple{typeof(S),Sys,Int})
    @test !hasmethod(Peridynamics.update_sym_tensor!,
                     Tuple{typeof(S),Sys,Int,SMatrix{3,3,Float64,9}})
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
