# `CMaterial`: the classic correspondence formulation of `src/physics/correspondence.jl`.

@testmodule CorrespondenceCase begin
    using Peridynamics, Test

    # the origin, the three unit vectors and an isolated point at (2,2,2), unit volumes
    function tetra_with_isolated_point(mat; kwargs...)
        ref_position = [0.0 1.0 0.0 0.0 2.0
                        0.0 0.0 1.0 0.0 2.0
                        0.0 0.0 0.0 1.0 2.0]
        body = Body(mat, ref_position, fill(1.0, 5))
        material!(body; horizon=1.5, kwargs...)
        no_failure!(body)
        return body
    end

    # the force density of `body` after point 2 was moved by `Δx` in x within one step `Δt` at
    # time `t` (with the velocities the move implies, which the rotated formulations need for the
    # velocity gradient; `Δt = 0` leaves them zero)
    function force_after_shift(body; Δx=0.0015, t=0.0, Δt=0.0)
        dh = Peridynamics.threads_data_handler(body, VelocityVerlet(steps=1), 1)
        chunk = dh.chunks[1]
        @test chunk.storage.position == body.position
        @test iszero(chunk.storage.b_int)
        chunk.storage.position[1, 2] += Δx
        if Δt > 0
            chunk.storage.velocity[1, 2] = Δx / Δt
            chunk.storage.velocity_half[1, 2] = Δx / Δt
        end
        Peridynamics.calc_force_density!(chunk, t, Δt)
        return chunk.storage.b_int
    end
end

@testitem "CMaterial: force density of a displaced point" setup=[CorrespondenceCase] begin
    C = CorrespondenceCase
    body = C.tetra_with_isolated_point(CMaterial(model=NeoHookePenalty()); rho=1, E=1, nu=0.25, Gc=1.0)
    b_int = C.force_after_shift(body)
    # the isolated point has no bonds and feels nothing
    @test iszero(b_int[:, 5])
    # linear momentum is conserved (unit volumes)
    @test all(abs.(sum(b_int; dims=2)) .< 1e-12 * maximum(abs, b_int))
    # points 3 and 4 are symmetric with respect to the displaced bond 1-2: the same x-force,
    # and the force of 3 in y equals the force of 4 in z
    @test b_int[1, 3] ≈ b_int[1, 4] atol=1e-12 * maximum(abs, b_int)
    @test b_int[2, 3] ≈ b_int[3, 4]
    @test abs(b_int[3, 3]) < 1e-12 * maximum(abs, b_int)
    @test abs(b_int[2, 4]) < 1e-12 * maximum(abs, b_int)
    # the displaced point is pulled back, the origin pulled towards it
    @test b_int[1, 2] < 0 < b_int[1, 1]
    # regression values of the full force state, so that a change of the formulation shows
    @test b_int[:, 1] ≈ [0.007185432432164514, 0.0023974081843325646, 0.0023974081843347313]
    @test b_int[1, 2] ≈ -0.007185432432123352
    @test b_int[2, 3] ≈ -0.002397408184361381
end

@testitem "CMaterial: the zero-energy mode stabilizations agree on a tetrahedron" setup=[CorrespondenceCase] begin
    # Four points with complete families and one affine shift: the deformation is exactly
    # homogeneous within every family, so there is no zero-energy mode to stabilize and both
    # stabilizations give the force density of the bare formulation.
    C = CorrespondenceCase
    kwargs = (; rho=1, E=1, nu=0.25, Gc=1.0)
    b_silling = C.force_after_shift(C.tetra_with_isolated_point(CMaterial(model=NeoHookePenalty()); kwargs...))
    b_wan = C.force_after_shift(C.tetra_with_isolated_point(CMaterial(model=NeoHookePenalty(), zem=ZEMWan()); kwargs...))
    b_none = C.force_after_shift(C.tetra_with_isolated_point(CMaterial(model=NeoHookePenalty(), zem=ZEMSilling(Cs=0.0)); kwargs...))
    @test b_silling ≈ b_wan
    @test b_silling ≈ b_none
end

@testitem "CMaterial: one material on two point sets equals one material on all" setup=[CorrespondenceCase] begin
    # `material!` per point set with identical parameters must give the same forces as a single
    # `material!` call — the interface between the sets is invisible to the force density.
    C = CorrespondenceCase
    single = C.tetra_with_isolated_point(CMaterial(model=NeoHookePenalty()); rho=1, E=1, nu=0.25, Gc=1.0)
    ref_position = single.position
    split = Body(CMaterial(model=NeoHookePenalty()), ref_position, fill(1.0, 5))
    point_set!(split, :a, [1, 2, 3])
    point_set!(split, :b, [4, 5])
    material!(split, :a; horizon=1.5, rho=1, E=1, nu=0.25, Gc=1.0)
    material!(split, :b; horizon=1.5, rho=1, E=1, nu=0.25, Gc=1.0)
    no_failure!(split)
    @test C.force_after_shift(split) ≈ C.force_after_shift(single)
end

@testitem "CMaterial: exported von Mises stress" begin
    using Peridynamics.StaticArrays
    body = Body(CMaterial(), [0.0 1.0; 0.0 0.0; 0.0 0.0], [1.0, 1.0])
    material!(body, horizon=1.5, rho=8000, E=210e9, nu=0.25)
    dh = Peridynamics.threads_data_handler(body, VelocityVerlet(steps=1), 1)
    (; mat, storage, system, paramsetup) = dh.chunks[1]
    # σ = 100 everywhere: the deviator has no diagonal and the von Mises stress is √3 ⋅ 100 ⋅ √3
    σ = @SMatrix fill(100.0, 3, 3)
    for i in 1:2
        Peridynamics.update_tensor!(storage.cauchy_stress, i, σ)
    end
    σvm = Peridynamics.export_field(Val(:von_mises_stress), mat, system, storage, paramsetup, 0.0)
    @test all(σvm .≈ 300.0)
end

@testitem "CMaterial: show" begin
    @test contains(sprint(show, CMaterial()), "CMaterial")
    @test contains(sprint(show, CMaterial(maxdmg=0.5)), "maxdmg=0.5")
    @test contains(sprint(show, MIME("text/plain"), CMaterial()), "CMaterial")
end
