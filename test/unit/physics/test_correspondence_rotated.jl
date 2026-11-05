@testitem "CRMaterial" begin
    mat = CRMaterial()
    @test mat.constitutive_model isa SaintVenantKirchhoff

    @test_throws ArgumentError CRMaterial(model=NeoHookePenalty())

    io = IOBuffer()
    show(IOContext(io, :compact=>true), MIME("text/plain"), mat)
    msg = String(take!(io))
    @test startswith(msg, "CRMaterial")
end

@testitem "CRMaterial: force density of a displaced point" setup=[CorrespondenceCase] begin
    # The rotated formulation integrates the unrotated stress rate from the velocity gradient
    # over the step (so the shift needs the velocities it implies). For a first step from the
    # reference configuration the rotation is the identity, and with the linear elastic model
    # the force density equals that of the classic correspondence formulation with the same
    # model up to the linearization in the strain: 0.25 % for a strain of 0.1 %, 20 % for 10 %.
    C = CorrespondenceCase
    kwargs = (; rho=8000, E=210e9, nu=0.25, Gc=100.0)
    Δx, Δt = 0.001, 1e-7
    b_cr = C.force_after_shift(C.tetra_with_isolated_point(CRMaterial(model=LinearElastic()); kwargs...);
                               Δx, t=Δt, Δt)
    b_c = C.force_after_shift(C.tetra_with_isolated_point(CMaterial(model=LinearElastic()); kwargs...);
                              Δx, t=Δt, Δt)
    @test iszero(b_cr[:, 5]) # the isolated point
    @test all(abs.(sum(b_cr; dims=2)) .< 1e-12 * maximum(abs, b_cr)) # linear momentum
    @test b_cr[1, 2] < 0 < b_cr[1, 1] # the displaced point is pulled back
    @test b_cr ≈ b_c rtol=1e-2
    @test !(b_cr ≈ b_c) # but not exactly: the rate form is linear in the strain
end

@testitem "CRMaterial: the zero-energy mode stabilizations agree on a tetrahedron" setup=[CorrespondenceCase] begin
    # no zero-energy mode exists for four points under an affine shift, see the CMaterial item
    C = CorrespondenceCase
    kwargs = (; rho=8000, E=210e9, nu=0.25, Gc=100.0)
    Δx, Δt = 0.1, 1e-6
    b_silling = C.force_after_shift(C.tetra_with_isolated_point(CRMaterial(model=LinearElastic()); kwargs...);
                                    Δx, t=Δt, Δt)
    b_wan = C.force_after_shift(C.tetra_with_isolated_point(CRMaterial(model=LinearElastic(), zem=ZEMWan()); kwargs...);
                                Δx, t=Δt, Δt)
    @test b_silling ≈ b_wan
end

@testitem "CRMaterial: exported von Mises stress from the unrotated stress" begin
    using Peridynamics.StaticArrays, Peridynamics.LinearAlgebra
    body = Body(CRMaterial(), [0.0 1.0; 0.0 0.0; 0.0 0.0], [1.0, 1.0])
    material!(body, horizon=1.5, rho=8000, E=210e9, nu=0.25)
    dh = Peridynamics.threads_data_handler(body, VelocityVerlet(steps=1), 1)
    (; mat, storage, system, paramsetup) = dh.chunks[1]
    # the rotation is the identity at initialization
    @test all(Peridynamics.get_tensor(storage.rotation, i) ≈ I for i in 1:2)
    σ = @SMatrix fill(100.0, 3, 3)
    for i in 1:2
        Peridynamics.update_tensor!(storage.unrotated_stress, i, σ)
    end
    σvm = Peridynamics.export_field(Val(:von_mises_stress), mat, system, storage, paramsetup, 0.0)
    @test all(σvm .≈ 300.0)
end
